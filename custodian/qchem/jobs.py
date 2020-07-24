# coding: utf-8

from __future__ import unicode_literals, division
import math

# New QChem job module

import os
import shutil
import copy
import subprocess
import numpy as np
from pymatgen.core import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.fragmenter import metal_edge_extender
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import (QCOutput,
                                       ScratchFileParser,
                                       check_for_structure_changes)

from custodian.custodian import Job
from custodian.qchem.utils import perturb_coordinates, vector_list_diff


__author__ = "Samuel Blau, Brandon Wood, Shyam Dwaraknath, Evan Spotte-Smith"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "3/20/18"
__credits__ = "Xiaohui Qu"


class QCJob(Job):
    """
    A basic QChem Job.
    """

    def __init__(self,
                 qchem_command,
                 max_cores,
                 multimode="openmp",
                 input_file="mol.qin",
                 output_file="mol.qout",
                 qclog_file="mol.qclog",
                 suffix="",
                 calc_loc=None,
                 save_scratch=False,
                 backup=True):
        """
        Args:
            qchem_command (str): Command to run QChem.
            max_cores (int): Maximum number of cores to parallelize over.
            multimode (str): Parallelization scheme, either openmp or mpi.
            input_file (str): Name of the QChem input file.
            output_file (str): Name of the QChem output file.
            qclog_file (str): Name of the file to redirect the standard output
                to. None means not to record the standard output.
            suffix (str): String to append to the file in postprocess.
            calc_loc (str): Path where Q-Chem should run. Defaults to None, in
                which case Q-Chem will run in the system-defined QCLOCALSCR.
            save_scratch (bool): Whether to save full scratch directory contents.
                Defaults to False.
            backup (bool): Whether to backup the initial input file. If True, the
                input will be copied with a ".orig" appended. Defaults to True.
        """
        self.qchem_command = qchem_command.split(" ")
        self.multimode = multimode
        self.input_file = input_file
        self.output_file = output_file
        self.max_cores = max_cores
        self.qclog_file = qclog_file
        self.suffix = suffix
        self.calc_loc = calc_loc
        self.save_scratch = save_scratch
        self.backup = backup

    @property
    def current_command(self):
        multi = {"openmp": "-nt", "mpi": "-np"}
        if self.multimode not in multi:
            raise RuntimeError("ERROR: Multimode should only be set to openmp or mpi")
        command = [multi[self.multimode], str(self.max_cores), self.input_file, self.output_file, "scratch"]
        command = self.qchem_command + command
        com_str = ""
        for part in command:
            com_str = com_str + " " + part
        return com_str

    def setup(self):
        if self.backup:
            shutil.copy(self.input_file, "{}.orig".format(self.input_file))
        if self.multimode == 'openmp':
            os.environ['QCTHREADS'] = str(self.max_cores)
            os.environ['OMP_NUM_THREADS'] = str(self.max_cores)
        os.environ["QCSCRATCH"] = os.getcwd()
        print("Current QCSCRATCH is: {}".format(os.environ["QCSCRATCH"]))
        if self.calc_loc is not None:
            os.environ["QCLOCALSCR"] = self.calc_loc

    def postprocess(self):
        scratch_dir = os.path.join(os.environ["QCSCRATCH"], "scratch")
        for file in ["HESS", "GRAD", "plots/dens.0.cube"]:
            file_path = os.path.join(scratch_dir, file)
            if os.path.exists(file_path):
                shutil.copy(file_path, os.getcwd())
        if self.suffix != "":
            shutil.move(self.input_file, self.input_file + self.suffix)
            shutil.move(self.output_file, self.output_file + self.suffix)
            shutil.move(self.qclog_file, self.qclog_file + self.suffix)
            for file in ["HESS", "GRAD", "dens.0.cube"]:
                if os.path.exists(file):
                    shutil.move(file, file + self.suffix)
        if not self.save_scratch:
            shutil.rmtree(scratch_dir)

    def run(self):
        """
        Perform the actual QChem run.

        Returns:
            (subprocess.Popen) Used for monitoring.
        """
        local_scratch = os.path.join(os.environ["QCLOCALSCR"], "scratch")
        if os.path.exists(local_scratch):
            shutil.rmtree(local_scratch)
        qclog = open(self.qclog_file, 'w')
        p = subprocess.Popen(self.current_command, stdout=qclog, shell=True)
        return p

    @classmethod
    def opt_with_frequency_flattener(cls,
                                     qchem_command,
                                     multimode="openmp",
                                     input_file="mol.qin",
                                     output_file="mol.qout",
                                     qclog_file="mol.qclog",
                                     max_iterations=10,
                                     max_molecule_perturb_scale=0.3,
                                     check_connectivity=True,
                                     linked=True,
                                     transition_state=False,
                                     freq_before_opt=False,
                                     save_final_scratch=False,
                                     **QCJob_kwargs):
        """
        Optimize a structure and calculate vibrational frequencies to check if the
        structure is in a true minima. If a frequency is negative, iteratively
        perturb the geometry, optimize, and recalculate frequencies until all are
        positive, aka a true minima has been found.

        Args:
            qchem_command (str): Command to run QChem.
            multimode (str): Parallelization scheme, either openmp or mpi.
            input_file (str): Name of the QChem input file.
            output_file (str): Name of the QChem output file.
            max_iterations (int): Number of perturbation -> optimization -> frequency
                iterations to perform. Defaults to 10.
            linked (bool): Whether or not to use the linked flattener. Defaults to True.
            save_final_scratch (bool): Whether to save full scratch directory contents
                at the end of the flattening. Defaults to False.
            transition_state (bool): If True (default False), use a ts
                optimization (search for a saddle point instead of a minimum)
            freq_before_opt (bool): If True (default False), run a frequency
                calculation before any opt/ts searches to improve understanding
                of the local potential energy surface.
            **QCJob_kwargs: Passthrough kwargs to QCJob. See
                :class:`custodian.qchem.jobs.QCJob`.
        """
        if not os.path.exists(input_file):
            raise AssertionError('Input file must be present!')

        if transition_state:
            opt_method = "ts"
        else:
            opt_method = "opt"

        energy_diff_cutoff = 0.0000001

        orig_input = QCInput.from_file(input_file)
        freq_rem = copy.deepcopy(orig_input.rem)
        freq_rem["job_type"] = "freq"
        opt_rem = copy.deepcopy(orig_input.rem)
        opt_rem["scf_guess_always"] = True
        if "geom_opt_max_cycles" not in opt_rem:
            opt_rem["geom_opt_max_cycles"] = 250
        opt_rem["job_type"] = opt_method
        opt_rem["geom_opt_hessian"] = "read"
        first = True
        energy_history = list()

        if freq_before_opt:
            yield (QCJob(qchem_command=qchem_command,
                         multimode=multimode,
                         input_file=input_file,
                         output_file=output_file,
                         qclog_file=qclog_file,
                         suffix=".freq_pre",
                         save_scratch=True,
                         backup=first,
                         **QCJob_kwargs))

            opt_QCInput = QCInput(molecule=orig_input.molecule,
                                  rem=opt_rem,
                                  opt=orig_input.opt,
                                  pcm=orig_input.pcm,
                                  solvent=orig_input.solvent,
                                  smx=orig_input.smx)
            opt_QCInput.write_file(input_file)
            first = False

        for ii in range(max_iterations):
            yield (QCJob(
                qchem_command=qchem_command,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                qclog_file=qclog_file,
                suffix=".{}_".format(opt_method) + str(ii),
                save_scratch=True,
                backup=first,
                **QCJob_kwargs))
            opt_outdata = QCOutput(output_file + ".{}_".format(opt_method) + str(ii)).data
            opt_indata = QCInput.from_file(input_file + ".{}_".format(opt_method) + str(ii))
            if opt_indata.rem["scf_algorithm"] != freq_rem["scf_algorithm"]:
                freq_rem["scf_algorithm"] = opt_indata.rem["scf_algorithm"]
                opt_rem["scf_algorithm"] = opt_indata.rem["scf_algorithm"]
            first = False
            if opt_outdata["structure_change"] == "unconnected_fragments" and not opt_outdata["completion"]:
                # TODO: Do we want this behavior for TS searching?
                break
            else:
                energy_history.append(opt_outdata.get("final_energy"))
                freq_QCInput = QCInput(
                    molecule=opt_outdata.get("molecule_from_optimized_geometry"),
                    rem=freq_rem,
                    opt=orig_input.opt,
                    pcm=orig_input.pcm,
                    solvent=orig_input.solvent,
                    smx=orig_input.smx)
                freq_QCInput.write_file(input_file)
                yield (QCJob(
                    qchem_command=qchem_command,
                    multimode=multimode,
                    input_file=input_file,
                    output_file=output_file,
                    qclog_file=qclog_file,
                    suffix=".freq_" + str(ii),
                    save_scratch=True,
                    backup=first,
                    **QCJob_kwargs))
                outdata = QCOutput(output_file + ".freq_" + str(ii)).data
                indata = QCInput.from_file(input_file + ".freq_" + str(ii))
                if indata.rem["scf_algorithm"] != freq_rem["scf_algorithm"]:
                    freq_rem["scf_algorithm"] = indata.rem["scf_algorithm"]
                    opt_rem["scf_algorithm"] = indata.rem["scf_algorithm"]
                errors = outdata.get("errors")
                if len(errors) != 0:
                    raise AssertionError('No errors should be encountered while flattening frequencies!')
                if transition_state:
                    if outdata.get('frequencies')[0] < 0.0 and outdata.get('frequencies')[1] > 0.0:
                        print("Saddle point found!")
                        break
                    elif abs(outdata.get("frequencies")[1]) < 15.0 and outdata.get("frequencies")[2] > 0.0:
                        print("Second small imaginary frequency (smaller than 15.0) - not worth further flattening!")
                        break
                    else:
                        opt_QCInput = QCInput(
                            molecule=opt_outdata.get("molecule_from_optimized_geometry"),
                            rem=opt_rem,
                            opt=orig_input.opt,
                            pcm=orig_input.pcm,
                            solvent=orig_input.solvent,
                            smx=orig_input.smx)
                        opt_QCInput.write_file(input_file)
                else:
                    if outdata.get('frequencies')[0] > 0.0:
                        print("All frequencies positive!")
                        break
                    elif abs(outdata.get('frequencies')[0]) < 15.0 and outdata.get('frequencies')[1] > 0.0:
                        print("One negative frequency smaller than 15.0 - not worth further flattening!")
                        break
                    else:
                        if len(energy_history) > 1:
                            if abs(energy_history[-1] - energy_history[-2]) < energy_diff_cutoff:
                                print("Energy change below cutoff!")
                                break
                        opt_QCInput = QCInput(
                            molecule=opt_outdata.get("molecule_from_optimized_geometry"),
                            rem=opt_rem,
                            opt=orig_input.opt,
                            pcm=orig_input.pcm,
                            solvent=orig_input.solvent,
                            smx=orig_input.smx)
                        opt_QCInput.write_file(input_file)
        if not save_final_scratch:
            shutil.rmtree(os.path.join(os.getcwd(), "scratch"))
