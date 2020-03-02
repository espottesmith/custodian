# coding: utf-8

from __future__ import unicode_literals, division
import os
import shutil
import subprocess

from custodian.custodian import Job


__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "02/17/20"


class GSMJob(Job):
    """
    A basic job using the Growing String Method as implemented in pyGSM.
    """

    def __init__(self,
                 command,
                 package,
                 charge,
                 xyz_file="input_geometry.xyz",
                 input_file="gsm.inp",
                 output_file="gsm.out",
                 mode="DE_GSM",
                 suffix="",
                 string_id=0,
                 calc_loc=None,
                 save_scratch=False,
                 multimode=True,
                 max_cores=32,
                 ends_fixed=False,
                 num_nodes=9,
                 additional_options=None):
        """
        Args:
            command (str): Command to run pyGSM.
            package (str): Electronic structure package to use with pyGSM.
                Technically, a wide range of packages can be used, but within
                custodian, only "QChem" and "XTB" are supported.
            xyz_file (str): Name of the input geometry file for GSM.
            input_file (str): Name of the electronic structure package input file
                for the calculation (format depends on the package).
            output_file (str): Name of the file where all output from pyGSM
                should be piped.
            mode (str): One of DE_GSM (for the double-ended growing-string
                method), SE_GSM (for the single-ended growing-string method), or
                SE_Cross (for the crossing single-ended growing-string method).
                Default is "DE_GSM".
            suffix (str): String to append to the file in postprocess.
            string_id (int): ID for scratch file directory
            calc_loc (str): Path where pyGSM should run. Defaults to None, in
                which case the program will run in the system-defined SCRATCH.
            save_scratch (bool): Whether to save full scratch directory contents.
                Defaults to False.
            multimode (bool): If True (default), then use multiprocessing.
            max_cores (int): Number of processes to use for parallelization
            ends_fixed (bool): If True (default False), do not optimize the ends
                of the string (the reactants and products).
            additional_options (dict, or None): If not None (default), this
                contains key-value pairs for all additional command-line options
                for pyGSM.
        """
        self.command = command.split(" ")
        if package in ["QChem", "XTB"]:
            self.package = package
        else:
            raise ValueError("Invalid electronic structure package given. "
                             "Options include QChem and XTB.")
        if mode in ["DE_GSM", "SE_GSM", "SE_Cross"]:
            self.mode = mode
        else:
            raise ValueError("Invalid GSM mode given. "
                             "Options include DE_GSM, SE_GSM, and SE_Cross.")
        self.charge = charge
        self.xyz_file = xyz_file
        self.input_file = input_file
        self.output_file = output_file
        self.suffix = suffix
        self.string_id = string_id
        self.calc_loc = calc_loc
        self.save_scratch = save_scratch
        self.multimode = multimode
        self.max_cores = max_cores
        self.ends_fixed = ends_fixed
        self.num_nodes = num_nodes

        if additional_options is None:
            self.additional_options = dict()
        else:
            self.additional_options = additional_options

    @property
    def current_command(self):
        command = self.command + ["-mode", self.mode, "-xyzfile", self.xyz_file,
                                  "-lot_inp_file", self.input_file, "-package",
                                  self.package, "-ID", str(self.string_id),
                                  "-num_nodes", str(self.num_nodes), "-charge",
                                  str(self.charge)]

        if self.ends_fixed:
            if self.mode == "DE_GSM":
                com = ["-reactant_geom_fixed", "-product_geom_fixed"]
            else:
                com = ["-reactant_geom_fixed"]
            command += com

        if self.multimode:
            com = ["-use_multiprocessing", "-nproc", str(self.max_cores)]
            command += com

        for key, value in self.additional_options.items():
            command.append("-" + key)
            if value is not None:
                command.append(value)

        com_str = " ".join(command)

        return com_str

    def setup(self):
        if self.multimode:
            os.environ['OMP_NUM_THREADS'] = str(self.max_cores)
        os.environ["SCRATCH"] = os.getcwd()
        print("Current SCRATCH is: {}".format(os.environ["SCRATCH"]))
        if self.calc_loc is not None:
            os.environ["SCRATCH"] = self.calc_loc

    def postprocess(self):
        scratch_dir = os.path.join(os.environ["SCRATCH"], "scratch")
        string_dir = os.path.join(os.environ["SCRATCH"],
                                  "string_{:03d}".format(self.string_id))
        if self.suffix != "":
            shutil.move(self.input_file, self.input_file + self.suffix)
            shutil.move(self.output_file, self.output_file + self.suffix)
        if not self.save_scratch:
            shutil.rmtree(scratch_dir, ignore_errors=True)
            shutil.rmtree(string_dir, ignore_errors=True)

    def run(self):
        """
        Perform the actual GSM run.

        Returns:
            (subprocess.Popen) Used for monitoring.
        """
        local_scratch = os.path.join(os.environ["SCRATCH"], "scratch")
        local_string = os.path.join(os.environ["SCRATCH"], "string_{:03d}".format(self.string_id))
        if os.path.exists(local_scratch):
            shutil.rmtree(local_scratch)
            shutil.rmtree(local_string)
        p = subprocess.Popen(self.current_command,
                             stdout=open(self.output_file, 'w'),
                             shell=True)
        return p
