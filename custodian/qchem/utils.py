import math

import numpy as np

from pymatgen.io.qchem.outputs import QCOutput


def perturb_coordinates(old_coords, negative_freq_vecs, molecule_perturb_scale,
                        reversed_direction=False):
    max_dis = max(
        [math.sqrt(sum([x ** 2 for x in vec])) for vec in negative_freq_vecs])
    scale = molecule_perturb_scale / max_dis
    normalized_vecs = [[x * scale for x in vec] for vec in negative_freq_vecs]
    direction = 1.0
    if reversed_direction:
        direction = -1.0
    return [[c + v * direction for c, v in zip(coord, vec)]
            for coord, vec in zip(old_coords, normalized_vecs)]


def vector_list_diff(vecs1, vecs2):
    diff = 0.0
    if len(vecs1) != len(vecs2):
        raise AssertionError("ERROR: Vectors must be of equal length! Exiting...")
    for ii, vec1 in enumerate(vecs1):
        diff += np.linalg.norm(vecs2[ii] - vec1)
    return diff
