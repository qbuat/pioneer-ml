import awkward as ak
import numpy as np


def get_pion_rz(arr):
    _decay_x = ak.firsts(arr['decay.posX'])
    _decay_y = ak.firsts(arr['decay.posY'])
    _decay_r = np.sqrt(_decay_x * _decay_x + _decay_y * _decay_y)
    _decay_z = ak.firsts(arr['decay.posZ'])
    return _decay_r, _decay_z
