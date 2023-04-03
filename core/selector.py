import os
import awkward as ak
import numpy as np
import tabulate

from . import log; log = log.getChild(__name__)

_pion_mass = 139.570 * 1e3 # keV
_muon_mass = 105.658 * 1e3 # keV


def select_atar_dar(arr, evt_type='pidar_mudar', decay_in_atar=True, verbose=False):
    """
    """
    cutflow = [('init', len(arr))]
    if evt_type in ('pidar_mudar', 'pidar_mudif', 'pidif_mudar',):
        decay_volume = 2
    elif evt_type == 'pienu':
        decay_volume = 1
    else:
        raise ValueError
    
    
    # ---
    arr = arr[ak.num(arr['decay.volume']) == decay_volume]
    cutflow += [('decays-recorded', len(arr))]

    if evt_type == 'pidar_mudif':
        # arr = arr[np.fabs(1e3 * arr['decay.motherEnergy'][:,1] - _muon_mass) < 2]
        arr = arr[np.fabs(1e3 * arr['decay.motherEnergy'][:,1]) > 1]
        cutflow += [('muon-decay-in-flight', len(arr))]

    # ---
    # arr = arr[np.fabs(1e3 * arr['decay.motherEnergy'][:,0] - _pion_mass) < 2]
    if evt_type in ('pidar_mudar', 'pidar_mudif', 'pienu'):
        arr = arr[np.fabs(1e3 * arr['decay.motherEnergy'][:,0]) < 2]
        cutflow += [('pion-decay-at-rest', len(arr))]

    # ---
    if evt_type in ('pidar_mudar', 'pidif_mudar',):
        # arr = arr[np.fabs(1e3 * arr['decay.motherEnergy'][:,1] - _muon_mass) < 2]
        arr = arr[np.fabs(1e3 * arr['decay.motherEnergy'][:,1]) < 2]
        cutflow += [('muon-decay-at-rest', len(arr))]


    # # ---
    if decay_in_atar:
        arr = arr[np.fabs(arr['decay.volume'][:,0] - 100000) < 50000]
        cutflow += [('pion-decay-in-atar', len(arr))]

    # ---
    if evt_type in ('pidar_mudar', 'pidar_mudif', 'pidif_mudar',):
        if decay_in_atar:
            arr = arr[np.fabs(arr['decay.volume'][:,1] - 100000) < 50000]
            cutflow += [('muon-decay-in-atar', len(arr))]

    


        # ---
    arr = arr[np.fabs(ak.max(arr['atar.time'], axis=1) - ak.min(arr['atar.time'], axis=1)) > 3]
    cutflow += [('atar hits: |first - last| > 3 ns', len(arr))]

    if verbose:
        print(tabulate.tabulate(cutflow, headers=['cut', 'yields']))

    return arr



def select_pimunu_atar_pidar(arr, mu_dar=True, verbose=False):
    """
    """
    cutflow = [('init', len(arr))]
    decay_volume = 2
    
    # ---
    arr = arr[ak.num(arr['decay.volume']) == decay_volume]
    cutflow += [('decays-recorded', len(arr))]

    # ---
    arr = arr[np.fabs(arr['decay.volume'][:,0] - 100000) < 50000]
    cutflow += [('pion-decay-in-atar', len(arr))]

    # ---
    arr = arr[np.fabs(arr['decay.volume'][:,1] - 100000) < 50000]
    cutflow += [('muon-decay-in-atar', len(arr))]

    # ---
    # arr = arr[np.fabs(1e3 * arr['decay.motherEnergy'][:,0] - _pion_mass) < 2]
    arr = arr[np.fabs(1e3 * arr['decay.motherEnergy'][:,0]) < 2]
    cutflow += [('pion-decay-at-rest', len(arr))]
    
    if mu_dar:
        arr = arr[np.fabs(1e3 * arr['decay.motherEnergy'][:,1]) < 2]
        cutflow += [('muon-decay-at-rest', len(arr))]


    if verbose:
        print(tabulate.tabulate(cutflow, headers=['cut', 'yields']))

    return arr
