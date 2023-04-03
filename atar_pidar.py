#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import uproot
import numpy as np
import awkward as ak
import particle
from matplotlib import pyplot as plt
import argparse
import logging
import tabulate
import ROOT
import array
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

from atar.utils import geometry, compute_atar_ids, energy_centrality
_atar_geo = geometry('/Users/qbuat/pioneer/main/MonteCarlo/geometry/generator/test_output.gdml')
# _atar_geo = geometry('/Users/quentin/pioneer/PIONEER/MonteCarlo/geometry/generator/test_output.gdml')
def get_atar_pos(pix, variable='x'):
    if pix == -999:
        return -999
    return _atar_geo[pix][variable]
v_get_atar_pos = np.vectorize(get_atar_pos)

FIELDS = [
    'atar.pdgid',  'atar.pxlID', 'atar.edep', 'atar.time',
    # 'calo.Etot', 'calo.theta',
    'decay.motherPDGID', 'decay.time', 'decay.motherEnergy', 'decay.volume',
    'decay.posX', 'decay.posY', 'decay.posZ',
    'decay.daughterPDGID', 'decay.daughterMomX', 'decay.daughterMomY', 'decay.daughterMomZ',
    'init.mom_x', 'init.mom_y',
    ]


def fit():
    pass

    # need to implement pxlID uniqueness check
def reco_data(arr):
    # pion = arr['atar.pxlID'][arr['atar.pdgid'] == 211]
    pxl = arr['atar.pxlID']
    pxl_x = pxl[arr['atar.planeid']%2 == 0]
    pxl_y = pxl[arr['atar.planeid']%2 != 0]
    # x = v_get_atar_pos(arr['atar.pxlID'], 'x')
    # y = v_get_atar_pos(arr['atar.pxlID'], 'y')
    # z = v_get_atar_pos(arr['atar.pxlID'], 'z')
    
    x = v_get_atar_pos(pxl_x, 'x')
    y = v_get_atar_pos(pxl_y, 'y')
    z_x = v_get_atar_pos(pxl_x, 'z')
    z_y = v_get_atar_pos(pxl_y, 'z')
    e_x = []
    e_y = []
    e_z_x = np.array([0.006 for i in range(len(z_x))])
    e_z_y = np.array([0.006 for i in range(len(z_y))])
    e_x = np.array([0.01 for i in range(len(x))])
    e_y = np.array([0.01 for i in range(len(y))])
    return x, y, z_x, z_y, e_x, e_y, e_z_x, e_z_y
    

def plot_atar_energy(pienu, pimunu, plot_suffix):
    """
    """


    _anti_electron_pimunu_px =  pimunu['decay.daughterMomX'][ pimunu['decay.daughterPDGID']==-11]
    _anti_electron_pimunu_py =  pimunu['decay.daughterMomY'][ pimunu['decay.daughterPDGID']==-11]
    _anti_electron_pimunu_pz =  pimunu['decay.daughterMomZ'][ pimunu['decay.daughterPDGID']==-11] 
    _anti_electron_pimunu_p = np.sqrt(_anti_electron_pimunu_px*_anti_electron_pimunu_px + _anti_electron_pimunu_py*_anti_electron_pimunu_py +_anti_electron_pimunu_pz*_anti_electron_pimunu_pz)

    _anti_electron_pienu_px =  pienu['decay.daughterMomX'][ pienu['decay.daughterPDGID']==-11]
    _anti_electron_pienu_py =  pienu['decay.daughterMomY'][ pienu['decay.daughterPDGID']==-11]
    _anti_electron_pienu_pz =  pienu['decay.daughterMomZ'][ pienu['decay.daughterPDGID']==-11]
    _anti_electron_pienu_p = np.sqrt(_anti_electron_pienu_px*_anti_electron_pienu_px + _anti_electron_pienu_py*_anti_electron_pienu_py +_anti_electron_pienu_pz*_anti_electron_pienu_pz)

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Anti-electron True Momentum [MeV]')
    ax.set_ylabel('Entries')
    plt.hist(ak.flatten(_anti_electron_pienu_p, axis=None), bins=100, range=(0, 100), color='blue', label=r'$\pi^{+}\to e^{+}\nu_{e}$')
    plt.hist(ak.flatten(_anti_electron_pimunu_p, axis=None), bins=100, range=(0, 100), color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
    plt.yscale('log')
    plt.legend()
    plt.savefig(os.path.join(
        'plots', 'antielectron_momentum_{}.pdf'.format(plot_suffix)))


    
    _pimunu_last_atar_hit_time = ak.max(pimunu['atar.time'], axis=1)
    _pienu_last_atar_hit_time = ak.max(pienu['atar.time'], axis=1)

    # remove atar hits within 1ns of the last one
    _pimunu_atar_energy_wo_pos = pimunu['atar.edep'][np.fabs(pimunu['atar.time'] - _pimunu_last_atar_hit_time) > 1]
    _pienu_atar_energy_wo_pos = pienu['atar.edep'][np.fabs(pienu['atar.time'] - _pienu_last_atar_hit_time) > 1]

    _pimunu_energy_sum = ak.sum(_pimunu_atar_energy_wo_pos, axis=1)
    _pienu_energy_sum = ak.sum(_pienu_atar_energy_wo_pos, axis=1)

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('ATAR Energy [MeV]')
    ax.set_ylabel('Entries')
    plt.hist(_pimunu_energy_sum, bins=100, range=(0, 20), color='blue', label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
    plt.hist(_pienu_energy_sum, bins=100, range=(0, 20), color='red', alpha=0.5,  label=r'$\pi^{+}\to e^{+}\nu_{e}$')
    plt.yscale('log')
    plt.legend()
    plt.savefig(os.path.join(
        'plots', 'atar_energy_{}.pdf'.format(plot_suffix)))

    if  len(_pienu_energy_sum) != 0:
        efficiency = len(_pienu_energy_sum[_pienu_energy_sum<14]) / len(_pienu_energy_sum)
    else:
        efficiency = -1
    if len(_pimunu_energy_sum[_pimunu_energy_sum<14]) != 0:
        rejection = len(_pimunu_energy_sum) / len(_pimunu_energy_sum[_pimunu_energy_sum<14])
    else:
        rejection = -1
    return efficiency, rejection
    
_pion_mass = 139.570 * 1e3
_muon_mass = 105.658 * 1e3

if __name__ == '__main__':


    do_plots = True
    

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    log.addHandler(ch)
    log.info('Starting basic analysis!')


    # _data_file = './data_v3/'
    # _data_file = '/Users/quentin/Desktop/for_vincent/'
    # _data_file = '/Users/qbuat/Downloads/noDegrader/'
    # _data_file = '/Users/qbuat/pioneer/main/data_55_2.5/'

    # _data_file = '/Users/qbuat/pioneer/data_60_2.5'
    # _data_file = '/Users/qbuat/pioneer/data_65_2.5_degrader0.5'
    _data_file = '/Users/qbuat/pioneer/data_55_2.5_nodegrader'
    _tree_name = 'sim'
    _pienu_files = []
    for _f in os.listdir(os.path.join(
            _data_file, 'pienu')):
        _file_abspath = os.path.join(
            _data_file, 'pienu', _f)
        _pienu_files += [_file_abspath + ':' + _tree_name]
        # break

    _pimunu_files = []
    for _f in os.listdir(os.path.join(
            _data_file, 'pimunu')):
        _file_abspath = os.path.join(
            _data_file, 'pimunu', _f)
        _pimunu_files += [_file_abspath + ':' + _tree_name]
        # break
    
    log.info('Creating pienu arrays')
    _pienu_arrays = uproot.concatenate(_pienu_files, FIELDS, library='ak')
    log.info('Creating pimunu arrays')
    _pimunu_arrays = uproot.concatenate(_pimunu_files, FIELDS, library='ak')

    
    cutflow_pimunu = [('init', len(_pimunu_arrays))]
    cutflow_pienu = [('init', len(_pienu_arrays))]

    # ---
    _pimunu_arrays = _pimunu_arrays[ak.num(_pimunu_arrays['decay.volume']) == 2]
    _pienu_arrays = _pienu_arrays[ak.num(_pienu_arrays['decay.volume']) == 1]
    cutflow_pimunu += [('two-decays-recorded', len(_pimunu_arrays))]
    cutflow_pienu += [('one-decay-recorded', len(_pienu_arrays))]
    efficiency, rejection = plot_atar_energy(_pienu_arrays, _pimunu_arrays, '1_after_decay_cut')
    eff_rej_table = [('1_after_decay_cut', efficiency, rejection)]

    # ---
    _pimunu_arrays = _pimunu_arrays[np.fabs(_pimunu_arrays['decay.volume'][:,0] - 100000) < 50000]
    _pienu_arrays =_pienu_arrays[np.fabs(_pienu_arrays['decay.volume'][:,0] - 100000) < 50000]
    cutflow_pimunu += [('pion-decay-in-atar', len(_pimunu_arrays))]
    cutflow_pienu += [('pion-decay-in-atar', len(_pienu_arrays))]
    efficiency, rejection = plot_atar_energy(_pienu_arrays, _pimunu_arrays, '2_pion_decay_in_atar')
    eff_rej_table += [('2_pion_decay_in_atar', efficiency, rejection)]

    # ---
    _pimunu_arrays = _pimunu_arrays[np.fabs(_pimunu_arrays['decay.volume'][:,1] - 100000) < 50000]
    cutflow_pimunu += [('muon-decay-in-atar', len(_pimunu_arrays))]
    efficiency, rejection = plot_atar_energy(_pienu_arrays, _pimunu_arrays, '3_muon_decay_in_atar')
    eff_rej_table += [('3_muon_decay_in_atar', efficiency, rejection)]

    # ---
    _pimunu_arrays = _pimunu_arrays[np.fabs(1e3*_pimunu_arrays['decay.motherEnergy'][:,0] - _pion_mass) < 2]
    _pienu_arrays = _pienu_arrays[np.fabs(1e3*_pienu_arrays['decay.motherEnergy'][:,0] - _pion_mass) < 2]
    cutflow_pimunu += [('pion-decay-at-rest', len(_pimunu_arrays))]
    cutflow_pienu += [('pion-decay-at-rest', len(_pienu_arrays))]
    efficiency, rejection = plot_atar_energy(_pienu_arrays, _pimunu_arrays, '4_pion_decay_at_rest')
    eff_rej_table += [('4_pion_decay_at_rest', efficiency, rejection)]
    
    # ---
    _pimunu_arrays = _pimunu_arrays[np.fabs(1e3*_pimunu_arrays['decay.motherEnergy'][:,1] - _muon_mass) < 2]
    cutflow_pimunu += [('muon-decay-at-rest', len(_pimunu_arrays))]
    efficiency, rejection = plot_atar_energy(_pienu_arrays, _pimunu_arrays, '5_muon_decay_at_rest')
    eff_rej_table += [('5_muon_decay_at_rest', efficiency, rejection)]

    # --- 
    _pimunu_arrays = _pimunu_arrays[np.fabs(ak.max(_pimunu_arrays['atar.time'], axis=1) - ak.min(_pimunu_arrays['atar.time'], axis=1)) > 3]
    _pienu_arrays = _pienu_arrays[np.fabs(ak.max(_pienu_arrays['atar.time'], axis=1) - ak.min(_pienu_arrays['atar.time'], axis=1)) > 3]
    cutflow_pimunu += [('atar hits: |first - last| > 3 ns', len(_pimunu_arrays))]
    cutflow_pienu += [('atar hits: |first - last| > 3 ns', len(_pienu_arrays))]
    efficiency, rejection = plot_atar_energy(_pienu_arrays, _pimunu_arrays, '6_at_least_3ns_total_time_in_atar')
    eff_rej_table += [('6_at_least_3ns_total_time_in_atar', efficiency, rejection)]


    # -- atar computations
    compute_atar_ids(_pimunu_arrays)
    compute_atar_ids(_pienu_arrays)

    _pimunu_atar_arrays = ak.zip({
        _f: _pimunu_arrays[_f] for _f in filter(lambda n: 'atar.' in n, _pimunu_arrays.fields)})
    _pienu_atar_arrays = ak.zip({
        _f: _pienu_arrays[_f] for _f in filter(lambda n: 'atar.' in n, _pienu_arrays.fields)})


    _pimunu_ecentral, _pimunu_etot = energy_centrality(_pimunu_atar_arrays)
    _pienu_ecentral, _pienu_etot = energy_centrality(_pienu_atar_arrays)


    # --- 
    _pimunu_arrays = _pimunu_arrays[_pimunu_ecentral / _pimunu_etot > 0.75]
    _pienu_arrays = _pienu_arrays[_pienu_ecentral / _pienu_etot > 0.75]
    cutflow_pimunu += [('central / total energy > 0.75', len(_pimunu_arrays))]
    cutflow_pienu += [('central / total energy > 0.75', len(_pienu_arrays))]
    efficiency, rejection = plot_atar_energy(_pienu_arrays, _pimunu_arrays, '7_central_energy_0.75')
    eff_rej_table += [('7_central_energy_0.75', efficiency, rejection)]

    print(tabulate.tabulate(cutflow_pimunu, headers=['cut', 'pimunu yields']))
    print()
    print(tabulate.tabulate(cutflow_pienu, headers=['cut', 'pienu yields']))
    print()
    print(tabulate.tabulate(eff_rej_table, headers=['step', 'efficiency', 'rejection']))
    



    
    # _pimunu_first_atar_hit_time = _pimunu_arrays['atar.time'][:,0]
    # _pienu_first_atar_hit_time = _pienu_arrays['atar.time'][:,0]
    # _pimunu_atar_intime = _pimunu_atar_arrays[np.fabs(_pimunu_atar_arrays['atar.time'] - _pimunu_first_atar_hit_time) < 1]
    # _pimunu_atar_intime_central =  _pimunu_atar_intime[_pimunu_atar_intime['atar.pixelid'] > 4]
    # _pimunu_atar_intime_central =  _pimunu_atar_intime_central[_pimunu_atar_intime_central['atar.pixelid'] < 95]
    # _pimunu_atar_intime_central =  _pimunu_atar_intime_central[_pimunu_atar_intime_central['atar.planeid'] > 7]
    # _pimunu_atar_intime_central =  _pimunu_atar_intime_central[_pimunu_atar_intime_central['atar.planeid'] < 40]
    # _pimunu_ecentral =  ak.sum(_pimunu_atar_intime_central['atar.edep'], axis=1)
    # _pimunu_etot =  ak.sum(_pimunu_atar_intime['atar.edep'], axis=1)
    # _pimunu_arrays = _pimunu_arrays[_pimunu_ecentral / _pimunu_etot > 0.75]
    # cutflow_pimunu += [('central / total energy > 0.75', len(_pimunu_arrays))]
    # print(tabulate.tabulate(cutflow_pimunu))

    
    # use last ATAR hit as a proxy for tracker hit
    # # _pimunu_last_atar_hit_time = _pimunu_arrays['atar.time'][:,-1]
    # # _pienu_last_atar_hit_time = _pienu_arrays['atar.time'][:,-1]
    # _pimunu_last_atar_hit_time = ak.max(_pimunu_arrays['atar.time'], axis=1)
    # _pienu_last_atar_hit_time = ak.max(_pienu_arrays['atar.time'], axis=1)


    # # remove atar hits within 1ns of the last one
    # _pimunu_atar_energy_wo_pos = _pimunu_arrays['atar.edep'][np.fabs(_pimunu_arrays['atar.time'] - _pimunu_last_atar_hit_time) > 1]
    # _pienu_atar_energy_wo_pos = _pienu_arrays['atar.edep'][np.fabs(_pienu_arrays['atar.time'] - _pienu_last_atar_hit_time) > 1]

    # _pimunu_energy_sum = ak.sum(_pimunu_atar_energy_wo_pos, axis=1)
    # _pienu_energy_sum = ak.sum(_pienu_atar_energy_wo_pos, axis=1)

    # fig = plt.figure()
    # ax = fig.add_subplot()
    # ax.set_xlabel('ATAR Energy [MeV]')
    # ax.set_ylabel('Entries')
    # plt.hist(_pimunu_energy_sum, bins=100, range=(0, 20), color='blue', label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
    # plt.hist(_pienu_energy_sum, bins=100, range=(0, 20), color='red', alpha=0.5,  label=r'$\pi^{+}\to e^{+}\nu_{e}$')
    # plt.yscale('log')
    # plt.legend()
    # plt.savefig(os.path.join(
    #     'plots', 'atar_energy_after_cut.pdf'))







    # _pimunu_last_atar_hit_time = ak.max(_pimunu_arrays['atar.time'], axis=1)
    # _pimunu_atar_energy_wo_pos = _pimunu_arrays['atar.edep'][np.fabs(_pimunu_arrays['atar.time'] - _pimunu_last_atar_hit_time) > 2]
    # _pimunu_energy_sum = ak.sum(_pimunu_atar_energy_wo_pos, axis=1)


    # _pimunu_arrays = _pimunu_arrays[_pimunu_energy_sum < 12.5]
    
    # _pimunu_pion_decay_x = ak.firsts(_pimunu_arrays['decay.posX'])
    # _pimunu_pion_decay_y = ak.firsts(_pimunu_arrays['decay.posY'])
    # _pimunu_pion_decay_r = np.sqrt(_pimunu_pion_decay_x*_pimunu_pion_decay_x + _pimunu_pion_decay_y*_pimunu_pion_decay_y)
    # _pimunu_pion_decay_z = ak.firsts(_pimunu_arrays['decay.posZ'])

    # fig = plt.figure()
    # ax = fig.add_subplot()
    # ax.set_xlabel('Pion Longitudinal Decay Position [mm]')
    # ax.set_ylabel('Pion Radial Decay Position [mm]')
    # plt.hist2d(
    #     ak.flatten(_pimunu_pion_decay_z, axis=0).to_numpy(),
    #     ak.flatten(_pimunu_pion_decay_r, axis=0).to_numpy(),
    #     bins=(100, 100),
    #     range=((0, 8), (0, 15)),
    #     label=r'$\pi^{+}\to \mu^{+}\nu_{\mu} \to e^{+}\nu_{\mu}\nu_{e}$')
