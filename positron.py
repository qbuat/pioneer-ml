#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import uproot
import numpy as np
import awkward as ak
from matplotlib import pyplot as plt
import argparse
import logging
import tabulate

from atar.utils import geometry, compute_atar_ids, energy_centrality
_atar_geo = geometry('/Users/qbuat/pioneer/main/MonteCarlo/geometry/generator/test_output.gdml')
# _atar_geo = geometry('/Users/quentin/pioneer/PIONEER/MonteCarlo/geometry/generator/test_output.gdml')

FIELDS = [
    'atar.pdgid',  'atar.pxlID', 'atar.edep', 'atar.time',
    # 'calo.Etot', 'calo.theta',
    'decay.motherPDGID', 'decay.time', 'decay.motherEnergy', 'decay.volume',
    'decay.posX', 'decay.posY', 'decay.posZ',
    'decay.daughterPDGID', 'decay.daughterMomX', 'decay.daughterMomY', 'decay.daughterMomZ',
    'init.mom_x', 'init.mom_y',
    ]



def plot_positron(
        _pienu_arrays,
        _pimunu_arrays,
        plot_suffix):
    """
    """
    print ('plotting ' + plot_suffix)
    _pimunu_pidar = _pimunu_arrays[np.fabs(1e3*_pimunu_arrays['decay.motherEnergy'][:,0] - _pion_mass) < 2]
    _pimunu_pidif = _pimunu_arrays[np.fabs(1e3*_pimunu_arrays['decay.motherEnergy'][:,0] - _pion_mass) > 2]

    pimunu_pidar_mudar = _pimunu_pidar[np.fabs(1e3*_pimunu_pidar['decay.motherEnergy'][:,1] - _muon_mass) < 2]
    pimunu_pidif_mudar = _pimunu_pidif[np.fabs(1e3*_pimunu_pidif['decay.motherEnergy'][:,1] - _muon_mass) < 2]

    pimunu_pidar_mudif = _pimunu_pidar[np.fabs(1e3*_pimunu_pidar['decay.motherEnergy'][:,1] - _muon_mass) > 2]
    pimunu_pidif_mudif = _pimunu_pidif[np.fabs(1e3*_pimunu_pidif['decay.motherEnergy'][:,1] - _muon_mass) > 2]
    
    pienu_pidar = _pienu_arrays[np.fabs(1e3*_pienu_arrays['decay.motherEnergy'][:,0] - _pion_mass) < 2]
    pienu_pidif = _pienu_arrays[np.fabs(1e3*_pienu_arrays['decay.motherEnergy'][:,0] - _pion_mass) > 2]

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Positron True Momentum [MeV]')
    ax.set_ylabel('Entries')
    plt.hist(
        ak.flatten(pienu_pidar['positron'][:,0], axis=1),
        bins=100,
        range=(0, 100),
        color='blue',
        label='pi [dar] - e',
        histtype='step')
    plt.hist(
        ak.flatten(pienu_pidif['positron'][:,0], axis=1),
        bins=100,
        range=(0, 100),
        color='royalblue',
        label='pi [dif] - e',
        histtype='step')
    plt.hist(
        ak.flatten(pimunu_pidar_mudar['positron'][:,1], axis=1),
        bins=100,
        range=(0, 100),
        color='darkred',
        alpha=0.5,
        label='pi [dar] - mu [dar] - e',
        histtype='step')
    plt.hist(
        ak.flatten(pimunu_pidif_mudar['positron'][:,1], axis=1),
        bins=100,
        range=(0, 100),
        color='red',
        alpha=0.5,
        label='pi [dif] - mu [dar] - e',
        histtype='step')
    plt.hist(
        ak.flatten(pimunu_pidar_mudif['positron'][:,1], axis=1),
        bins=100,
        range=(0, 100),
        color='green',
        alpha=0.5,
        label='pi [dar] - mu [dif] - e',
        histtype='step')
    plt.hist(
        ak.flatten(pimunu_pidif_mudif['positron'][:,1], axis=1),
        bins=100,
        range=(0, 100),
        color='goldenrod',
        alpha=0.5,
        label='pi [dif] - mu [dif] - e',
        histtype='step')

    plt.yscale('log')
    plt.legend()
    plt.savefig(os.path.join(
        'plots', 'positron_momentum_{}.pdf'.format(plot_suffix)))

    _yields = [
        ('pi [dar] - e', len(pienu_pidar), len(pienu_pidar) / len(_pienu_arrays)),
        ('pi [dif] - e', len(pienu_pidif), len(pienu_pidif) / len(_pienu_arrays)),
        ('pi [dar] - mu [dar] - e', len(pimunu_pidar_mudar), len(pimunu_pidar_mudar) / len(_pimunu_arrays)),
        ('pi [dar] - mu [dif] - e', len(pimunu_pidar_mudif), len(pimunu_pidar_mudif) / len(_pimunu_arrays)),
        ('pi [dif] - mu [dar] - e', len(pimunu_pidif_mudar), len(pimunu_pidif_mudar) / len(_pimunu_arrays)),
        ('pi [dif] - mu [dif] - e', len(pimunu_pidif_mudif), len(pimunu_pidif_mudif) / len(_pimunu_arrays)),
        ]
    print (tabulate.tabulate(_yields, headers=['category', 'yields', 'fraction']))
    
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
    _data_file = '/Users/qbuat/Downloads/noDegrader/'
    # _data_file = '/Users/qbuat/pioneer/main/data_55_2.5/'
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


    _positron_pimunu_px =  _pimunu_arrays['decay.daughterMomX'][ _pimunu_arrays['decay.daughterPDGID']==-11]
    _positron_pimunu_py =  _pimunu_arrays['decay.daughterMomY'][ _pimunu_arrays['decay.daughterPDGID']==-11]
    _positron_pimunu_pz =  _pimunu_arrays['decay.daughterMomZ'][ _pimunu_arrays['decay.daughterPDGID']==-11] 
    _positron_pimunu_p = np.sqrt(_positron_pimunu_px*_positron_pimunu_px + _positron_pimunu_py*_positron_pimunu_py +_positron_pimunu_pz*_positron_pimunu_pz)
    _pimunu_arrays['positron'] = _positron_pimunu_p

    _positron_pienu_px =  _pienu_arrays['decay.daughterMomX'][ _pienu_arrays['decay.daughterPDGID']==-11]
    _positron_pienu_py =  _pienu_arrays['decay.daughterMomY'][ _pienu_arrays['decay.daughterPDGID']==-11]
    _positron_pienu_pz =  _pienu_arrays['decay.daughterMomZ'][ _pienu_arrays['decay.daughterPDGID']==-11]
    _positron_pienu_p = np.sqrt(_positron_pienu_px*_positron_pienu_px + _positron_pienu_py*_positron_pienu_py +_positron_pienu_pz*_positron_pienu_pz)
    _pienu_arrays['positron'] = _positron_pienu_p


    


    plot_positron(
        _pienu_arrays,
        _pimunu_arrays,
        'after_decay_cut')
    
    # # ---
    # _pimunu_arrays = _pimunu_arrays[np.fabs(_pimunu_arrays['decay.volume'][:,0] - 100000) < 50000]
    # _pienu_arrays =_pienu_arrays[np.fabs(_pienu_arrays['decay.volume'][:,0] - 100000) < 50000]
    # cutflow_pimunu += [('pion-decay-in-atar', len(_pimunu_arrays))]
    # cutflow_pienu += [('pion-decay-in-atar', len(_pienu_arrays))]
    # efficiency, rejection = plot_atar_energy(_pienu_arrays, _pimunu_arrays, '2_pion_decay_in_atar')
    # eff_rej_table += [('2_pion_decay_in_atar', efficiency, rejection)]

    # # ---
    # _pimunu_arrays = _pimunu_arrays[np.fabs(_pimunu_arrays['decay.volume'][:,1] - 100000) < 50000]
    # cutflow_pimunu += [('muon-decay-in-atar', len(_pimunu_arrays))]
    # efficiency, rejection = plot_atar_energy(_pienu_arrays, _pimunu_arrays, '3_muon_decay_in_atar')
    # eff_rej_table += [('3_muon_decay_in_atar', efficiency, rejection)]

    # ---
    # _pimunu_arrays = _pimunu_arrays[np.fabs(1e3*_pimunu_arrays['decay.motherEnergy'][:,0] - _pion_mass) < 2]
    # _pienu_arrays = _pienu_arrays[np.fabs(1e3*_pienu_arrays['decay.motherEnergy'][:,0] - _pion_mass) < 2]
    # plot_positron(_pienu_arrays, _pimunu_arrays, '2_pion_decay_at_rest')
    
    # # ---
    # _pimunu_arrays = _pimunu_arrays[np.fabs(1e3*_pimunu_arrays['decay.motherEnergy'][:,1] - _muon_mass) < 2]
    # plot_positron(_pienu_arrays, _pimunu_arrays, '3_muon_decay_at_rest')

    # # --- 
    # _pimunu_arrays = _pimunu_arrays[np.fabs(ak.max(_pimunu_arrays['atar.time'], axis=1) - ak.min(_pimunu_arrays['atar.time'], axis=1)) > 3]
    # _pienu_arrays = _pienu_arrays[np.fabs(ak.max(_pienu_arrays['atar.time'], axis=1) - ak.min(_pienu_arrays['atar.time'], axis=1)) > 3]
    # cutflow_pimunu += [('atar hits: |first - last| > 3 ns', len(_pimunu_arrays))]
    # cutflow_pienu += [('atar hits: |first - last| > 3 ns', len(_pienu_arrays))]
    # efficiency, rejection = plot_atar_energy(_pienu_arrays, _pimunu_arrays, '6_at_least_3ns_total_time_in_atar')
    # eff_rej_table += [('6_at_least_3ns_total_time_in_atar', efficiency, rejection)]


    # # -- atar computations
    # compute_atar_ids(_pimunu_arrays)
    # compute_atar_ids(_pienu_arrays)

    # _pimunu_atar_arrays = ak.zip({
    #     _f: _pimunu_arrays[_f] for _f in filter(lambda n: 'atar.' in n, _pimunu_arrays.fields)})
    # _pienu_atar_arrays = ak.zip({
    #     _f: _pienu_arrays[_f] for _f in filter(lambda n: 'atar.' in n, _pienu_arrays.fields)})


    # _pimunu_ecentral, _pimunu_etot = energy_centrality(_pimunu_atar_arrays)
    # _pienu_ecentral, _pienu_etot = energy_centrality(_pienu_atar_arrays)


    # # --- 
    # _pimunu_arrays = _pimunu_arrays[_pimunu_ecentral / _pimunu_etot > 0.75]
    # _pienu_arrays = _pienu_arrays[_pienu_ecentral / _pienu_etot > 0.75]
    # cutflow_pimunu += [('central / total energy > 0.75', len(_pimunu_arrays))]
    # cutflow_pienu += [('central / total energy > 0.75', len(_pienu_arrays))]
    # efficiency, rejection = plot_atar_energy(_pienu_arrays, _pimunu_arrays, '7_central_energy_0.75')
    # eff_rej_table += [('7_central_energy_0.75', efficiency, rejection)]

    # print(tabulate.tabulate(cutflow_pimunu, headers=['cut', 'pimunu yields']))
    # print()
    # print(tabulate.tabulate(cutflow_pienu, headers=['cut', 'pienu yields']))
    # print()
    # print(tabulate.tabulate(eff_rej_table, headers=['step', 'efficiency', 'rejection']))
    



    
    # # _pimunu_first_atar_hit_time = _pimunu_arrays['atar.time'][:,0]
    # # _pienu_first_atar_hit_time = _pienu_arrays['atar.time'][:,0]
    # # _pimunu_atar_intime = _pimunu_atar_arrays[np.fabs(_pimunu_atar_arrays['atar.time'] - _pimunu_first_atar_hit_time) < 1]
    # # _pimunu_atar_intime_central =  _pimunu_atar_intime[_pimunu_atar_intime['atar.pixelid'] > 4]
    # # _pimunu_atar_intime_central =  _pimunu_atar_intime_central[_pimunu_atar_intime_central['atar.pixelid'] < 95]
    # # _pimunu_atar_intime_central =  _pimunu_atar_intime_central[_pimunu_atar_intime_central['atar.planeid'] > 7]
    # # _pimunu_atar_intime_central =  _pimunu_atar_intime_central[_pimunu_atar_intime_central['atar.planeid'] < 40]
    # # _pimunu_ecentral =  ak.sum(_pimunu_atar_intime_central['atar.edep'], axis=1)
    # # _pimunu_etot =  ak.sum(_pimunu_atar_intime['atar.edep'], axis=1)
    # # _pimunu_arrays = _pimunu_arrays[_pimunu_ecentral / _pimunu_etot > 0.75]
    # # cutflow_pimunu += [('central / total energy > 0.75', len(_pimunu_arrays))]
    # # print(tabulate.tabulate(cutflow_pimunu))

    
    # # use last ATAR hit as a proxy for tracker hit
    # # # _pimunu_last_atar_hit_time = _pimunu_arrays['atar.time'][:,-1]
    # # # _pienu_last_atar_hit_time = _pienu_arrays['atar.time'][:,-1]
    # # _pimunu_last_atar_hit_time = ak.max(_pimunu_arrays['atar.time'], axis=1)
    # # _pienu_last_atar_hit_time = ak.max(_pienu_arrays['atar.time'], axis=1)


    # # # remove atar hits within 1ns of the last one
    # # _pimunu_atar_energy_wo_pos = _pimunu_arrays['atar.edep'][np.fabs(_pimunu_arrays['atar.time'] - _pimunu_last_atar_hit_time) > 1]
    # # _pienu_atar_energy_wo_pos = _pienu_arrays['atar.edep'][np.fabs(_pienu_arrays['atar.time'] - _pienu_last_atar_hit_time) > 1]

    # # _pimunu_energy_sum = ak.sum(_pimunu_atar_energy_wo_pos, axis=1)
    # # _pienu_energy_sum = ak.sum(_pienu_atar_energy_wo_pos, axis=1)

    # # fig = plt.figure()
    # # ax = fig.add_subplot()
    # # ax.set_xlabel('ATAR Energy [MeV]')
    # # ax.set_ylabel('Entries')
    # # plt.hist(_pimunu_energy_sum, bins=100, range=(0, 20), color='blue', label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
    # # plt.hist(_pienu_energy_sum, bins=100, range=(0, 20), color='red', alpha=0.5,  label=r'$\pi^{+}\to e^{+}\nu_{e}$')
    # # plt.yscale('log')
    # # plt.legend()
    # # plt.savefig(os.path.join(
    # #     'plots', 'atar_energy_after_cut.pdf'))







