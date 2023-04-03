#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import awkward as ak
import particle
from matplotlib import pyplot as plt
import argparse
import logging

from atar.utils import geometry
from core.io import get_array
from core.selector import select_atar_dar
from core.utils import get_pion_rz
from atar.utils import compute_atar_ids


def plot_rz_map(r, z, suffix='default'):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Pion Longitudinal Decay Position [mm]')
    ax.set_ylabel('Pion Radial Decay Position [mm]')
    plt.hist2d(
        ak.flatten(z, axis=0).to_numpy(),
        ak.flatten(r, axis=0).to_numpy(),
        bins=(100, 100),
        range=((0, 8), (0, 15)),
        label=r'$\pi^{+}\to \mu^{+}\nu_{\mu} \to e^{+}\nu_{\mu}\nu_{e}$')
    plt.savefig(os.path.join(
        'plots', 'pion_decay_rz_position_{}.pdf'.format(suffix)))

def plot_1d(arrs, names, plot_name='default', x_label='x [unit]', title='Degrader Thickness', log_x=False, log_y=False, **kwargs):
    """
    """
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel(x_label)
    ax.set_ylabel('Entries')
    # ax.set_ylabel('Arbitrary Units')
    # for arr, name in zip(arrs, names):
    #     plt.hist(arr, label=name, **kwargs)
    plt.hist(arrs, label=names, **kwargs)
    if log_y:
        plt.yscale('log')
    if log_x:
        plt.xscale('log')
    plt.legend()
    plt.title(title)
    plt.savefig(os.path.join(
        'plots', '{}.pdf'.format(plot_name)))
    

if __name__ == '__main__':


    debug = False

    from core import log; log.getChild(__name__)
    if debug:
        log.setLevel(logging.DEBUG)
        from core import ch;
        ch.setLevel(logging.DEBUG)
        
    log.info('Starting basic analysis!')

    _pimunu_nodegrade = get_array(
        '/Users/qbuat/pioneer/data_55_2.5_nodegrader/',
        'pimunu',
        debug=debug)

    _pimunu_degrade025 = get_array(
        '/Users/qbuat/pioneer/data_65_2.5_degrader0.25/',
        'pimunu',
        debug=debug)

    _pimunu_degrade05 = get_array(
        '/Users/qbuat/pioneer/data_65_2.5_degrader0.5/',
        'pimunu',
        debug=debug)

    _pimunu_degrade075 = get_array(
        '/Users/qbuat/pioneer/data_65_2.5/',
        'pimunu',
        debug=debug)
    
    compute_atar_ids(_pimunu_degrade025)
    compute_atar_ids(_pimunu_degrade05)
    compute_atar_ids(_pimunu_degrade075)
    _pimunu_nodegrade = select_atar_dar(_pimunu_nodegrade, verbose=True)
    _pimunu_degrade025 = select_atar_dar(_pimunu_degrade025, verbose=True)
    _pimunu_degrade05 = select_atar_dar(_pimunu_degrade05, verbose=True)
    _pimunu_degrade075 = select_atar_dar(_pimunu_degrade075, verbose=True)

    # # atar only stuff
    # _pimunu_atar_arrays = ak.zip({
    #     _f: _pimunu_arrays[_f] for _f in filter(lambda n: 'atar.' in n, _pimunu_arrays.fields)})
    # _pienu_atar_arrays = ak.zip({
    #     _f: _pienu_arrays[_f] for _f in filter(lambda n: 'atar.' in n, _pienu_arrays.fields)})


    log.info('plotting decay stuff')

    _nodegrade_r, _nodegrade_z = get_pion_rz(_pimunu_nodegrade)
    plot_rz_map(_nodegrade_r, _nodegrade_z, suffix='nodegrader')

    _degrade025_r, _degrade025_z = get_pion_rz(_pimunu_degrade025)
    plot_rz_map(_degrade025_r, _degrade025_z, suffix='degrader025')
    
    _degrade05_r, _degrade05_z = get_pion_rz(_pimunu_degrade05)
    plot_rz_map(_degrade05_r, _degrade05_z, suffix='degrader05')
    
    _degrade075_r, _degrade075_z = get_pion_rz(_pimunu_degrade075)
    plot_rz_map(_degrade075_r, _degrade075_z, suffix='degrader075')
    
    
    # plot_1d(
    #     [
    #         ak.flatten(_degrade025_r, axis=0).to_numpy(),
    #         ak.flatten(_degrade05_r, axis=0).to_numpy(),
    #         ak.flatten(_degrade075_r, axis=0).to_numpy()
    #     ],
    #     ['0.25', '0.5', '0.75'],
    #     plot_name='pion_decay_radial_positon_degraders',
    #     x_label='Pion Decay Radial Position [mm]',
    #     bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 30, 100, 1000, 10000, 20000],
    #     # bins=310,
    #     # range=(-100, 3000),
    #     log_x=True,
    #     log_y=True,
    #     histtype='step')
    
    plot_1d(
        [
            ak.flatten(_degrade025_r, axis=0).to_numpy(),
            ak.flatten(_degrade05_r, axis=0).to_numpy(),
            ak.flatten(_degrade075_r, axis=0).to_numpy()
        ],
        ['0.25', '0.5', '0.75'],
        plot_name='pion_decay_radial_positon_degraders_zoomedin',
        x_label='Pion Decay Radial Position [mm]',
        bins=[-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        # bins=310,
        # range=(-100, 3000),
        # log_y=True,
        histtype='step')

    # plot_1d(
    #     [
    #         ak.flatten(_degrade025_z, axis=0).to_numpy(),
    #         ak.flatten(_degrade05_z, axis=0).to_numpy(),
    #         ak.flatten(_degrade075_z, axis=0).to_numpy()
    #     ],
    #     ['0.25', '0.5', '0.75'],
    #     plot_name='pion_decay_longitudinal_positon_degraders',
    #     x_label='Pion Decay Longitudinal Position [mm]',
    #     bins=[0, 1, 2, 10, 20, 100, 200, 1000, 2000, 10000, 100000],
    #     log_x=True,
    #     log_y=True,
    #     # bins=310,
    #     # range=(-100, 3000),
    #     histtype='step')

    plot_1d(
        [
            ak.flatten(_degrade025_z, axis=0).to_numpy(),
            ak.flatten(_degrade05_z, axis=0).to_numpy(),
            ak.flatten(_degrade075_z, axis=0).to_numpy()
        ],
        ['0.25', '0.5', '0.75'],
        plot_name='pion_decay_longitudinal_positon_degraders_zoomedin',
        x_label='Pion Decay Longitudinal Position [mm]',
        bins=[-1, 0, 1, 2, 3, 4, 5, 6, 7],
        log_y=True,
        # bins=310,
        # range=(-100, 3000),
        histtype='step')

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Energy Deposited by charged pion [MeV]')
    ax.set_ylabel('A. U.')
    sum_nodegrade = ak.sum(_pimunu_nodegrade['atar.edep'][_pimunu_nodegrade['atar.pdgid'] == 211], axis = 1)
    sum_degrade05 = ak.sum(_pimunu_degrade05['atar.edep'][_pimunu_degrade05['atar.pdgid'] == 211], axis = 1)
    plt.hist([sum_nodegrade, sum_degrade05], label=['no degrader', 'Degrader (0.5)'], density=True, histtype='step')
    plt.legend()
    plt.savefig(os.path.join(
        'plots', 'pion_energy_degrader_setup.pdf'))

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Energy Deposited by anti-muon [MeV]')
    ax.set_ylabel('A. U.')
    sum_nodegrade = ak.sum(_pimunu_nodegrade['atar.edep'][_pimunu_nodegrade['atar.pdgid'] == -13], axis = 1)
    sum_degrade05 = ak.sum(_pimunu_degrade05['atar.edep'][_pimunu_degrade05['atar.pdgid'] == -13], axis = 1)
    plt.hist([sum_nodegrade, sum_degrade05], label=['no degrader', 'Degrader (0.5)'], bins=[3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5], density=True, histtype='step')
    plt.legend()
    plt.savefig(os.path.join(
        'plots', 'antimuon_energy_degrader_setup.pdf'))
