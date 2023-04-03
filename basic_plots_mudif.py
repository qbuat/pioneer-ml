#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import logging

from core.io import get_array
from core.selector import select_atar_dar
from atar.utils import compute_atar_ids
from core import log; log.getChild(__name__)


import awkward as ak
from matplotlib import pyplot as plt
import numpy as np

def plot_1d(arr, xlabel='', plot_name='plot.pdf', log_y=False):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Entries')
    plt.hist(arr, bins=100, color='red', label=r'$\pi^{+}\to \mu^{+}\nu_{\mu} \to e^{+}\nu_{\mu}\nu_{e}$')
    plt.legend()
    if log_y:
        plt.yscale('log')
    
    plt.savefig(os.path.join(
        'plots', plot_name))
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('directory')
    parser.add_argument('--event-type', choices=['pienu', 'pimunu', 'pimudif'], default='pimudif')
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--verbose',action='store_true', default=False)
    args = parser.parse_args()

    if args.debug:
        log.setLevel(logging.DEBUG)
        from core import ch;
        ch.setLevel(logging.DEBUG)
        
    log.info('Starting basic analysis!')
    _array = get_array(args.directory, args.event_type, debug=args.debug)
    # _array = _array[ak.num(_array['decay.volume']) == 2]
    # compute_atar_ids(_array)
    _array = select_atar_dar(_array, verbose=args.verbose, evt_type=args.event_type)


    _vol_pi = _array['decay.volume'][:,0] 
    _x_pi = _array['decay.posX'][:,0] 
    _y_pi = _array['decay.posY'][:,0] 
    _z_pi = _array['decay.posZ'][:,0] 
    _r_pi = np.sqrt(_x_pi*_x_pi + _y_pi*_y_pi)
    _e_pi = _array['decay.motherEnergy'][:,1]
    
    _vol_mu = _array['decay.volume'][:,1] 
    _x_mu = _array['decay.posX'][:,1] 
    _y_mu = _array['decay.posY'][:,1] 
    _z_mu = _array['decay.posZ'][:,1] 
    _r_mu = np.sqrt(_x_mu*_x_mu + _y_mu*_y_mu)


    plot_1d(_e_pi, xlabel='Pion Energy Before Decay [MeV]', plot_name='pidif_pion_ene.pdf')
    plot_1d(_r_pi, xlabel='Pion Radial Decay Position [mm]', plot_name='pidif_pion_r.pdf')
    plot_1d(_r_mu, xlabel='Muon Radial Decay Position [mm]', plot_name='pidif_muon_r.pdf')
    plot_1d(_z_pi, xlabel='Pion Longitudinal Decay Position [mm]', plot_name='pidif_pion_z.pdf')
    plot_1d(_z_mu, xlabel='Muon Longitudinal Decay Position [mm]', plot_name='pidif_muon_z.pdf')
    plot_1d(_z_mu - _z_pi, xlabel='Muon - Pion Longitudinal Decay Position [mm]', plot_name='pidif_muon_minus_pion_z.pdf')
    plot_1d(_vol_pi, xlabel='Pion Decay Volume', plot_name='pion_decay_volume.pdf', log_y=True)
    plot_1d(_vol_mu, xlabel='Muon Decay Volume', plot_name='muon_decay_volume.pdf', log_y=True)

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Pion Longitudinal Decay Position [mm]')
    ax.set_ylabel('Pion Radial Decay Position [mm]')
    plt.hist2d(
        _z_pi.to_numpy(),
        _r_pi.to_numpy(),
        bins=(100, 100),
        range=((0, 8), (0, 15)),
        label=r'$\pi^{+}\to \mu^{+}\nu_{\mu} \to e^{+}\nu_{\mu}\nu_{e}$')
    # plt.zscale('log')
    # plt.legend()
    plt.savefig(os.path.join(
        'plots', 'pion_decay_rz_position.pdf'))

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Muon Longitudinal Decay Position [mm]')
    ax.set_ylabel('Muon Radial Decay Position [mm]')
    plt.hist2d(
        _z_mu.to_numpy(),
        _r_mu.to_numpy(),
        bins=(100, 100),
        range=((0, 8), (0, 15)),
        label=r'$\pi^{+}\to \mu^{+}\nu_{\mu} \to e^{+}\nu_{\mu}\nu_{e}$')
    # plt.zscale('log')
    # plt.legend()
    plt.savefig(os.path.join(
        'plots', 'muon_decay_rz_position.pdf'))

    select_atar_dar(_array, verbose=args.verbose, evt_type=args.event_type)


    
