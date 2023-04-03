#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import uproot
import re
import xml.etree.ElementTree as ET
import numpy as np
import awkward as ak
import particle
from matplotlib import pyplot as plt
import argparse
import logging

from atar.utils import geometry

FIELDS = [
        'atar.pdgid',  'atar.pxlID', 'atar.edep', 'atar.time',
        'calo.Etot', 
        'decay.motherPDGID', 'decay.time', 'decay.motherEnergy',
        'decay.daughterPDGID', 'decay.daughterMomX', 'decay.daughterMomY', 'decay.daughterMomZ',
        # 'init.vtx_x', 'init.mom_y',
    ]




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


    _data_file = './data_v2/'
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
    print ('number of events = ', len(_pienu_arrays))
    log.info('Creating pimunu arrays')
    _pimunu_arrays = uproot.concatenate(_pimunu_files, FIELDS, library='ak')
    print ('number of events = ', len(_pimunu_arrays))


    ATAR_BASE_ID = 100001
    numPixels = 100
    _pimunu_arrays['atar.planeid'] = ak.values_astype((_pimunu_arrays['atar.pxlID'] - ATAR_BASE_ID) / numPixels, 'int64')
    _pimunu_arrays['atar.pixelid'] = (_pimunu_arrays['atar.pxlID'] - ATAR_BASE_ID) % numPixels
    _pienu_arrays['atar.planeid'] = ak.values_astype((_pienu_arrays['atar.pxlID'] - ATAR_BASE_ID) / numPixels, 'int64')
    _pienu_arrays['atar.pixelid'] = (_pienu_arrays['atar.pxlID'] - ATAR_BASE_ID) % numPixels
    

    # atar only stuff
    _pimunu_atar_arrays = ak.zip({
        _f: _pimunu_arrays[_f] for _f in filter(lambda n: 'atar.' in n, _pimunu_arrays.fields)})
    _pienu_atar_arrays = ak.zip({
        _f: _pienu_arrays[_f] for _f in filter(lambda n: 'atar.' in n, _pienu_arrays.fields)})

    _pimunu_atar_arrays_core = _pimunu_atar_arrays[_pimunu_atar_arrays['atar.planeid'] > 8]        
    _pimunu_atar_arrays_core = _pimunu_atar_arrays_core[_pimunu_atar_arrays_core['atar.planeid'] < 39]        
    _pimunu_atar_arrays_core = _pimunu_atar_arrays_core[_pimunu_atar_arrays_core['atar.pixelid'] > 5]        
    _pimunu_atar_arrays_core = _pimunu_atar_arrays_core[_pimunu_atar_arrays_core['atar.pixelid'] < 94]        
    
    _pimunu_atar_totE = ak.sum(_pimunu_atar_arrays['atar.edep'], axis=1)
    _pimunu_atar_coreE = ak.sum(_pimunu_atar_arrays_core['atar.edep'], axis=1)
    _pimunu_atar_coreE_over_totE = _pimunu_atar_coreE / _pimunu_atar_totE

    _pienu_atar_arrays_core = _pienu_atar_arrays[_pienu_atar_arrays['atar.planeid'] > 8]        
    _pienu_atar_arrays_core = _pienu_atar_arrays_core[_pienu_atar_arrays_core['atar.planeid'] < 39]        
    _pienu_atar_arrays_core = _pienu_atar_arrays_core[_pienu_atar_arrays_core['atar.pixelid'] > 5]        
    _pienu_atar_arrays_core = _pienu_atar_arrays_core[_pienu_atar_arrays_core['atar.pixelid'] < 94]        
    
    _pienu_atar_totE = ak.sum(_pienu_atar_arrays['atar.edep'], axis=1)
    _pienu_atar_coreE = ak.sum(_pienu_atar_arrays_core['atar.edep'], axis=1)
    _pienu_atar_coreE_over_totE = _pienu_atar_coreE / _pienu_atar_totE

    # cut_pienu = _pienu_arrays['atar.planeid'] > 8
    # cut_pienu = cut_pienu == (_pienu_arrays['atar.planeid'] < 39)
    # cut_pienu = cut_pienu == (_pienu_arrays['atar.pixelid'] > 5)
    # cut_pienu = cut_pienu == (_pienu_arrays['atar.pixelid'] < 94)
    
    # _pienu_atar_totE = ak.sum(_pienu_arrays['atar.edep'], axis=1)
    # _pienu_atar_coreE = ak.sum(_pienu_arrays['atar.edep'][cut_pienu], axis=1)
    # _pienu_atar_coreE_over_totE = _pienu_atar_coreE/_pienu_atar_totE

    
    # _pimunu_arrays = _pimunu_arrays[_pimunu_arrays['atar.pxlID']!=0]
    # _pienu_arrays = _pienu_arrays[_pienu_arrays['atar.pxlID']!=0]
    
    _pi_munu_npix = ak.num(_pimunu_arrays['atar.pxlID'])
    _pi_enu_npix = ak.num(_pienu_arrays['atar.pxlID'])

    # _pi_munu_npix = np.array([len(a) for a in _pimunu_arrays['atar.pxlID']])
    # _pi_enu_npix = np.array([len(a) for a in _pienu_arrays['atar.pxlID']])

    _pi_munu_e_t40 = ak.sum(_pimunu_arrays['atar.edep'][_pimunu_arrays['atar.time'] < 40], axis=1)  
    _pi_enu_e_t40 = ak.sum(_pienu_arrays['atar.edep'][_pienu_arrays['atar.time'] < 40], axis=1)  

    
    log.info ('max number of pixel pimunu = {}, pienu = {}'.format(
        max(_pi_munu_npix),  max(_pi_enu_npix)))


    # log.info('Creating and loading ATAR geometry as a python dict')
    # _geo_file = "/Users/quentin/pioneer/PIONEER/MonteCarlo/geometry/generator/test_output.gdml"
    # _atar_geo = geometry(_geo_file)

    if do_plots:
        log.info('plotting decay stuff')
        _anti_muon_px =  _pimunu_arrays['decay.daughterMomX'][ _pimunu_arrays['decay.daughterPDGID']==-13]
        _anti_muon_py =  _pimunu_arrays['decay.daughterMomY'][ _pimunu_arrays['decay.daughterPDGID']==-13]
        _anti_muon_pz =  _pimunu_arrays['decay.daughterMomZ'][ _pimunu_arrays['decay.daughterPDGID']==-13]
        _anti_muon_p = np.sqrt(_anti_muon_px*_anti_muon_px + _anti_muon_py*_anti_muon_py + _anti_muon_pz*_anti_muon_pz)
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Anti-muon True Momentum [MeV]')
        ax.set_ylabel('Entries')
        plt.hist(ak.flatten(_anti_muon_p, axis=None), bins=100, range=(0, 100), color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'antimuon_momentum.pdf'))

        _anti_electron_pimunu_px =  _pimunu_arrays['decay.daughterMomX'][ _pimunu_arrays['decay.daughterPDGID']==-11]
        _anti_electron_pimunu_py =  _pimunu_arrays['decay.daughterMomY'][ _pimunu_arrays['decay.daughterPDGID']==-11]
        _anti_electron_pimunu_pz =  _pimunu_arrays['decay.daughterMomZ'][ _pimunu_arrays['decay.daughterPDGID']==-11] 
        # _anti_electron_pimunu_time         =  _pimunu_arrays['decay.time'][ _pimunu_arrays['decay.daughterPDGID']==-11]
        # _anti_electron_pimunu_motherenergy =  _pimunu_arrays['decay.motherEnergy'][ _pimunu_arrays['decay.daughterPDGID']==-11]
        _anti_muon_px = _anti_muon_px[_anti_electron_pimunu_px != None]
        _anti_muon_py = _anti_muon_py[_anti_electron_pimunu_px != None]
        _anti_muon_pz = _anti_muon_pz[_anti_electron_pimunu_px != None]
        _anti_muon_p = _anti_muon_p[_anti_electron_pimunu_px != None]
        _anti_electron_pimunu_px =  _anti_electron_pimunu_px[_anti_electron_pimunu_px != None]
        _anti_electron_pimunu_py =  _anti_electron_pimunu_py[_anti_electron_pimunu_px != None]
        _anti_electron_pimunu_pz =  _anti_electron_pimunu_pz[_anti_electron_pimunu_px != None]
        _anti_electron_pimunu_p = np.sqrt(_anti_electron_pimunu_px*_anti_electron_pimunu_px + _anti_electron_pimunu_py*_anti_electron_pimunu_py +_anti_electron_pimunu_pz*_anti_electron_pimunu_pz)

        _anti_electron_pienu_px =  _pienu_arrays['decay.daughterMomX'][ _pienu_arrays['decay.daughterPDGID']==-11]
        _anti_electron_pienu_py =  _pienu_arrays['decay.daughterMomY'][ _pienu_arrays['decay.daughterPDGID']==-11]
        _anti_electron_pienu_pz =  _pienu_arrays['decay.daughterMomZ'][ _pienu_arrays['decay.daughterPDGID']==-11]
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
            'plots', 'antielectron_momentum.pdf'))


        fig = plt.figure()
        cos_theta_pimunu = (ak.flatten(_anti_muon_px, axis=None)*ak.flatten(_anti_electron_pimunu_px, axis=None) + ak.flatten(_anti_muon_py, axis=None)*ak.flatten(_anti_electron_pimunu_py, axis=None) + ak.flatten(_anti_muon_pz, axis=None) * ak.flatten(_anti_electron_pimunu_pz, axis=None)) / (ak.flatten(_anti_muon_p, axis=None) * ak.flatten(_anti_electron_pimunu_p, axis=None)) 

        ax = fig.add_subplot()
        ax.set_xlabel('Cos(angle between positron and anti-muon)')
        ax.set_ylabel('Entries')
        plt.hist(ak.flatten(cos_theta_pimunu, axis=None), bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        # plt.yscale('log')
        plt.ylim(bottom=0)
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'cosine_positron_antimuon.pdf'))
        


    if do_plots:
        log.info('plotting number of pixels')
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Total Energy Deposited in ATAR')
        ax.set_ylabel('Entries')
        plt.hist(ak.sum(_pienu_arrays['atar.edep'], axis=1), bins=100, color='blue', label=r'$\pi^{+}\to e^{+}\nu_{e}$')
        plt.hist(ak.sum(_pimunu_arrays['atar.edep'], axis=1), bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'atar_energy.pdf'))


        log.info('plotting number of pixels')
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Total Energy Deposited in ATAR (for events with 75% of energy in core of the ATAR)')
        ax.set_ylabel('Entries')
        plt.hist(ak.sum(_pienu_arrays['atar.edep'][_pienu_atar_coreE_over_totE>0.75], axis=1), bins=100, color='blue', label=r'$\pi^{+}\to e^{+}\nu_{e}$')
        plt.hist(ak.sum(_pimunu_arrays['atar.edep'][_pimunu_atar_coreE_over_totE>0.75], axis=1), bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'atar_energy_coreE_over_0.75totE.pdf'))

        log.info('plotting number of pixels')
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Number of ATAR pixels')
        ax.set_ylabel('Entries')
        plt.hist(_pi_enu_npix, bins=100, color='blue', label=r'$\pi^{+}\to e^{+}\nu_{e}$')
        plt.hist(_pi_munu_npix, bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'number_of_pixels.pdf'))

        log.info('plotting pixel energy in time window < 40ns')
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Sum of Strip energies (with t<40ns)')
        ax.set_ylabel('Entries')
        plt.hist(_pi_enu_e_t40, bins=100, color='blue', label=r'$\pi^{+}\to e^{+}\nu_{e}$')
        plt.hist(_pi_munu_e_t40, bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'pixel_energy_below40ns.pdf'))

        log.info('plotting pixel energy in a 2ns window around the highest energy deposit')
        _time_window = 1 # ns
        _pimunu_atar_hottesthit_time = _pimunu_arrays['atar.time'][ak.argmax(_pimunu_arrays['atar.edep'], axis=1, keepdims=True)]
        _pimunu_atar_time_diff =  _pimunu_arrays['atar.time'] - ak.flatten(_pimunu_atar_hottesthit_time)
        _pimunu_energy_1ns_window = _pimunu_arrays['atar.edep'][np.abs(_pimunu_atar_time_diff) < _time_window]

        _pienu_atar_hottesthit_time = _pienu_arrays['atar.time'][ak.argmax(_pienu_arrays['atar.edep'], axis=1, keepdims=True)]
        _pienu_atar_time_diff =  _pienu_arrays['atar.time'] - ak.flatten(_pienu_atar_hottesthit_time)
        _pienu_energy_1ns_window = _pienu_arrays['atar.edep'][np.abs(_pienu_atar_time_diff) < _time_window]

        
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Sum of Strip Energies witin |\\Delta t| < {0} ns around the most energic strip'.format(_time_window))
        ax.set_ylabel('Entries')
        plt.hist(ak.sum(_pienu_energy_1ns_window, axis=1), bins=100, color='blue', label=r'$\pi^{+}\to e^{+}\nu_{e}$')
        plt.hist(ak.sum(_pimunu_energy_1ns_window, axis=1), bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'pixel_energy_{0}ns.pdf'.format(_time_window)))

        _pimunu_atar_arrays_1ns = _pimunu_atar_arrays[np.abs(_pimunu_atar_time_diff) < _time_window]
        _pimunu_atar_arrays_1ns_core = _pimunu_atar_arrays_1ns[_pimunu_atar_arrays_1ns['atar.planeid'] > 8]        
        _pimunu_atar_arrays_1ns_core = _pimunu_atar_arrays_1ns_core[_pimunu_atar_arrays_1ns_core['atar.planeid'] < 39]        
        _pimunu_atar_arrays_1ns_core = _pimunu_atar_arrays_1ns_core[_pimunu_atar_arrays_1ns_core['atar.pixelid'] > 5]        
        _pimunu_atar_arrays_1ns_core = _pimunu_atar_arrays_1ns_core[_pimunu_atar_arrays_1ns_core['atar.pixelid'] < 94]        
    
        _pimunu_atar_1ns_totE  = ak.sum(_pimunu_atar_arrays_1ns['atar.edep'], axis=1)
        _pimunu_atar_1ns_coreE = ak.sum(_pimunu_atar_arrays_1ns_core['atar.edep'], axis=1)
        _pimunu_atar_1ns_coreE_over_totE = _pimunu_atar_1ns_coreE / _pimunu_atar_1ns_totE

        _pienu_atar_arrays_1ns = _pienu_atar_arrays[np.abs(_pienu_atar_time_diff) < _time_window]
        _pienu_atar_arrays_1ns_core = _pienu_atar_arrays_1ns[_pienu_atar_arrays_1ns['atar.planeid'] > 8]        
        _pienu_atar_arrays_1ns_core = _pienu_atar_arrays_1ns_core[_pienu_atar_arrays_1ns_core['atar.planeid'] < 39]        
        _pienu_atar_arrays_1ns_core = _pienu_atar_arrays_1ns_core[_pienu_atar_arrays_1ns_core['atar.pixelid'] > 5]        
        _pienu_atar_arrays_1ns_core = _pienu_atar_arrays_1ns_core[_pienu_atar_arrays_1ns_core['atar.pixelid'] < 94]        
    
        _pienu_atar_1ns_totE = ak.sum(_pienu_atar_arrays_1ns['atar.edep'], axis=1)
        _pienu_atar_1ns_coreE = ak.sum(_pienu_atar_arrays_1ns_core['atar.edep'], axis=1)
        _pienu_atar_1ns_coreE_over_totE = _pienu_atar_1ns_coreE / _pienu_atar_1ns_totE

        
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Core / Total Energy (witin |\\Delta t| < {0} ns around the most energic strip)'.format(_time_window))
        ax.set_ylabel('Entries')
        plt.hist(_pienu_atar_1ns_coreE_over_totE, bins=100, color='blue', label=r'$\pi^{+}\to e^{+}\nu_{e}$')
        plt.hist(_pimunu_atar_1ns_coreE_over_totE, bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'pixel_energy_ratio_{0}ns.pdf'.format(_time_window)))

        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Total Energy (with 0.75 of energy in core)'.format(_time_window))
        ax.set_ylabel('Entries')
        plt.hist(_pienu_atar_1ns_totE[_pienu_atar_1ns_coreE_over_totE>0.75], bins=100, color='blue', label=r'$\pi^{+}\to e^{+}\nu_{e}$')
        plt.hist(_pimunu_atar_1ns_totE[_pimunu_atar_1ns_coreE_over_totE>0.75], bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'pixel_energy_total_1ns_{0}ns.pdf'.format(_time_window)))

        log.info('plotting number of muon pixels')
        _pi_munu_muonhits = _pimunu_arrays['atar.pxlID'][_pimunu_arrays['atar.pdgid'] == -13]
        _p_munu_uniquehits =  ak.Array([np.unique(_pi_munu_muonhits[i]) for i in range(len(_pi_munu_muonhits))])
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Number of ATAR pixels associated to a true muon')
        ax.set_ylabel('Entries')
        plt.hist(ak.num(_p_munu_uniquehits), bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        # plt.hist(ak.num(_pi_munu_muonhits), bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'number_of_muon_pixels.pdf'))


        
    if do_plots:
        log.info('plotting calo result')
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Energy [MeV]')
        ax.set_ylabel('Entries')
        e_pienu = ak.flatten(_pienu_arrays['calo.Etot'])
        e_pimunu = ak.flatten(_pimunu_arrays['calo.Etot'])
        
        plt.hist(e_pienu, bins=100, color='blue', label=r'$\pi^{+}\to e^{+}\nu_{e}$')
        plt.hist(e_pimunu, bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'calo_energy.pdf'))

        _pienu_muecut = ak.any([_pi_enu_e_t40 > 4.5, _pi_enu_e_t40 < 3.5], axis=0)
        _pimunu_muecut = ak.any([_pi_munu_e_t40 > 4.5,  _pi_munu_e_t40 < 3.5], axis=0)
        e_pienu_muecut = _pienu_arrays['calo.Etot'][_pienu_muecut]
        e_pimunu_muecut = _pimunu_arrays['calo.Etot'][_pimunu_muecut]
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Energy [MeV]')
        ax.set_ylabel('Entries')
        
        plt.hist(ak.flatten(e_pienu_muecut), bins=100, color='blue', label=r'$\pi^{+}\to e^{+}\nu_{e}$')
        plt.hist(ak.flatten(e_pimunu_muecut), bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'calo_energy_after_muon_energy_cut.pdf'))

        _energy_cut = 0.5 #MeV
        e_pienu_muecut = _pienu_arrays['calo.Etot'][np.abs(ak.sum(_pienu_energy_1ns_window, axis=1) - 4.12) > _energy_cut]
        e_pimunu_muecut = _pimunu_arrays['calo.Etot'][np.abs(ak.sum(_pimunu_energy_1ns_window, axis=1) - 4.12) > _energy_cut]
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Energy [MeV]')
        ax.set_ylabel('Entries')
        
        plt.hist(ak.flatten(e_pienu_muecut), bins=100, color='blue', label=r'$\pi^{+}\to e^{+}\nu_{e}$')
        plt.hist(ak.flatten(e_pimunu_muecut), bins=100, color='red', alpha=0.5,  label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'calo_energy_after_muon_energy_cut_v2.pdf'))



        # plt.show()

    if do_plots:
        plt.show()

# print ('mother particles of the decay object')
# pi_plus = ak.flatten(a['decay.motherEnergy'][a['decay.motherPDGID'] == 211])
# print ('Number of charge pions:', len(pi_plus))
# plt.figure()
# plt.hist(pi_plus, bins=100)
#  # plt.show()

# mu_plus = ak.flatten(a['decay.motherEnergy'][a['decay.motherPDGID'] == -13])
# print ('Number of anti muons', len(mu_plus))
# plt.figure()
# plt.hist(mu_plus, bins=100)
# # plt.show()


# # In[13]:


# muonic_neutrino = a['decay.daughterPDGID'][a['decay.daughterPDGID'] == 14]
# print (len(muonic_neutrino))

# mu_plus_x = a['decay.daughterMomX'][a['decay.daughterPDGID'] == -13]
# mu_plus_y = a['decay.daughterMomY'][a['decay.daughterPDGID'] == -13]
# mu_plus_z = a['decay.daughterMomZ'][a['decay.daughterPDGID'] == -13]
# print(mu_plus_x)
# print(ak.flatten(mu_plus_x))
# mu_plus_mom = np.sqrt(mu_plus_x*mu_plus_x + mu_plus_y*mu_plus_y + mu_plus_z*mu_plus_z)
# print (len(mu_plus_mom))


# e_plus = a['decay.daughterPDGID'][a['decay.daughterPDGID'] == -11]
# print (len(e_plus))



# # In[14]:


# ievt = 26

# b = a['atar.pxlID'][ievt] / 100
# b = np.int32(b)
# c = a['atar.pxlID'][ievt][b % 2 ==0]
# d = a['atar.pxlID'][ievt][b % 2 !=0]
# print(len(c))
# plt.figure()
# plt.scatter(
#     c %100,
#     np.int32((c - 100000)/100),
# )
# plt.scatter(
#     d %100,
#     np.int32((d - 100000)/100),
# )
# # plt.show()



# # In[15]:


# b


# # In[16]:

# file = "/Users/qbuat/pioneer/main/MonteCarlo/geometry/generator/test_output.gdml"
# gdmlTree = ET.parse(file)
# gdmlRoot = gdmlTree.getroot()
# atar_geo = {}
# pattern = r'target_pixel_(?P<pixid>10[0-4][0-9][0-9][0-9])intarget_containerpos'
# prog = re.compile(pattern)
# for child in gdmlRoot:
#     #print (child.tag, child.attrib)
#     if child.tag == 'define':
#         for gc in child:
#             if not gc.tag == 'position':
#                 continue
#             if not 'name' in gc.attrib.keys():
#                 raise ValueError
#             if not 'intarget' in gc.attrib['name']:
#                 continue
#             result = prog.match(gc.attrib['name'])
#             atar_geo[int(result.group('pixid'))] = {
#                 'x': float(gc.attrib['x']),
#                 'y': float(gc.attrib['y']),
#                 'z': float(gc.attrib['z']),
#             }
# # print (atar_geo)


# ievt = 26

# atar_pix = a['atar.pxlID'][ievt]
# atar_x = [atar_geo[ID]['x'] for ID in atar_pix]
# atar_y = [atar_geo[ID]['y'] for ID in atar_pix]
# atar_z = [atar_geo[ID]['z'] for ID in atar_pix]
# print (len(atar_x), len(atar_y), len(atar_z))

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')

# ax.scatter(
#     atar_x,
#     atar_y,
#     atar_z
# )


# plt.show()


# # In[18]:


# p = particle.Particle.from_pdgid(-11)
# p

