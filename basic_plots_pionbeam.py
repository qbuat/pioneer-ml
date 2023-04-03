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
    'calo.Etot', 'calo.theta',
    'decay.motherPDGID', 'decay.time', 'decay.motherEnergy', 'decay.volume',
    'decay.posX', 'decay.posY', 'decay.posZ',
    'decay.daughterPDGID', 'decay.daughterMomX', 'decay.daughterMomY', 'decay.daughterMomZ',
    'init.mom_x', 'init.mom_y',
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


    # _data_file = './data_v3/'
    # _data_file = '/Users/quentin/Desktop/for_vincent/'
    # _data_file = '/Users/qbuat/pioneer/main/data_55_2.5/'
    # _data_file = '/Users/qbuat/Downloads/noDegrader/'

    # _data_file = '/Users/qbuat/pioneer/data_62_2.5'
    # _data_file = '/Users/qbuat/pioneer/data_65_2.5'
    _data_file = '/Users/qbuat/pioneer/test_uniform_beam'
    _tree_name = 'sim'
    _pienu_files = []
    for _f in os.listdir(os.path.join(
            _data_file, 'pienu')):
        _file_abspath = os.path.join(
            _data_file, 'pienu', _f)
        _pienu_files += [_file_abspath + ':' + _tree_name]
        break

    _pimunu_files = []
    for _f in os.listdir(os.path.join(
            _data_file, 'pimunu')):
        _file_abspath = os.path.join(
            _data_file, 'pimunu', _f)
        _pimunu_files += [_file_abspath + ':' + _tree_name]
        break
    
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


    if do_plots:
        log.info('plotting decay stuff')


        _pimunu_pion_decay_x = ak.firsts(_pimunu_arrays['decay.posX'])
        _pimunu_pion_decay_y = ak.firsts(_pimunu_arrays['decay.posY'])
        _pimunu_pion_decay_r = np.sqrt(_pimunu_pion_decay_x*_pimunu_pion_decay_x + _pimunu_pion_decay_y*_pimunu_pion_decay_y)
        _pimunu_pion_decay_z = ak.firsts(_pimunu_arrays['decay.posZ'])

        
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Pion Radial Decay Position [mm]')
        ax.set_ylabel('Entries')
        plt.hist(_pimunu_pion_decay_r, bins=100, color='red', alpha=0.5, label=r'$\pi^{+}\to \mu^{+}\nu_{\mu} \to e^{+}\nu_{\mu}\nu_{e}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'pion_decay_radial_position.pdf'))

        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Pion Longitudinal Decay Position [mm]')
        ax.set_ylabel('Entries')
        plt.hist(_pimunu_pion_decay_z, bins=100, color='red', alpha=0.5, label=r'$\pi^{+}\to \mu^{+}\nu_{\mu} \to e^{+}\nu_{\mu}\nu_{e}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'pion_decay_z_position.pdf'))


        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Pion Longitudinal Decay Position [mm]')
        ax.set_ylabel('Pion Radial Decay Position [mm]')
        plt.hist2d(
            ak.flatten(_pimunu_pion_decay_z, axis=0).to_numpy(),
            ak.flatten(_pimunu_pion_decay_r, axis=0).to_numpy(),
            bins=(100, 100),
            range=((0, 8), (0, 15)),
            label=r'$\pi^{+}\to \mu^{+}\nu_{\mu} \to e^{+}\nu_{\mu}\nu_{e}$')
        # plt.zscale('log')
        # plt.legend()
        plt.savefig(os.path.join(
            'plots', 'pion_decay_rz_position.pdf'))


        
        _pimunu_arrays = _pimunu_arrays[np.fabs(ak.firsts(_pimunu_arrays['decay.volume']) - 100000) < 50000]
        _pienu_arrays  = _pienu_arrays[np.fabs(ak.firsts(_pienu_arrays['decay.volume']) - 100000) < 50000]
        _pion_mass = 139.570 * 1e3
        _pion_ene = ak.firsts(_pimunu_arrays['decay.motherEnergy'])
        # split based on pion mass +/- 2 keV
        _pimunu_arrays_pidar = _pimunu_arrays[np.fabs(1e3*_pion_ene - _pion_mass) < 2]
        _pimunu_arrays_pidif = _pimunu_arrays[np.fabs(1e3*_pion_ene - _pion_mass) > 2]

        
        _anti_electron_pimunu_px =  _pimunu_arrays['decay.daughterMomX'][ _pimunu_arrays['decay.daughterPDGID']==-11]
        _anti_electron_pimunu_py =  _pimunu_arrays['decay.daughterMomY'][ _pimunu_arrays['decay.daughterPDGID']==-11]
        _anti_electron_pimunu_pz =  _pimunu_arrays['decay.daughterMomZ'][ _pimunu_arrays['decay.daughterPDGID']==-11] 
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



        _anti_electron_pimunu_pidar_px =  _pimunu_arrays_pidar['decay.daughterMomX'][ _pimunu_arrays_pidar['decay.daughterPDGID']==-11]
        _anti_electron_pimunu_pidar_py =  _pimunu_arrays_pidar['decay.daughterMomY'][ _pimunu_arrays_pidar['decay.daughterPDGID']==-11]
        _anti_electron_pimunu_pidar_pz =  _pimunu_arrays_pidar['decay.daughterMomZ'][ _pimunu_arrays_pidar['decay.daughterPDGID']==-11] 
        _anti_electron_pimunu_pidar_p = np.sqrt(_anti_electron_pimunu_pidar_px*_anti_electron_pimunu_pidar_px + _anti_electron_pimunu_pidar_py*_anti_electron_pimunu_pidar_py +_anti_electron_pimunu_pidar_pz*_anti_electron_pimunu_pidar_pz)

        _anti_electron_pimunu_pidif_px =  _pimunu_arrays_pidif['decay.daughterMomX'][ _pimunu_arrays_pidif['decay.daughterPDGID']==-11]
        _anti_electron_pimunu_pidif_py =  _pimunu_arrays_pidif['decay.daughterMomY'][ _pimunu_arrays_pidif['decay.daughterPDGID']==-11]
        _anti_electron_pimunu_pidif_pz =  _pimunu_arrays_pidif['decay.daughterMomZ'][ _pimunu_arrays_pidif['decay.daughterPDGID']==-11] 
        _anti_electron_pimunu_pidif_p = np.sqrt(_anti_electron_pimunu_pidif_px*_anti_electron_pimunu_pidif_px + _anti_electron_pimunu_pidif_py*_anti_electron_pimunu_pidif_py +_anti_electron_pimunu_pidif_pz*_anti_electron_pimunu_pidif_pz)

        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Anti-electron True Momentum [MeV]')
        ax.set_ylabel('Entries')
        plt.hist(ak.flatten(_anti_electron_pimunu_pidar_p, axis=None), bins=100, range=(0, 100), color='green', label=r'$\pi^{+} [DAR] \to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.hist(ak.flatten(_anti_electron_pimunu_pidif_p, axis=None), bins=100, range=(0, 100), color='purple', alpha=0.5,  label=r'$\pi^{+} [DIF] \to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'antielectron_momentum_pimunu.pdf'))

        
        
        _anti_muon_pimunu_pidar_px =  _pimunu_arrays_pidar['decay.daughterMomX'][ _pimunu_arrays_pidar['decay.daughterPDGID']==-13]
        _anti_muon_pimunu_pidar_py =  _pimunu_arrays_pidar['decay.daughterMomY'][ _pimunu_arrays_pidar['decay.daughterPDGID']==-13]
        _anti_muon_pimunu_pidar_pz =  _pimunu_arrays_pidar['decay.daughterMomZ'][ _pimunu_arrays_pidar['decay.daughterPDGID']==-13] 
        _anti_muon_pimunu_pidar_p = np.sqrt(_anti_muon_pimunu_pidar_px*_anti_muon_pimunu_pidar_px + _anti_muon_pimunu_pidar_py*_anti_muon_pimunu_pidar_py +_anti_muon_pimunu_pidar_pz*_anti_muon_pimunu_pidar_pz)

        _anti_muon_pimunu_pidif_px =  _pimunu_arrays_pidif['decay.daughterMomX'][ _pimunu_arrays_pidif['decay.daughterPDGID']==-13]
        _anti_muon_pimunu_pidif_py =  _pimunu_arrays_pidif['decay.daughterMomY'][ _pimunu_arrays_pidif['decay.daughterPDGID']==-13]
        _anti_muon_pimunu_pidif_pz =  _pimunu_arrays_pidif['decay.daughterMomZ'][ _pimunu_arrays_pidif['decay.daughterPDGID']==-13] 
        _anti_muon_pimunu_pidif_p = np.sqrt(_anti_muon_pimunu_pidif_px*_anti_muon_pimunu_pidif_px + _anti_muon_pimunu_pidif_py*_anti_muon_pimunu_pidif_py +_anti_muon_pimunu_pidif_pz*_anti_muon_pimunu_pidif_pz)


        # _anti_muon_pimunu_p_pidif = _anti_muon_pimunu_p[np.fabs(_pion_ene - 140) > 1]
        # _anti_muon_pimunu_p_pidar = _anti_muon_pimunu_p[np.fabs(_pion_ene - 140) < 1]
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Anti-muon True Momentum [MeV]')
        ax.set_ylabel('Entries')
        plt.hist(ak.flatten(_anti_muon_pimunu_pidar_p, axis=None), bins=100, range=(0, 100), color='green', label=r'$\pi^{+} [DAR] \to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.hist(ak.flatten(_anti_muon_pimunu_pidif_p, axis=None), bins=100, range=(0, 100), color='purple', alpha=0.5,  label=r'$\pi^{+} [DIF] \to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'antimuon_momentum_pimunu.pdf'))


        _calo_etot_pimunu_pidif = _pimunu_arrays_pidif['calo.Etot']
        _calo_etot_pimunu_pidar = _pimunu_arrays_pidar['calo.Etot']
        _calo_etot_pienu = _pienu_arrays['calo.Etot']

        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Calo Energy [MeV]')
        ax.set_ylabel('Entries')
        plt.hist(ak.flatten(_calo_etot_pimunu_pidar, axis=None), bins=100, range=(0, 100), color='green', label=r'$\pi^{+} [DAR] \to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.hist(ak.flatten(_calo_etot_pimunu_pidif, axis=None), bins=100, range=(0, 100),  color='purple', alpha=0.5,  label=r'$\pi^{+} [DIF] \to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        # plt.hist(ak.flatten(_calo_etot_pienu, axis=None), bins=100, range=(0, 100), color='blue', alpha=0.5, label=r'$\pi^{+} \to e^{+}\nu_{e}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'caloetot_pimunu.pdf'))

        
        _atar_etot_pimunu_pidif = ak.sum(_pimunu_arrays_pidif['atar.edep'], axis=1)
        _atar_etot_pimunu_pidar = ak.sum(_pimunu_arrays_pidar['atar.edep'], axis=1)
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Total ATAR Energy [MeV]')
        ax.set_ylabel('Entries')
        plt.hist(ak.flatten(_atar_etot_pimunu_pidar, axis=None), bins=100, range=(0, 100), color='green', label=r'$\pi^{+} [DAR] \to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.hist(ak.flatten(_atar_etot_pimunu_pidif, axis=None), bins=100, range=(0, 100), color='purple', alpha=0.5,  label=r'$\pi^{+} [DIF] \to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'ataretot_pimunu.pdf'))

        _decay_volume_pimunu_pidif = ak.firsts(_pimunu_arrays_pidif['decay.volume'])
        _decay_volume_pimunu_pidar = ak.firsts(_pimunu_arrays_pidar['decay.volume'])
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('Decay Volume')
        ax.set_ylabel('Entries')
        plt.hist(ak.flatten(_decay_volume_pimunu_pidar, axis=None), bins=100, color='green', label=r'$\pi^{+} [DAR] \to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.hist(ak.flatten(_decay_volume_pimunu_pidif, axis=None), bins=100, color='purple', alpha=0.5,  label=r'$\pi^{+} [DIF] \to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
        plt.yscale('log')
        plt.legend()
        plt.savefig(os.path.join(
            'plots', 'decay_volume_pimunu.pdf'))
