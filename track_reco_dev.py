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
import tabulate
import ROOT
import array
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

from atar.utils import geometry
_atar_geo = geometry('/Users/qbuat/pioneer/main/MonteCarlo/geometry/generator/test_output.gdml')
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
    

def atar_2d_projection(arr, ievt):
    # fit two graphs: (z, x) and (z, y)
    # ievt = 0
    if len(arr['atar.pxlID'][ievt]) == 0:
        return
    x, y, z_x,z_y,  e_x, e_y, e_z_x, e_z_y = reco_data(arr[ievt])

    pdgid_x = arr['atar.pdgid'][ievt][arr['atar.planeid'][ievt] %2 == 0]
    x_x_pion = x[pdgid_x == 211]
    z_x_pion = z_x[pdgid_x == 211]
    e_x_x_pion = e_x[pdgid_x == 211]
    e_z_x_pion = e_z_x[pdgid_x == 211]
    gr_pion_z_x = ROOT.TGraphErrors(
        len(x_x_pion),
         z_x_pion,
         x_x_pion,
         e_z_x_pion,
         e_x_x_pion)
    gr_pion_z_x.SetMarkerColor(ROOT.kRed)
    gr_pion_z_x.SetLineColor(ROOT.kRed)
    
    x_x_antimuon = x[pdgid_x == -13]
    z_x_antimuon = z_x[pdgid_x == -13]
    e_x_x_antimuon = e_x[pdgid_x == -13]
    e_z_x_antimuon = e_z_x[pdgid_x == -13]
    gr_antimuon_z_x = ROOT.TGraphErrors(len(x_x_antimuon), z_x_antimuon, x_x_antimuon, e_z_x_antimuon, e_x_x_antimuon)
    gr_antimuon_z_x.SetMarkerColor(ROOT.kBlue)
    gr_antimuon_z_x.SetLineColor(ROOT.kBlue)

    x_x_positron = x[pdgid_x == -11]
    z_x_positron = z_x[pdgid_x == -11]
    e_x_x_positron = e_x[pdgid_x == -11]
    e_z_x_positron = e_z_x[pdgid_x == -11]
    gr_positron_z_x = ROOT.TGraphErrors(len(x_x_positron), z_x_positron, x_x_positron, e_z_x_positron, e_x_x_positron)
    gr_positron_z_x.SetMarkerColor(ROOT.kBlack)
    gr_positron_z_x.SetLineColor(ROOT.kBlack)

    h_template = ROOT.TH1F("", "", 1, -0.3, 0.3)
    h_template.SetBinContent(1, -9999)
    h_template.GetYaxis().SetRangeUser(-1, 1)
    h_template.GetXaxis().SetTitle('Z [cm] (beam direction)')
    h_template.GetYaxis().SetTitle('X [cm] ')
    c = ROOT.TCanvas()
    h_template.Draw()
    gr_pion_z_x.Draw('same PE')
    gr_antimuon_z_x.Draw('same PE')
    gr_positron_z_x.Draw('same PE')
    c.SaveAs('plots/atar_event_display_z_x_event{}.pdf'.format(ievt))

    
    pdgid_y = arr['atar.pdgid'][ievt][arr['atar.planeid'][ievt] %2 != 0]

    y_y_pion = y[pdgid_y == 211]
    z_y_pion = z_y[pdgid_y == 211]
    e_y_y_pion = e_y[pdgid_y == 211]
    e_z_y_pion = e_z_y[pdgid_y == 211]

    y_y_antimuon = y[pdgid_y == -13]
    z_y_antimuon = z_y[pdgid_y == -13]
    e_y_y_antimuon = e_y[pdgid_y == -13]
    e_z_y_antimuon = e_z_y[pdgid_y == -13]

    y_y_positron = y[pdgid_y == -11]
    z_y_positron = z_y[pdgid_y == -11]
    e_y_y_positron = e_y[pdgid_y == -11]
    e_z_y_positron = e_z_y[pdgid_y == -11]

    c1 = ROOT.TCanvas()
    h_template.GetYaxis().SetTitle('X [cm] ')
    h_template.Draw()
    gr_pion = ROOT.TGraphErrors(len(y_y_pion), z_y_pion, y_y_pion, e_z_y_pion, e_y_y_pion)
    gr_pion.SetMarkerColor(ROOT.kRed)
    gr_pion.SetLineColor(ROOT.kRed)
    gr_pion.Draw('same PE')
    if len(y_y_antimuon) !=0:
        gr_antimuon = ROOT.TGraphErrors(len(y_y_antimuon), z_y_antimuon, y_y_antimuon, e_z_y_antimuon, e_y_y_antimuon)
        gr_antimuon.SetMarkerColor(ROOT.kBlue)
        gr_antimuon.SetLineColor(ROOT.kBlue)
        gr_antimuon.Draw('same PE')
    if len(y_y_positron) !=0:
        gr_positron = ROOT.TGraphErrors(len(y_y_positron), z_y_positron, y_y_positron, e_z_y_positron, e_y_y_positron)
        gr_positron.SetMarkerColor(ROOT.kBlack)
        gr_positron.SetLineColor(ROOT.kBlack)
        gr_positron.Draw('same PE')
    c1.SaveAs('plots/atar_event_display_z_y_event{}.pdf'.format(ievt))

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


    _data_file = './data_v3/'
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
    log.info('Creating pimunu arrays')
    _pimunu_arrays = uproot.concatenate(_pimunu_files, FIELDS, library='ak')
    
    cutflow = [('init', len(_pimunu_arrays))]
    _pimunu_arrays = _pimunu_arrays[ak.num(_pimunu_arrays['decay.volume']) == 2]
    cutflow += [('two-decays-recorded', len(_pimunu_arrays))]
    _pimunu_arrays = _pimunu_arrays[np.fabs(_pimunu_arrays['decay.volume'][:,0] - 100000) < 50000]
    cutflow += [('pion-decay-in-atar', len(_pimunu_arrays))]
    _pimunu_arrays = _pimunu_arrays[np.fabs(_pimunu_arrays['decay.volume'][:,1] - 100000) < 50000]
    cutflow += [('muon-decay-in-atar', len(_pimunu_arrays))]

    print(tabulate.tabulate(cutflow))

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



    
    # fit two graphs: (z, x) and (z, y)
    atar_2d_projection(_pimunu_arrays, 0)
    atar_2d_projection(_pimunu_arrays, 1)
    atar_2d_projection(_pimunu_arrays, 2)
    atar_2d_projection(_pimunu_arrays, 3)

    
    # chi2_list = []
    # best_chi2 = -9999
    # best_fit = None
    # best_graph = None
    # for i in range(len(x_x)):
    #     if i  < 1:
    #         continue

    #     gr = ROOT.TGraphErrors(i, z_x[0:i], x_x[0:i], e_z_x[0:i], e_x_x[0:i])
    #     fit = ROOT.TF1("line", '[0]*x+[1]', z_x[0],z_x[i])
    #     res = gr.Fit(fit, 'RSFMQ')
    #     print (i, res.Chi2(), res.Ndf())
    #     c = ROOT.TCanvas()
    #     gr.Draw("AP err0")
    #     fit.SetLineColor(2)
    #     fit.Draw('same')
    #     c.SaveAs('fit_{}.png'.format(i))
