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

# ML stuff
from keras.models import Model
from keras.layers import Input, TimeDistributed
from keras.layers.merging import concatenate
from keras.layers.core import Dense, Dropout, Flatten, Reshape, Masking
from keras.layers import LSTM
from keras.optimizers import Adam
from sklearn.model_selection import train_test_split
import tensorflow as tf

FIELDS = [
        'atar.pdgid',  'atar.pxlID', 'atar.edep', 'atar.time',
        'calo.Etot', 
        # 'decay.motherPDGID', 'decay.motherEnergy', 
        # 'decay.daughterPDGID', 'decay.daughterMomX', 'decay.daughterMomY', 'decay.daughterMomZ',
        # 'init.vtx_x', 'init.mom_y',
    ]



def rnn_model(
        n_pix_odd=300,
        n_pix_even=300,
        n_vars=4,
        masking_val=-999):
    """
    """
    # inputs
    odd_input  = Input(shape=(n_pix_odd, n_vars,))
    even_input = Input(shape=(n_pix_even, n_vars,))

    # odd 
    odd_x = Masking(-999)(odd_input)
    odd_x = TimeDistributed(Dense(32, activation='relu'))(odd_input)
    odd_x = LSTM(128)(odd_x)

    # even 
    even_x = Masking(-999)(even_input)
    even_x = TimeDistributed(Dense(32, activation='relu'))(even_input)
    even_x = LSTM(128)(even_x)

    merged_x = concatenate([odd_x, even_x])
    out_x = Dense(128, activation='relu')(merged_x)
    out_x = Dense(32, activation='relu')(out_x)
    out_x = Dense(1, activation='sigmoid', name='out')(out_x)

    # final model output
    model_input = [
        odd_input,
        even_input,
        ]

    outputs = [
        out_x,
    ]

    _model = Model(inputs=model_input, outputs=outputs)
    return _model


def geometry(gdml_file):
    """
    """
    # file = "/Users/qbuat/pioneer/main/MonteCarlo/geometry/generator/test_output.gdml"
    # file = "/Users/quentin/pioneer/PIONEER/MonteCarlo/geometry/generator/test_output.gdml"
    gdmlTree = ET.parse(gdml_file)
    gdmlRoot = gdmlTree.getroot()

    atar_geo = {}
    pattern = r'target_pixel_(?P<pixid>10[0-9][0-9][0-9][0-9])intarget_containerpos'
    prog = re.compile(pattern)
    for child in gdmlRoot:
        if child.tag != 'define':
            continue
        for gc in child:
            if gc.tag != 'position':
                continue
            if not 'name' in gc.attrib.keys():
                raise ValueError
            if not 'intarget' in gc.attrib['name']:
                continue
            result = prog.match(gc.attrib['name'])
            atar_geo[int(result.group('pixid'))] = {
                'x': float(gc.attrib['x']),
                'y': float(gc.attrib['y']),
                'z': float(gc.attrib['z']),
            }
    # print (atar_geo)
    return atar_geo

def get_atar_pos(pix, variable='x'):
    if pix == -999:
        return -999
    return _atar_geo[pix][variable]

v_get_atar_pos = np.vectorize(get_atar_pos)

if __name__ == '__main__':


    do_plots = False
    do_fit   = False
    

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    log.addHandler(ch)
    log.info('Starting ATAR analysis!')

    log.info('Creating and loading ATAR geometry as a python dict')
    _file = "/Users/quentin/pioneer/PIONEER/MonteCarlo/geometry/generator/test_output.gdml"
    _atar_geo = geometry(_file)

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

    log.info('Creating pimunu arrays')
    _pimunu_arrays = uproot.concatenate(_pimunu_files, FIELDS, library='ak')

    # _pimunu_arrays = _pimunu_arrays[_pimunu_arrays['atar.pxlID']!=0]
    # _pienu_arrays = _pienu_arrays[_pienu_arrays['atar.pxlID']!=0]
    
    _pi_munu_npix = np.array([len(a) for a in _pimunu_arrays['atar.pxlID']])
    _pi_enu_npix = np.array([len(a) for a in _pienu_arrays['atar.pxlID']])
    log.info ('max number of pixel pimunu = {}, pienu = {}'.format(
        max(_pi_munu_npix),  max(_pi_enu_npix)))
    if do_plots:
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

    log.info('Retrieve x, y, z arrays for pimunu')

    _pi_munu_pix_100 = (_pimunu_arrays['atar.pxlID'] - 1)/ 100
    _pi_munu_pix_100 = ak.values_astype(_pi_munu_pix_100, 'int64')
    
    _pi_munu_pix_odd  = _pimunu_arrays['atar.pxlID'][_pi_munu_pix_100 % 2 != 0]
    _pi_munu_pix_even = _pimunu_arrays['atar.pxlID'][_pi_munu_pix_100 % 2 == 0]

    _pi_munu_edep_odd  = _pimunu_arrays['atar.edep'][_pi_munu_pix_100 % 2 != 0]
    _pi_munu_edep_even = _pimunu_arrays['atar.edep'][_pi_munu_pix_100 % 2 == 0]
    _pi_munu_time_odd  = _pimunu_arrays['atar.time'][_pi_munu_pix_100 % 2 != 0]
    _pi_munu_time_even = _pimunu_arrays['atar.time'][_pi_munu_pix_100 % 2 == 0]


    _pi_munu_pix_odd  = ak.fill_none(ak.pad_none(_pi_munu_pix_odd, 300, clip=True), -999)
    _pi_munu_pix_even = ak.fill_none(ak.pad_none(_pi_munu_pix_even, 300, clip=True), -999)
    _pi_munu_edep_odd  = ak.fill_none(ak.pad_none(_pi_munu_edep_odd, 300, clip=True), -999)
    _pi_munu_edep_even = ak.fill_none(ak.pad_none(_pi_munu_edep_even, 300, clip=True), -999)
    _pi_munu_time_odd  = ak.fill_none(ak.pad_none(_pi_munu_time_odd, 300, clip=True), -999)
    _pi_munu_time_even = ak.fill_none(ak.pad_none(_pi_munu_time_even, 300, clip=True), -999)

    _pi_munu_odd_y  = v_get_atar_pos(_pi_munu_pix_odd, variable='y')
    _pi_munu_odd_z  = v_get_atar_pos(_pi_munu_pix_odd, variable='z')
    _pi_munu_even_x = v_get_atar_pos(_pi_munu_pix_even, variable='x')
    _pi_munu_even_z = v_get_atar_pos(_pi_munu_pix_even, variable='z')
    
    log.info('Retrieve x, y, z arrays for pienu')
    _pi_enu_pix_100 = (_pienu_arrays['atar.pxlID'] - 1)/ 100
    _pi_enu_pix_100 = ak.values_astype(_pi_enu_pix_100, 'int64')
    
    _pi_enu_pix_odd  = _pienu_arrays['atar.pxlID'][_pi_enu_pix_100 % 2 != 0]
    _pi_enu_pix_even = _pienu_arrays['atar.pxlID'][_pi_enu_pix_100 % 2 == 0]

    _pi_enu_edep_odd  = _pienu_arrays['atar.edep'][_pi_enu_pix_100 % 2 != 0]
    _pi_enu_edep_even = _pienu_arrays['atar.edep'][_pi_enu_pix_100 % 2 == 0]
    _pi_enu_time_odd  = _pienu_arrays['atar.time'][_pi_enu_pix_100 % 2 != 0]
    _pi_enu_time_even = _pienu_arrays['atar.time'][_pi_enu_pix_100 % 2 == 0]

    _pi_enu_pix_odd  = ak.fill_none(ak.pad_none(_pi_enu_pix_odd, 300, clip=True), -999)
    _pi_enu_pix_even = ak.fill_none(ak.pad_none(_pi_enu_pix_even, 300, clip=True), -999)
    _pi_enu_edep_odd  = ak.fill_none(ak.pad_none(_pi_enu_edep_odd, 300, clip=True), -999)
    _pi_enu_edep_even = ak.fill_none(ak.pad_none(_pi_enu_edep_even, 300, clip=True), -999)
    _pi_enu_time_odd  = ak.fill_none(ak.pad_none(_pi_enu_time_odd, 300, clip=True), -999)
    _pi_enu_time_even = ak.fill_none(ak.pad_none(_pi_enu_time_even, 300, clip=True), -999)

    _pi_enu_odd_y  = v_get_atar_pos(_pi_enu_pix_odd, variable='y')
    _pi_enu_odd_z  = v_get_atar_pos(_pi_enu_pix_odd, variable='z')
    _pi_enu_even_x = v_get_atar_pos(_pi_enu_pix_even, variable='x')
    _pi_enu_even_z = v_get_atar_pos(_pi_enu_pix_even, variable='z')

    log.info('Ready for analysis')

    log.info('creating and compiling the neural network')
    _nn_model = rnn_model()
    optimizer = Adam(learning_rate=1e-5)
    _nn_model.compile(
        loss='binary_crossentropy',
        optimizer=optimizer,
        metrics=['accuracy']
    )

    log.info('format data for training')
    X_odd_enu = np.concatenate(
        [
            ak.to_numpy(_pi_enu_odd_y).reshape(len(_pi_enu_odd_y), 300, 1),
            ak.to_numpy(_pi_enu_odd_z).reshape(len(_pi_enu_odd_z), 300, 1),
            ak.to_numpy(_pi_enu_edep_odd).reshape(len(_pi_enu_pix_odd), 300, 1),
            ak.to_numpy(_pi_enu_time_odd).reshape(len(_pi_enu_pix_odd), 300, 1),
        ],
        axis=2
    )
    X_odd_munu = np.concatenate(
        [
            # ak.to_numpy(_pi_munu_pix_odd).reshape(len(_pi_enu_pix_odd), 300, 1),
            ak.to_numpy(_pi_munu_odd_y).reshape(len(_pi_enu_odd_y), 300, 1),
            ak.to_numpy(_pi_munu_odd_z).reshape(len(_pi_enu_odd_z), 300, 1),
            ak.to_numpy(_pi_munu_edep_odd).reshape(len(_pi_enu_pix_odd), 300, 1),
            ak.to_numpy(_pi_munu_time_odd).reshape(len(_pi_enu_pix_odd), 300, 1),
        ],
        axis=2
    )
    X_odd = np.concatenate([X_odd_enu, X_odd_munu])
        
    
    X_even_enu = np.concatenate(
        [
            # ak.to_numpy(_pi_enu_pix_even).reshape(len(_pi_enu_pix_even), 300, 1),
            ak.to_numpy(_pi_enu_even_x).reshape(len(_pi_enu_even_x), 300, 1),
            ak.to_numpy(_pi_enu_even_z).reshape(len(_pi_enu_even_z), 300, 1),
            ak.to_numpy(_pi_enu_edep_even).reshape(len(_pi_enu_pix_even), 300, 1),
            ak.to_numpy(_pi_enu_time_even).reshape(len(_pi_enu_pix_even), 300, 1),
        ],
        axis=2)
    X_even_munu = np.concatenate(
        [
            ak.to_numpy(_pi_munu_even_x).reshape(len(_pi_enu_even_x), 300, 1),
            ak.to_numpy(_pi_munu_even_z).reshape(len(_pi_enu_even_z), 300, 1),
            ak.to_numpy(_pi_munu_edep_even).reshape(len(_pi_enu_pix_even), 300, 1),
            ak.to_numpy(_pi_munu_time_even).reshape(len(_pi_enu_pix_even), 300, 1)
        ],
        axis=2,
    )
    X_even = np.concatenate([X_even_enu, X_even_munu])
    
    _target = np.concatenate(
        [
            np.zeros(len(_pi_enu_pix_odd)),
            np.ones(len(_pi_munu_pix_even)),
        ]
    )
    
    if not (len(X_odd) == len(X_even) == len(_target)):
        raise ValueError

    _train_idx, _test_idx = train_test_split(
        np.arange(len(X_odd)),
        test_size=0.2,
        random_state=42)

    X_odd_train = X_odd[_train_idx]
    X_odd_test  = X_odd[_test_idx]
    X_even_train = X_even[_train_idx]
    X_even_test  = X_even[_test_idx]
    y_train = _target[_train_idx]
    y_test = _target[_test_idx]

    _model_file = 'rnn.h5'
    if do_fit:
        try:
            _history = _nn_model.fit(
                [X_odd_train, X_even_train],
                y_train,
                epochs=10,
                validation_split=0.1,
                callbacks=[
                    tf.keras.callbacks.EarlyStopping(verbose=True, patience=20, monitor='val_loss'),
                    tf.keras.callbacks.ModelCheckpoint(
                        _model_file,
                        monitor='val_loss',
                        verbose=True, 
                        save_best_only=True)])
            _nn_model.save(_model_file)
        # Allow keyboard interupt to not go over all epochs
        except KeyboardInterrupt:
            print('Ended early...')
    else:
        _nn_model = tf.keras.models.load_model(_model_file)
        
    # _pi_enu_pix = ak.fill_none(ak.pad_none(_pienu_arrays['atar.pxlID'], 600, clip=True), -999)
    # _pi_enu_x = v_get_atar_pos(_pi_enu_pix, variable='x')
    # _pi_enu_y = v_get_atar_pos(_pi_enu_pix, variable='y')
    # _pi_enu_z = v_get_atar_pos(_pi_enu_pix, variable='z')

    #     for f in fields:
    #         print(f, len(a[f][0]), a[f][0])
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

