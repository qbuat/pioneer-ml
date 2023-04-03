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
from core.identification import rnn_model


import tensorflow as tf
from keras.optimizers import Adam
from sklearn.model_selection import train_test_split
from sklearn import metrics

if __name__ == '__main__':


    debug    = False
    do_fit   = False
    do_plots = True

    from core import log; log.getChild(__name__)
    if debug:
        log.setLevel(logging.DEBUG)
        from core import ch;
        ch.setLevel(logging.DEBUG)
        
    log.info('Starting basic analysis!')

    _pimunu = get_array(
        'data_v5/test_uniform_beam',
        # '/Users/qbuat/pioneer/data_65_2.5_degrader0.5/',
        'pimudif',
        debug=debug)

    _pienu = get_array(
        'data_v5/test_uniform_beam',
        # '/Users/qbuat/pioneer/data_65_2.5_degrader0.5/',
        'pienu',
        debug=debug)

    compute_atar_ids(_pimunu)
    compute_atar_ids(_pienu)
    _pimunu = select_atar_dar(_pimunu, evt_type='pimudif', verbose=True)
    _pienu = select_atar_dar(_pienu, evt_type='pienu', verbose=True)

    # atar only stuff
    _pimunu_atar = ak.zip({
        _f: _pimunu[_f] for _f in filter(lambda n: 'atar.' in n, _pimunu.fields)})
    _pienu_atar = ak.zip({
        _f: _pienu[_f] for _f in filter(lambda n: 'atar.' in n, _pienu.fields)})



    _fields = ['atar.planeid', 'atar.pixelid', 'atar.time', 'atar.edep']

    #----
    X_pimunu_odd = []
    X_pimunu_eve = []
    _pimunu_atar_odd = _pimunu_atar[_pimunu_atar['atar.planeid'] % 2 != 0]
    _pimunu_atar_eve = _pimunu_atar[_pimunu_atar['atar.planeid'] % 2 == 0]
    for field in _fields:
        _odd =  ak.fill_none(ak.pad_none(_pimunu_atar_odd[field], 300, clip=True), -999)
        _odd = ak.to_numpy(_odd).reshape(len(_odd), 300, 1)        
        X_pimunu_odd += [_odd]
        _eve =  ak.fill_none(ak.pad_none(_pimunu_atar_eve[field], 300, clip=True), -999)
        _eve = ak.to_numpy(_eve).reshape(len(_eve), 300, 1)        
        X_pimunu_eve += [_eve]
    X_pimunu_odd = np.concatenate(X_pimunu_odd, axis=2)
    X_pimunu_eve = np.concatenate(X_pimunu_eve, axis=2)

    #----
    X_pienu_odd = []
    X_pienu_eve = []
    _pienu_atar_odd = _pienu_atar[_pienu_atar['atar.planeid'] % 2 != 0]
    _pienu_atar_eve = _pienu_atar[_pienu_atar['atar.planeid'] % 2 == 0]
    for field in _fields:
        _odd =  ak.fill_none(ak.pad_none(_pienu_atar_odd[field], 300, clip=True), -999)
        _odd = ak.to_numpy(_odd).reshape(len(_odd), 300, 1)        
        X_pienu_odd += [_odd]
        _eve =  ak.fill_none(ak.pad_none(_pienu_atar_eve[field], 300, clip=True), -999)
        _eve = ak.to_numpy(_eve).reshape(len(_eve), 300, 1)        
        X_pienu_eve += [_eve]
    X_pienu_odd = np.concatenate(X_pienu_odd, axis=2)
    X_pienu_eve = np.concatenate(X_pienu_eve, axis=2)


    X_odd = np.concatenate([X_pienu_odd, X_pimunu_odd])
    X_eve = np.concatenate([X_pienu_eve, X_pimunu_odd])
    
    _target = np.concatenate(
        [
            np.ones(len(_pienu)),
            np.zeros(len(_pimunu)),
        ]
    )
    
    if not (len(X_odd) == len(X_eve) == len(_target)):
        raise ValueError

    _train_idx, _test_idx = train_test_split(
        np.arange(len(X_odd)),
        test_size=0.2,
        random_state=42)

    X_odd_train = X_odd[_train_idx]
    X_odd_test  = X_odd[_test_idx]
    X_eve_train = X_eve[_train_idx]
    X_eve_test  = X_eve[_test_idx]
    y_train = _target[_train_idx]
    y_test = _target[_test_idx]

    log.info('Ready for analysis')

    log.info('creating and compiling the neural network')
    _nn_model = rnn_model()
    optimizer = Adam(learning_rate=1e-5)
    _nn_model.compile(
        loss='binary_crossentropy',
        optimizer=optimizer,
        metrics=['accuracy']
    )

    _model_file = 'cache/rnn_pimudif.h5'
    if do_fit:
        try:
            _history = _nn_model.fit(
                [X_odd_train, X_eve_train],
                y_train,
                epochs=100,
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
        
    if do_plots:
        y_pred = _nn_model.predict([X_odd_test, X_eve_test])
        

        pienu_score = y_pred[y_test == 1]
        pimunu_score = y_pred[y_test == 0]
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('RNN Discriminant Score')
        ax.set_ylabel('Entries')
        
        plt.hist(
            pienu_score,
            bins=100,
            range=(0.8, 1),
            color='blue',
            label='pi-e')
        plt.hist(
            pimunu_score,
            bins=100,
            range=(0.8, 1),
            color='red',
            alpha=0.5,
            label='pi - mu [dif] - e')
        plt.yscale('log')
        plt.legend()


        
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
        fpr = fpr[fpr>0]
        tpr = tpr[np.where(fpr>0)]
        rej = 1 / fpr
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_xlabel('pi-e-nu efficiency')
        ax.set_ylabel('pi [dar] - mu [dif] - e rejection')
        plt.plot(tpr, rej)
        # plt.yscale('log')
        
        plt.show()
        
        
