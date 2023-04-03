#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import logging
import awkward as ak
from core.io import get_array, _array_to_hdf5



    

if __name__ == '__main__':



    data_dir = '/Users/qbuat/pioneer/rocks_data/prod_atar_merged5'
    debug=False
    
    from core import log; log.getChild(__name__)
    if debug:
        log.setLevel(logging.DEBUG)
        from core import ch;
        ch.setLevel(logging.DEBUG)
        
    log.info('Starting basic analysis!')

    _pidar_mudar = get_array(
        data_dir, 'pidar_mudar', debug=debug)
    ak.to_parquet(_pidar_mudar, os.path.join('cache', 'pidar_mudar.parquet'))

    _pidar_mudif = get_array(
        data_dir, 'pidar_mudif', debug=debug)
    ak.to_parquet(_pidar_mudif, os.path.join('cache', 'pidar_mudif.parquet'))
    
    _pidif_mudar = get_array(
        data_dir, 'pidif_mudar', debug=debug)
    ak.to_parquet(_pidif_mudar, os.path.join('cache', 'pidif_mudar.parquet'))

    _pienu = get_array(
        data_dir, 'pienu', debug=debug)
    ak.to_parquet(_pienu, os.path.join('cache', 'pienu.parquet'))
    
    
