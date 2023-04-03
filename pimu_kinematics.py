#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import logging

from core.io import get_array
from core.selector import select_pimunu_atar_pidar
from atar.utils import compute_atar_ids
from core import log; log.getChild(__name__)
from analysis.pimunu import plot_kin

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('directory')
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--verbose',action='store_true', default=False)
    args = parser.parse_args()

    if args.debug:
        log.setLevel(logging.DEBUG)
        from core import ch;
        ch.setLevel(logging.DEBUG)
        
    log.info('Starting basic analysis!')
    _array = get_array(args.directory, 'pimunu', debug=args.debug)
    compute_atar_ids(_array)
    _array = select_pimunu_atar_pidar(
        _array, mu_dar=True, verbose=args.verbose)

    # plot_kin(_array)
