#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import logging

from core.io import get_array
from core.selector import select_atar_dar
from atar.utils import compute_atar_ids
from core import log; log.getChild(__name__)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('directory')
    parser.add_argument('--event-type', choices=[
        'pienu',
        'pimunu',
        'pidar_mudar',
        'pidif_mudar',
        'pidar_mudif',
    ], default='pimunu')
    parser.add_argument('--decay-in-atar', default=False, action='store_true')
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--verbose',action='store_true', default=False)
    args = parser.parse_args()

    if args.debug:
        log.setLevel(logging.DEBUG)
        from core import ch;
        ch.setLevel(logging.DEBUG)
        
    log.info('Starting basic analysis!')
    if args.event_type == 'pienu':
        _decay_type = 'pienu'
    elif args.event_type == 'pidar_mudif':
        _decay_type = 'pimudif'
    else:
        _decay_type = 'pimunu'

    _array = get_array(args.directory, args.event_type, debug=args.debug)
    compute_atar_ids(_array)
    _array = select_atar_dar(
        _array,
        evt_type=_decay_type,
        decay_in_atar=args.decay_in_atar,
        verbose=args.verbose)


    # if _decay_type != 'pienu':
    #     # from analysis.pimunu import plot_mu_kin_energy_decay
    #     # plot_mu_kin_energy_decay(_array, args.event_type)
    #     from analysis.pimunu import plot_pion_decay_pos
    #     plot_pion_decay_pos(_array, args.event_type)

    # from atar.ed import plot_atar_event
    # plot_atar_event(_array[10], 'X')
    # plot_atar_event(_array[10], 'Y')

    # from atar.utils import unique_atar
    # unique_atar(_array[10])
