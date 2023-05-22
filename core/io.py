import ROOT
import os
import uproot
import h5py
import numpy as np
import awkward as ak
import json

from . import log; log = log.getChild(__name__)

FIELDS = [
    'atar.pdgid',  'atar.pxlID', 'atar.edep', 'atar.time',
    'calo.Etot', 'calo.theta',
    'decay.motherPDGID', 'decay.time', 'decay.motherEnergy', 'decay.volume',
    'decay.posX', 'decay.posY', 'decay.posZ',
    'decay.daughterPDGID', 'decay.daughterMomX', 'decay.daughterMomY', 'decay.daughterMomZ',
    'init.mom_x', 'init.mom_y',
    # 'track.pdgid', 'track.edep', 'track.post_x', 'track.post_y', 'track.post_z',
    ]


def _array_to_hdf5(ak_array, file_name):
    """
    """
    h5file = h5py.File(file_name, "w")
    group = h5file.create_group("awkward")
    form, length, container = ak.to_buffers(ak_array, container=group)
    group.attrs["form"] = form.to_json()
    group.attrs["length"] = json.dumps(length)


def load_from_cache(name):
    log.warning('{} -- loading awkward array from cache!'.format(name))
    log.warning('{} -- use at your own risk!'.format(name))
    
    _file_name = os.path.join(
        'cache',
        '{}.parquet'.format(name))
    if not os.path.exists(_file_name):
        log.error('{} does not exist!'.format(_file_name))
        return None

    return ak.from_parquet(_file_name)

def _list_files_and_dir(dir):
    _files = []
    _dirs = []
    for _f in os.listdir(dir):
        _fpath = os.path.join(dir, _f)
        _f_abspath = os.path.abspath(_fpath) 
        if os.path.isfile(_f_abspath):
            _files +=[_f_abspath]
        elif os.path.isdir(_f_abspath):
            _dirs += [_f_abspath]
        else:
            raise ValueError('What is {}?'.format(_f_abspath))
    return _files, _dirs

def get_files(dir, depth=3):
    _depth = 0
    _files = []
    _dirs = [dir]
    while (_depth < depth):
        _depth += 1
        _tmp_dirs = []
        for _dir in _dirs:
            _files_in_dir, _dirs_in_dir =  _list_files_and_dir(_dir)
            _files += _files_in_dir
            _tmp_dirs += _dirs_in_dir
        _dirs = _tmp_dirs
    return _files
            
def get_array(path, folder, tree_name='sim', fields=FIELDS, debug=False, use_cache=False):
    """
    """
    if use_cache:
        return load_from_cache(folder)
    
    _files = []
    for _f in sorted(get_files(os.path.join(
            path, folder))):
        if not  _f.endswith('.root'):
            continue
        _files += [_f + ':' + tree_name]

    if debug:
        _files = _files[:1]

    log.debug('List of files:')
    for _f in _files:
        # log.debug(_f)
        log.info(_f)
    log.debug(10 * '-')

    _array = uproot.concatenate(
        _files,
        fields,
        library='ak',
        num_workers=4)
    log.info('Opened {} with {} entries'.format(
        folder, len(_array)))
    return _array



def load(name, cache_folder='cache'):
    fname = os.path.join(cache_folder, "tr_hist{}.root".format(name))
    if not os.path.exists(fname):
        return None
    
    f = ROOT.TFile.Open(fname)
    hist = dict()
    for k in keys:
        h = f.Get(name + k)
        if h:
            print("found: ", k, " name: ", name + k)
            h.SetDirectory(0)
            hist[k] = h
        else:
            print("not found: ", k, " name: ", name + k)
            #if even one histo missing, read everything
            return None
    print("Loaded all histos for ", name)
    f.Close()
    return hist
    
def write(name, hist, cache_folder='cache'):
    rfile_name = 'tr_hist{}.root'.format(name)
    rfile = os.path.join(cache_folder, rfile_name)
    f = ROOT.TFile.Open(rfile, "RECREATE")
    for h in hist.values():
        h.Write()
        h.SetDirectory(0)
    f.Close()
    
