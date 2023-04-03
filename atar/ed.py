import matplotlib as mpl
from matplotlib import pyplot as plt
import awkward as ak
import numpy as np

def plot_atar_event(arr, axis='Y'):

    _atar = ak.zip({
         _f: arr[_f] for _f in filter(lambda n: 'atar.' in n, arr.fields)})

    if axis == 'Y':
        _atar_a = _atar[_atar['atar.planeid']%2 == 0]
    elif axis == 'X':
        _atar_a = _atar[_atar['atar.planeid']%2 != 0]
    else:
        raise ValueError

    _unique_pxl_list, _counts = np.unique(_atar_a['atar.pxlID'], return_counts=True)
    # for _pxl in _atar_y['atar.pxlID']:
    print(_unique_pxl_list)
    print(_counts)

    mpl.rcParams['image.cmap'] = 'cool'
    fig_y = plt.figure()
    ax = fig_y.add_subplot()
    ax.set_xlabel('Layer Number')
    ax.set_ylabel('Strip Number ({})'.format(axis))
    lines = []
    for i in range(50):
        if axis == 'Y' and i % 2 == 0:
            lines += [plt.axvline(i, 0, 100, color='gray', linestyle='--')]
        elif axis == 'X' and i % 2 != 0:
            lines += [plt.axvline(i, 0, 100, color='gray', linestyle='--')]

    colors = (_atar_a['atar.pdgid'] == 211) * 1 + (_atar_a['atar.pdgid'] == -13)*2 + (_atar_a['atar.pdgid'] ==	-11)*3
            # ax.set_ylim(0, 100)
    plt.scatter(_atar_a['atar.planeid'], _atar_a['atar.pixelid'], c=colors, s=100*_atar_a['atar.edep'])
    fig_y.savefig('plots/atar_ed_{}.png'.format(axis))
        
    mpl.rcParams['image.cmap'] = 'viridis'


