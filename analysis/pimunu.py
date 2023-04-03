import os
import numpy as np
import awkward as ak
from matplotlib import pyplot as plt

def get_label(evt_type):

    if evt_type == 'pienu':
        return r'$\pi^{+}\to e^{+}\nu_{e}$'
    elif evt_type == 'pimunu':
        return r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$'
    elif evt_type == 'pimudif':
        return 'pi - mu [dif] - e'
    else:
        raise ValueError
    
def get_color(evt_type):

    if evt_type == 'pienu':
        return 'blue'
    elif evt_type == 'pimunu':
        return 'red'
    elif evt_type == 'pimudif':
        return 'purple'
    else:
        raise ValueError
    
def plot_kin(arr):


    delta_time = arr['decay.time'][:,1] - arr['decay.time'][:, 0]
    delta_x = arr['decay.posX'][:,1] - arr['decay.posX'][:, 0]
    delta_y = arr['decay.posY'][:,1] - arr['decay.posY'][:, 0]
    delta_z = arr['decay.posZ'][:,1] - arr['decay.posZ'][:, 0]
    flight_distance = np.sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z)

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Anti-muon flight time')
    ax.set_ylabel('Entries')
    plt.hist(
        delta_time,
        bins=100,
        color='red',
        alpha=0.5,
        label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
    plt.yscale('log')
    plt.legend()
    plt.savefig(os.path.join(
        'plots', 'muon_flight_time.pdf'))
    
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Anti-muon flight distance')
    ax.set_ylabel('Entries')
    plt.hist(
        flight_distance,
        bins=100,
        color='red',
        alpha=0.5,
        label=r'$\pi^{+}\to \mu^{+}\nu_{\mu}\to e^{+}\nu_{e}\nu_{\mu}$')
    plt.yscale('log')
    plt.legend()
    plt.savefig(os.path.join(
        'plots', 'muon_flight_distance.pdf'))



def plot_pion_decay_pos(arr, evt_type):

    _pion_decay_x = arr['decay.posX'][:,0]
    _pion_decay_y = arr['decay.posY'][:,0]
    _pion_decay_r = np.sqrt(_pion_decay_x*_pion_decay_x + _pion_decay_y*_pion_decay_y)
    _pion_decay_z = arr['decay.posZ'][:,0]

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Pion Radial Decay Position [mm]')
    ax.set_ylabel('Entries')
    plt.hist(
        _pion_decay_r,
        bins=100,
        color=get_color(evt_type),
        alpha=0.5,
        label=get_label(evt_type))
    plt.yscale('log')
    plt.legend()
    plt.savefig(os.path.join(
        'plots', 'pion_decay_radial_position_{}.pdf'.format(
            evt_type)))

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Pion Longitudinal Decay Position [mm]')
    ax.set_ylabel('Entries')
    plt.hist(
        _pion_decay_z,
        bins=100,
        color=get_color(evt_type),
        alpha=0.5,
        label=get_label(evt_type))
    plt.yscale('log')
    plt.legend()
    plt.savefig(os.path.join(
        'plots', 'pion_decay_z_position_{}.pdf'.format(
            evt_type)))


    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Pion Longitudinal Decay Position [mm]')
    ax.set_ylabel('Pion Radial Decay Position [mm]')
    plt.hist2d(
        ak.flatten(_pion_decay_z, axis=0).to_numpy(),
        ak.flatten(_pion_decay_r, axis=0).to_numpy(),
        bins=(100, 100),
        range=((0, 8), (0, 15)),
        label=get_label(evt_type))
    # plt.zscale('log')
    # plt.legend()
    plt.savefig(os.path.join(
        'plots', 'pion_decay_rz_position_{}.pdf'.format(evt_type)))


    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Pion Longitudinal Decay Position [mm]')
    ax.set_ylabel('Pion Transverse (X) Decay Position [mm]')
    plt.scatter(
        ak.flatten(_pion_decay_z, axis=0).to_numpy(),
        ak.flatten(_pion_decay_x, axis=0).to_numpy(),
        label=get_label(evt_type))
    # plt.zscale('log')
    # plt.legend()
    plt.savefig(os.path.join(
        'plots', 'pion_decay_xz_position_{}.pdf'.format(evt_type)))

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Pion Longitudinal Decay Position [mm]')
    ax.set_ylabel('Pion Transverse (Y) Decay Position [mm]')
    plt.scatter(
        ak.flatten(_pion_decay_z, axis=0).to_numpy(),
        ak.flatten(_pion_decay_y, axis=0).to_numpy(),
        label=get_label(evt_type))
    # plt.zscale('log')
    # plt.legend()
    plt.savefig(os.path.join(
        'plots', 'pion_decay_yz_position_{}.pdf'.format(evt_type)))
    
    _muon_decay_x = arr['decay.posX'][:,1]
    _muon_decay_y = arr['decay.posY'][:,1]
    _muon_decay_r = np.sqrt(_muon_decay_x*_muon_decay_x + _muon_decay_y*_muon_decay_y)
    _muon_decay_z = arr['decay.posZ'][:,1]


    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Muon Longitudinal Decay Position [mm]')
    ax.set_ylabel('Entries')
    plt.hist(
        _muon_decay_z,
        bins=100,
        color=get_color(evt_type),
        alpha=0.5,
        label=get_label(evt_type))
    plt.yscale('log')
    plt.legend()
    plt.savefig(os.path.join(
        'plots', 'muon_decay_z_position_{}.pdf'.format(
            evt_type)))


    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Muon Longitudinal Decay Position [mm]')
    ax.set_ylabel('Muon Transverse (X) Decay Position [mm]')
    plt.scatter(
        ak.flatten(_muon_decay_z, axis=0).to_numpy(),
        ak.flatten(_muon_decay_x, axis=0).to_numpy(),
        label=get_label(evt_type))
    # plt.zscale('log')
    # plt.legend()
    plt.savefig(os.path.join(
        'plots', 'muon_decay_xz_position_{}.pdf'.format(evt_type)))

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Muon Longitudinal Decay Position [mm]')
    ax.set_ylabel('Muon Transverse (Y) Decay Position [mm]')
    plt.scatter(
        ak.flatten(_muon_decay_z, axis=0).to_numpy(),
        ak.flatten(_muon_decay_y, axis=0).to_numpy(),
        label=get_label(evt_type))
    # plt.zscale('log')
    # plt.legend()
    plt.savefig(os.path.join(
        'plots', 'muon_decay_yz_position_{}.pdf'.format(evt_type)))

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('Muon Longitudinal Decay Position [mm]')
    ax.set_ylabel('Muon Radial Decay Position [mm]')
    plt.hist2d(
        ak.flatten(_muon_decay_z, axis=0).to_numpy(),
        ak.flatten(_muon_decay_r, axis=0).to_numpy(),
        bins=(100, 100),
        range=((0, 8), (0, 15)),
        label=get_label(evt_type))
    # plt.zscale('log')
    # plt.legend()
    plt.savefig(os.path.join(
        'plots', 'muon_decay_rz_position_{}.pdf'.format(evt_type)))


def plot_mu_kin_energy_decay(array, evt_type):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel("Muon Kinetic Energy at Decay [MeV]")
    ax.set_ylabel('Entries')
    plt.hist(
        ak.flatten(array['decay.motherEnergy'][:,1], axis=0),
        bins=100,
        label=get_label(evt_type))
    plt.savefig(os.path.join(
        'plots', 'muon_kin_energy_at_decay_{}.pdf'.format(evt_type)))
