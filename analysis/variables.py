VARIABLES = {
    'ene_fivepix': {
        'name': 'ene_fivepix',
        'title': 'First Five Positron ATAR Hits',
        'unit': 'MeV',
        'binning': [100, 0, 5],
    },
    'ene_dE_dx': {
        'name': 'ene_dE_dx',
        'title': 'dE/dx (along z direction)',
        'unit': 'MeV/mm',
        'binning': [100, 0, 20],
        },

    'true_dE_dx': {
        'name': 'true_dE_dx',
        'title': 'True dE/dx (first muon step)',
        'unit': 'MeV/mm',
        'binning': [100, 0.01, 10],
        },

    'true_dE_dx_zonly': {
        'name': 'true_dE_dx_zonly',
        'title': 'True dE/dx (first muon step, only along Z)',
        'unit': 'MeV/mm',
        'binning': [100, 0.01, 10],
        },
    
    'delta_z_pion': {
        'name': 'delta_z_pion',
        'title': '#DeltaZ(reco, truth)',
        'unit': 'mm',
        'binning': [100, -0.2, 0.2],
    },

    'delta_dx_last': {
        'name': 'delta_dx_last',
        'title': 'dx in last hit: reco - truth',
        'unit': 'mm',
        'binning': [100, -0.2, 0.2],
    },

    'delta_dx_previous_to_last': {
        'name': 'delta_dx_previous_to_last',
        'title': 'dx in previous to last hit: reco - truth',
        'unit': 'mm',
        'binning': [100, -0.2, 0.2],
    },

    'delta_dz_first_muon': {
        'name': 'delta_dz_first_muon',
        'title': 'dZ(reco) - dZ(truth)',
        'unit': 'mm',
        'binning': [100, -0.4, 0.4],
    },

    'reco_dz_first_muon': {
        'name': 'reco_dz_first_muon',
        'title': 'dZ(reco)',
        'unit': 'mm',
        'binning': [100, -0.5, 0.5],
    },

    'truth_dz_first_muon': {
        'name': 'truth_dz_first_muon',
        'title': 'dZ(truth)',
        'unit': 'mm',
        'binning': [100, -0.5, 0.5],
    },

    'delta_id_pion': {
        'name': 'delta_id_pion',
        'title': '#Delta Pixel ID(reco, truth)',
        'binning': [400, -200, 200]
    },

    'total_energy': {
        'name': 'total_energy',
        'title': 'Energy',
        'unit': 'MeV',
        'binning': [80, 0, 80],
    },

}

VARIABLES_2D = {
    'pion_ene_last_travel': {
        'name': 'pion_ene_last_travel',
        'ylabel': '3D dx',
        'xlabel': 'Energy',
        'yunit' : 'mm',
        'xunit': 'MeV',
        'binning': [200, 0, 2, 200, 0, 0.200],
    },
    
    'pion_ene_last_and_prev_travel': {
        'name': 'pion_ene_last_and_prev_travel',
        'ylabel': '3D dx',
        'xlabel': 'Energy',
        'yunit' : 'mm',
        'xunit': 'MeV',
        'binning': [200, 0, 4, 200, 0, 0.4],
    },

    'pion_ene_previous_to_last_travel': {
        'name': 'pion_ene_previous_to_last_travel',
        'ylabel': '3D dx',
        'xlabel': 'Energy previous to last',
        'yunit' : 'mm',
        'xunit': 'MeV',
        'binning':  [200, 0, 2, 200, 0.12, 0.15],
    },
    'pion_ene_previous_to_last_angle': {
        'name': 'pion_ene_previous_to_last_angle',
        'ylabel': 'Polar Angle',
        'xlabel': 'Energy previous to last',
        'xunit': 'MeV',
        'binning':  [200, 0, 2, 200, 0, 1.5],
    },

    'pion_ene_last_previous_to_last': {
        'name': 'pion_ene_last_previous_to_last',
        'xlabel': 'Energy last hit',
        'ylabel': 'Energy previous-to-last hit',
        'xunit' : 'MeV',
        'yunit': 'MeV',
        'binning': [200, 0, 2, 200, 0, 2],
    },
    
    'pion_dr_last_previous_to_last': {
        'name': 'pion_dr_last_previous_to_last',
        'xlabel': 'dr last',
        'ylabel': 'dr previous-to-last',
        'xunit' : 'mm',
        'yunit': 'mm',
        'binning': [200, 0, .2, 200, 0, .2],
    },

    'pion_theta_last_previous_to_last': {
        'name': 'pion_theta_last_previous_to_last',
        'xlabel': 'last dx Polar Angle',
        'ylabel': 'previous-to-last Polar Angle',
        'binning': [100, 0, 2, 100, 0, 2],
    },
}
    
