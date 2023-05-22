


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
