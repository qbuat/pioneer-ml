# ML stuff
from keras.models import Model
from keras.layers import Input, TimeDistributed
from keras.layers.merging import concatenate
from keras.layers.core import Dense, Dropout, Flatten, Reshape, Masking
from keras.layers import LSTM
from keras.optimizers import Adam
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

