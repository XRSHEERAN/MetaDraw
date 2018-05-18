# These codes have been tested under WinPython-64bit-3.5.2.3Qt5

from model.featurize import Ion2Vector
from keras.models import load_model
from model import fragmentation_config as fconfig
import numpy as np

# init the ion types, max product ion charges and other basic parameters.
config = fconfig.HCD_Config()

# Ion2Vector converts the peptide into numeric features
pep2vec = Ion2Vector(modini = config.modini, ion_types = config.ion_types, cleave_site_num = config.time_step, max_charge = config.max_ion_charge, prev = 1, next = 1, from_CTerm = config.from_CTerm)

def output_predict_with_iontype(pred, peptide, charge):
    # pred is a 2D np.array
    pred_charge = charge-1 if charge <= config.max_ion_charge else config.max_ion_charge
    output = {}
    for i in range(len(config.ion_types)):
        it = config.ion_types[i]
        for ch in range(1, pred_charge+1):
            output['{}+{}'.format(ch, it)] = pred[:len(peptide)-1, i*config.max_ion_charge + ch-1]
    return output
    
pdeep = load_model('./h5/2-layer-BiLSTM-QE-M-M.h5')

# mod_info = '', because pdeep now just support peptides without modifications.
# charge is the precursor charge state
peptide = 'ACDEKFGK'
mod_info = ''
charge = 3

x, __ = pep2vec.FeaturizeOnePeptide(peptide, mod_info, charge)
pred = pdeep.predict(np.array([x]))[0,:,:]

print('from peptide N-term to C-term:\n')
for key, value in output_predict_with_iontype(pred, peptide, charge).items():
    print('{}: {}'.format(key, value))



