import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import sklearn.ensemble as skl
import argparse
import pickle
import time
from os.path import join, expanduser

parser = argparse.ArgumentParser(description='Load RandomForest model and use it to score input data')
parser.add_argument('--file','-f',type=str,required=True,help='Full path to file containing input data.')
parser.add_argument('--savedir','-sd',type=str,required=True,help='Path to temp directory in which to save outputs.')
parser.add_argument('--model','-mdl',type=str,required=True,help='Full path to file containing trained ML model.')
parser.add_argument('--verbose','-v',action='store_true',required=False,help='Logical toggle for verbose output.')
args = parser.parse_args()
if args.verbose:
	print('Verbose output is ON.')

# Retrieve file path from command line arguments
input_data_file = args.file
print('Loading training data from: {}'.format(input_data_file))

load_input = sio.loadmat(input_data_file)

input_features = load_input['feature_array']

# load model from pickle file
pkl_filename = args.model
with open(pkl_filename,'rb') as pkl_in:
	Mdl = pickle.load(pkl_in)

if args.verbose:
	t0 = time.time()

RF_predicted = Mdl.predict(input_features)
RF_proba = Mdl.predict_proba(input_features)
RF_labels_int = np.array(RF_predicted)
RF_labels = RF_labels_int.astype(float)

savedir = args.savedir
mat_outfile = join(savedir,'pyOut.mat')
sio.savemat(mat_outfile,mdict={'labels': RF_labels, 'confidence': RF_proba})

if args.verbose:
	t1 = time.time()
	t_el = t1-t0
	print('Time elapsed: {:.2f} seconds.'.format(t_el))
	print('Saved output to file: {}.\n'.format(mat_outfile))


#difftest = test_reals.astype(float) - test_predicted.astype(float)

#a=np.where(difftest!=0)[0]
#print('Found {} mismatched values out of {}.'.format(len(a),len(difftest)))