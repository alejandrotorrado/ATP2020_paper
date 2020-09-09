import numpy as np
import scipy.io as sio
import sklearn.ensemble as skl
import argparse
import pickle
from os.path import expanduser, join

parser = argparse.ArgumentParser(description='Train ML model using input training data from file.')
parser.add_argument('--file','-f',type=str,required=True,help='Full path to file containing training data.')
parser.add_argument('--ntrees','-nt',type=int,required=False,default=200,help='Number of classification trees to use.')
parser.add_argument('--savedir','-sd',type=str,required=True,help='Path to temp directory in which to save outputs.')
args = parser.parse_args()

# Retrieve file path from command line arguments
test_data_file = args.file
nTrees = args.ntrees
out_dir = args.savedir

print('Loading training data from: {}'.format(test_data_file))
print('Using {} classification trees for Random Forest model.'.format(nTrees))

load_test_data = sio.loadmat(test_data_file)

training_data = load_test_data['training']
training_scores = np.squeeze(load_test_data['scoredclass'])
print('feature array shape: {}'.format(training_data.shape))
print('class array shape: {}'.format(training_scores.shape))

RF = skl.RandomForestClassifier(nTrees,oob_score=True)
Mdl = RF.fit(training_data,training_scores)
oob_error = Mdl.oob_score_
feature_importances = Mdl.feature_importances_

print('RF OOB score is: {:.2f}%'.format(oob_error * 100))

pkl_filename = join(out_dir,'pyMdl.pkl')

with open(pkl_filename,'wb') as pkl_out:
	pickle.dump(Mdl,pkl_out,-1)

oob_matfile = join(out_dir,'OOB_error.mat')
print(oob_matfile)
sio.savemat(oob_matfile,mdict={'OOB_error': oob_error})

fi_matfile = join(out_dir,'feature_importances.mat')
print(fi_matfile)
sio.savemat(fi_matfile,mdict={'feature_importances': feature_importances})


