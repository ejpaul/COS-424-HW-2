'''
Created on Mar 23, 2015

@author: jonathanshor
@author: epaul626
'''
import numpy as np
import sys, time
from optparse import OptionParser

PRE_FNAME = "intersected_final_chr"
TRAIN_FNAME = "_cutoff_20_train_revised.bed"
SAMPLE_FNAME = "_cutoff_20_sample.bed"
TEST_FNAME = "_cutoff_20_test.bed"

def calc_r2_RMSE(preds, gTruth, intercept = 0):
# Expects two 'Beta' both of size nX1
# Accepts an intercept and calc's accordingly
# Returns (Coeff of Determination (r^2), RMSE)
	samps = len(preds)
	if samps != len(gTruth):
		print "Incompatible 'Beta's' passed to calc_r2"
		return (np.nan, np.nan)
	else:
		rss = 0.
		tss = 0.
		emean = np.mean(gTruth)
		for i in range(0, samps):
			rss += (gTruth[i] - preds[i] - intercept) ** 2
			tss += (gTruth[i] - emean) ** 2
		mse = rss / samps
		rmse = mse ** 0.5
		r2 = 1 - (rss / tss)
		return (r2, rmse)

def read_bed_dat_sample(mypath, chrom=1):
# Accepts path to location of PRE_FNAME + chrom + SAMPLE_FNAME
# chrom defaults to 1 
	return np.loadtxt(mypath + PRE_FNAME + str(chrom) + SAMPLE_FNAME, dtype=[('Chrom', np.str_, 4), ('Start', np.int32), \
									('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32), ('450k', np.int8)])

def read_bed_dat_test(mypath, chrom=1):
# Accepts path to location of PRE_FNAME + chrom + TEST_FNAME
# chrom defaults to 1 
	return np.loadtxt(mypath + PRE_FNAME + str(chrom) + TEST_FNAME, dtype=[('Chrom', np.str_, 4), ('Start', np.int32), \
									('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32), ('450k', np.int8)])

def read_bed_dat_train(mypath, chrom=1):
# Accepts path to location of PRE_FNAME + chrom + TRAIN_FNAME
# chrom defaults to 1 
	return np.loadtxt(mypath + PRE_FNAME + str(chrom) + TRAIN_FNAME, dtype=[('Chrom', np.str_, 4), ('Start', np.int32), \
									('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32, (33)), ('450k', np.int8)])

#def storePreds(path, yHats, paras, start_time):
# Stores data to disk in path+"predictions_"+paras+"_csec="+time.time()-start_time+".txt"
# path = directory file will be created
# yHats = array that will be entered as text, one row per line
# paras = string to be inserted as part of filename
# start_time = a time.time() to be subtracted at file create to indicate a running time
	#outfile= open(path+"predictions_"+str(paras)+"_csec="+str(int(time.time()-start_time))+".txt", 'w')
	#outfile.write(paras+"\n")
	#outfile.write("\n".join(yHats))
	#outfile.close()
def storePreds(path, yHats, paras, start_time):
# Stores ndarray to txt in path+"predictions_"+paras+"_csec="+time.time()-start_time+".txt"
# paras also written as a comment to first line
# yHats will be formatted as floats 
	np.savetxt(path+"predictions_"+str(paras)+"_csec="+str(int(time.time()-start_time))+".txt", yHats, fmt="%f", header=paras)

def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--path", dest="path", help='read bed data from PATH', metavar='PATH')
	(options, _args) = parser.parse_args()     
	path = options.path
	print "PATH = " + path
	train = read_bed_dat_train(path)
	sample = read_bed_dat_sample(path)

	print "sample[0] = %s" % sample[0]
	print "train[0] = %s" % train[0]
	print "len(train) = %s" % len(train)
	print "train['Beta'][0] = %s" % train['Beta'][0]

if __name__ == '__main__':
	main(sys.argv[1:])
