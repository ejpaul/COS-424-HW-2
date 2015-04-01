'''
Created on Mar 23, 2015

@author: jonathanshor
'''
import numpy as np
import sys, time
from optparse import OptionParser
import sklearn.linear_model as sklm
import UtilityFunctions
import FillData

#FEAT_LIST = ['Beta']
FEAT_LIST = 'Beta'

def get_X_y_Xstar_gTruth(train, sample, test, fillwith = "R"):
# Requires a train.bed and a sample.bed
# Accepts a string in ("M","R","K#") method to replace NAN in train. Default is "R".
# "M"=fill_mean, "R"=fill_rand, "K#"=fill_neighbors(#) with # parsed from string 
# Returns arrays X, y that can be used to train a regression model,
# i.e. are of length  = count of non-NAN rows in sample with 1 in last column
# first array of width = FEAT_LIST count, second array of width 1
# And returns arrays Xstar containing the remaining train rows to predict on 
    fillwith = fillwith.upper()
    if fillwith[0] == "K":
        train['Beta'] = FillData.fill_neighbors(np.transpose(np.array((train['Start'], train['End']))), \
                                                train['Beta'], int(fillwith[1:]))
    elif fillwith == "M":
        train['Beta'] = FillData.fill_mean(train['Beta'])
    else:
        train['Beta'] = FillData.fill_rand(train['Beta'])
    X = train[FEAT_LIST][sample['450k']==1]
    sampley = sample[FEAT_LIST][sample['450k']==1]
    X = X[~np.isnan(sampley)]
    y = sampley[~np.isnan(sampley)]

    
    Xstar = train[FEAT_LIST][sample['450k']==0]
    gTruth = test['Beta'][sample['450k']==0]
    Xstar = Xstar[~np.isnan(gTruth)]
    gTruth = gTruth[~np.isnan(gTruth)]


    print "X: %s ,y: %s, X*: %s, gT: %s" % (X.shape,y.shape, Xstar.shape, gTruth.shape)
    return (X, y, Xstar, gTruth)

# Requires a fitted model and an array of samples to test on
# Returns an array of the predictions for samples in xStars
# Perhaps include xStars in returned array for convenience? Would we need them?
#def predict_yHats(model, xStars):
    

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read bed data fom PATH', metavar='PATH')
    (options, args) = parser.parse_args()     
    path = options.path
    print "PATH = " + path
    start_time = time.time()
    
    train = UtilityFunctions.read_bed_dat_train(path)
    sample = UtilityFunctions.read_bed_dat_sample(path)
    test = UtilityFunctions.read_bed_dat_test(path)
    (trainX, sample_y, trainXstars, test_y) = get_X_y_Xstar_gTruth(train.copy(), sample.copy(), test.copy())
    
    lin_start_time = time.time()
    linRegr = sklm.LinearRegression()
    linRegr.fit(trainX, sample_y, n_jobs=-1)
    yHats = linRegr.predict(trainXstars)
    
    (r2, RMSE) = UtilityFunctions.calc_r2_RMSE(yHats, test_y)
    UtilityFunctions.storePreds(path, yHats, "LinearRegression_r2=%.3f_RMSE=%.3f" % (r2,RMSE), lin_start_time)
    print "r^2 = %s, RMSE = %s" % (r2, RMSE)
    print "sklearn.Linear_model.score = %s" % linRegr.score(trainXstars, test_y)
    
    #ridgeRegr = sklm.Ridge (alpha = .5)
    #ridgeRegr.fit(trainX, sample_y)
    
    #lassoRegr = sklm.Lasso(alpha = 0.1)
    
    print "Runtime: %f" % (time.time()-start_time)

if __name__ == '__main__':
    main(sys.argv[1:])