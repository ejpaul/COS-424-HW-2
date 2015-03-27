'''
Created on Mar 23, 2015

@author: jonathanshor
'''
import numpy as np
import sys
from optparse import OptionParser
import sklearn.linear_model as sklm

# Requires a train.bed and a sample.bed (or 'Beta'?)
# Returns an array that can be used to train a regression model,
# i.e. is of length  = count of rows in sample with 1 in last column
# omitting the positions we will later predict on
def get_trainable_Xy(train, sample):

# Requires a fitted model and an array of samples to test on
# Returns an array of the predictions for samples in xStars
# Perhaps include xStars in returned array for convenience? Would we need them?
def predict_yHats(model, xStars):
    

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read bed data fom PATH', metavar='PATH')
    (options, args) = parser.parse_args()     
    path = options.path
    print "PATH = " + path

    linRegr = sklm.LinearRegression()
    ridgeRegr = sklm.Ridge (alpha = .5)
    lassoRegr = linear_model.Lasso(alpha = 0.1)

if __name__ == '__main__':
    main(sys.argv[1:])