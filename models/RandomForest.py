'''
created on Mar 29, 2015

@author: epaul626
'''
import numpy as np
import sys, time
from optparse import OptionParser
from sklearn.ensemble import RandomForestRegressor
sys.path.append('../utils/')
import UtilityFunctions as uf
from FillData import fill_mean, fill_neighbors

NEIGH_FILL = 10

def feat_neighbors(sites, train_beta, train_extras, sample, test):
# Produces 'X' feature array of nearest neighbor beta values that will be fed
# into regressor
# Beta is an array of shape (# of CpG sites, nsamples)
# Sites is an array of shape (# of CpG sites, ) corresponding to 
# the start of each bp site
# train_extras: array of shape (# of CpG sites, 3) which
    # Adds 6 extra features: Exon, DHS, CGI, GC_100, GC_400, GC_1000
    num_extras = train_extras.shape[1]
    X = np.zeros((len(train_beta), 38+num_extras))
    for i in range(0, len(train_beta)):
        # Feature 1: CpG start site
        X[i,0] = sites[i]
        # Find 2 nearest neighbors to site i
        (indices,distance)  = find_neighbors(sites, i)
        index1 = indices[0]
        index2 = indices[1]
        # Feature 2: mean beta at neighbor 1
        X[i, 1] = train_beta[i+index1].sum() / 33.
        # Feature 3: distance to neighbor 1
        X[i, 2] = distance[0]
        # Feature 4: mean beta at neighbor 2
        X[i, 3] = train_beta[i+index2].sum() / 33.
        # Feature 5: distance to neighbor 2
        X[i, 4] = distance[1]
        for j in range(0, 33):
            # Features 6-38: 33 sample beta values at CpG site
            X[i, 5+j] = train_beta[i, j]
    # Features 39+: Extra features
    X[:,-num_extras:] = train_extras

    # Predict on feature set not on 450k chip
    Xstar = X[sample['450k']==0]
    gTruth = test['Beta'][sample['450k']==0]
    Xstar = Xstar[~np.isnan(gTruth)]
    gTruth = gTruth[~np.isnan(gTruth)]

    # Only train on non-NaN values of Y
    X = X[~np.isnan(sample['Beta'])]
    Y = sample['Beta'][~np.isnan(sample['Beta'])]
    return (X, Y, Xstar, gTruth)

def feat_neighbors_island(sites, train_beta, train_extras, sample, test, island=False):
# Produces 'X' feature array of nearest neighbor beta values that will be fed
# into regressor
# Beta is an array of shape (# of CpG sites, nsamples)
# Sites is an array of shape (# of CpG sites, ) corresponding to
# the start of each bp site
# train_extras: array of shape (# of CpG sites, 3) which
    # Adds 6 extra features: Exon, DHS, CGI, GC_100, GC_400, GC_1000
    num_extras = train_extras.shape[1]
    X = np.zeros((len(train_beta), 38+num_extras))
    for i in range(0, len(train_beta)):
        # Feature 1: CpG start site
        X[i,0] = sites[i]
        # Find 2 nearest neighbors to site i
        (indices,distance)  = find_neighbors(sites, i)
        index1 = indices[0]
        index2 = indices[1]
        # Feature 2: mean beta at neighbor 1
        X[i, 1] = train_beta[i+index1].sum() / 33.
        # Feature 3: distance to neighbor 1
        X[i, 2] = distance[0]
        # Feature 4: mean beta at neighbor 2
        X[i, 3] = train_beta[i+index2].sum() / 33.
        # Feature 5: distance to neighbor 2
        X[i, 4] = distance[1]
        for j in range(0, 33):
            # Features 6-38: 33 sample beta values at CpG site
            X[i, 5+j] = train_beta[i, j]
    # Features 39+: Extra features
    X[:,-num_extras:] = train_extras
    
    if island:
        # Cull to strictly CGI sites
        sel_CGI = [sample['CGI'] == 1]
        testc = test[sel_CGI]
        X = X[sel_CGI]
        samplec = sample[sel_CGI]
    else:
        testc = test.copy()
        samplec = sample.copy()

    # Predict on feature set not on 450k chip
    sel_450k = [samplec['450k']==0]
    Xstar = X[sel_450k]
    gTruth = testc['Beta'][sel_450k]
    
    notnan_gt = [~np.isnan(gTruth)]
    Xstar = Xstar[notnan_gt]
    gTruth = gTruth[notnan_gt]

    # Only train on non-NaN values of Y
    notnan_samp = [~np.isnan(samplec['Beta'])]
    X = X[notnan_samp]
    Y = samplec['Beta'][notnan_samp]
    return (X, Y, Xstar, gTruth)


def feat_neighbors_thresh(sites, train_beta, train_extras, sample, test, corr, thresh):
# Produces 'X' feature array of nearest neighbor beta values that will be fed
# into regressor
# Beta is an array of shape (# of CpG sites, nsamples)
# Sites is an array of shape (# of CpG sites, ) corresponding to
# the start of each bp site
# train_extras: array of shape (# of CpG sites, 3) which
    # Adds 6 extra features: Exon, DHS, CGI, GC_100, GC_400, GC_1000
    num_extras = train_extras.shape[1]
    X = np.zeros((len(train_beta), 38+num_extras))
    for i in range(0, len(train_beta)):
        # Feature 1: CpG start site
        X[i,0] = sites[i]
        # Find 2 nearest neighbors to site i
        (indices,distance)  = find_neighbors(sites, i)
        index1 = indices[0]
        index2 = indices[1]
        # Feature 2: mean beta at neighbor 1
        X[i, 1] = train_beta[i+index1].sum() / 33.
        # Feature 3: distance to neighbor 1
        X[i, 2] = distance[0]
        # Feature 4: mean beta at neighbor 2
        X[i, 3] = train_beta[i+index2].sum() / 33.
        # Feature 5: distance to neighbor 2
        X[i, 4] = distance[1]
        for j in range(0, 33):
            # Features 6-38: 33 sample beta values at CpG site
            X[i, 5+j] = train_beta[i, j]
    # Features 39+: Extra features
    X[:,-num_extras:] = train_extras

    # Train and test on sites with > thresh correlation value
    cul0 = [corr[:,0] > thresh]
    X = X[cul0]
    samplec = sample[cul0]
    testc = test[cul0]
    cul_corr = corr[cul0]
    cul1 = [cul_corr[:,1] > thresh]
    X = X[cul1]
    samplec = samplec[cul1]
    testc = testc[cul1]
#     corr = corr[corr[:,1] > thresh]

    # Predict on feature set not on 450k chip
    sel_450k = [samplec['450k']==0]
    Xstar = X[sel_450k]
    gTruth = testc['Beta'][sel_450k]
    
    notnan_gt = [~np.isnan(gTruth)]
    Xstar = Xstar[notnan_gt]
    gTruth = gTruth[notnan_gt]

    # Only train on non-NaN values of Y
    notnan_samp = [~np.isnan(samplec['Beta'])]
    X = X[notnan_samp]
    Y = samplec['Beta'][notnan_samp]
    return (X, Y, Xstar, gTruth)

# Find nearest 2 sites to current CpG site
# Returns list of indices of 2 nearest neighbors
# Where the index in 'sites' is i + (index1, index2)
def find_neighbors(sites, i):
    curr = sites[i]
    if (i > 0):
        d_l1 = curr - sites[i-1]
    if (i > 1):
        d_l2 = curr - sites[i-2]
    if (i < len(sites)-1):
        d_r1 = sites[i+1] - curr
    if (i < len(sites)-2):
        d_r2 = sites[i+2] - curr
    # Look to left
    if (i == 0):
        return ([1, 2], [d_r1, d_r2])
    if (i == 1):
        d_l2 = -1
    # Look to right
    if (i == len(sites)-1):
        return ([-1, -2], [d_l1, d_l2])
    if (i == len(sites)-2):
        d_r2 = -1
    # Compare distances
    if (d_l1 < d_r1 and d_l1 != -1):
        if (d_l2 < d_r1):
            return ([-1, -2], [d_l1, d_l2])
        else:
            return ([-1, 1], [d_l1, d_r1])
    if (d_r2 < d_l1 and d_r2 != -1):
        return ([1, 2], [d_r1, d_r2])
    else:
        return ([1, -1], [d_r1, d_l1])

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read bed data fom PATH', metavar='PATH')
    parser.add_option("-i", dest="island", action="store_true", default=False)
    parser.add_option("-n", dest="fill_neighb", action="store_true", default=False)
    parser.add_option("-m", dest="fill_mean", action="store_true", default=False)
    parser.add_option("-c", dest="chroms", type='int', default=0)
    parser.add_option("-g", dest="GCwindow", type='int', default=0, help='100,400, or 1000. Else all will be included. 0 to exclude feature.')
    parser.add_option("-t", dest="thresh", type='float', default=0.0, help='Correlation value threshold.')
    parser.add_option("-e", dest="estimators", type='int', default=35, help='Number of trees in the forest.')
    (options, _args) = parser.parse_args()
    path = options.path

#     fill_n = options.fill_neighb
#     fill_m = options.fill_mean
    chrom_range = [options.chroms]
    trees = options.estimators # Default 35
#     thresh = options.thresh
    island = options.island # Default False
#     gcwind = options.GCwindow
    start_time = time.time()
        
#     fill_m = False
    fill_vals = [False]
    if options.fill_mean:
        fill_vals += [True]
#     fill_vals = [False, True]
#     gc_vals = [100, 400, 1000, 1500, 0]
    gc_vals = [400, 0, 1500]
    thresh_vals = [0.25, 0.50, 0.75]
#             chrom_range = range(1,22) + [0]
#     chrom_range = [0]
    for chroms in chrom_range:
        # Read in full feature data
        train = uf.read_bed_dat_feat(path, chrom=chroms, ftype='train')
        sample = uf.read_bed_dat_feat(path, chrom=chroms, ftype='sample')
        test = uf.read_bed_dat_feat(path, chrom=chroms, ftype='test')
        print "Beds read %s" % (time.time() - start_time)
        sites = train['Start']

        for gcwind in gc_vals:
    #         options.GCwindow = g
            for fill_m in fill_vals:
                for thresh in thresh_vals:
                    if fill_m:
                        fill_str = 'mean'
                    else:
                        fill_str = 'neigh'
                    paras = "chr=%s_fill=%s_GC=%s_thrsh=%s_trees=%s_isl=%s" % (chroms,fill_str, gcwind, thresh, trees, island)
                    print paras

                    train_extras = np.c_[train['Exon'], train['DHS'], train['CGI']]
                    if gcwind:
                        try:
                            train_extras = np.c_[train_extras, train['GC_%s' % str(gcwind)]]
                        except ValueError:
                            train_extras = np.c_[train_extras, train['GC_100'], train['GC_400'], train['GC_1000']]
                    corrs = uf.read_corrs(path, chroms)
                    train_extras = np.c_[train_extras, corrs]

                    features_key = np.empty(38+ train_extras.shape[1],dtype='S12')
                    features_key[0] = 'Start'
                    features_key[1] = 'Nei_1_beta'
                    features_key[2] = 'Nei_1_dist'
                    features_key[3] = 'Nei_2_beta'
                    features_key[4] = 'Nei_2_dist'
                    for bnum in range(1,34):
                        features_key[bnum+4] = 'Beta' + str(bnum)
                    features_key[38] = 'Exon'
                    features_key[39] = 'DHS'
                    features_key[40] = 'CGI'
                    if gcwind:
                        features_key[41:] = 'GC'
                    features_key[-2:] = 'corrs'
                    print "train_extras[0]: %s, shape: %s" % (train_extras[1:-1][0], train_extras[1:-1].shape)
                    # Fill in NaNs in training data
                    if fill_m:
                        train_beta = fill_mean(train['Beta'])
                    else:
                        train_beta = fill_neighbors(sites, train['Beta'], NEIGH_FILL)
                # Trim first and last row off to eliminate nan values in corrs
                # Produce feature array 'X' and vector of beta values 'Y'
                # Produce feature array 'X*' to predict on and ground truth beta values 'Y*'
                    if thresh:
                        (X, Y, Xstar, Ystar) = feat_neighbors_thresh(sites[1:-1].copy(), train_beta[1:-1].copy(), \
                                            train_extras[1:-1].copy(), sample[1:-1].copy(), test[1:-1].copy(), corrs[1:-1], thresh)
                    elif island:
                        (X, Y, Xstar, Ystar) = feat_neighbors_island(sites[1:-1].copy(), train_beta[1:-1].copy(), \
                                                                        train_extras[1:-1].copy(), sample[1:-1].copy(), test[1:-1].copy(), train['CGI'])
                    else:
                        (X, Y, Xstar, Ystar) = feat_neighbors(sites[1:-1].copy(), train_beta[1:-1].copy(), \
                                            train_extras[1:-1].copy(), sample[1:-1].copy(), test[1:-1].copy())
                    print "X[0]: %s, X.shape: %s" % (np.c_[features_key, X[0]], X.shape)
                    full_feat_start = time.time()
                    # Initialize regressor with default parameters
                    model = RandomForestRegressor(oob_score=True, n_estimators=trees, n_jobs=-1)
                # Fit regressor using training data
                    model.fit(X, Y)
                # Predict on Xstar values 
                    Yhat = model.predict(Xstar)
                    print "Yhat.shape: %s" % Yhat.shape
                # Calculate r2 and RMSE
                    (r2, RMSE) = uf.calc_r2_RMSE(Yhat, Ystar)
                    importants = model.feature_importances_
    #                 print str(importants)
                    oob_str = str(model.oob_score_)
                    print "oob: " + oob_str
                    print "r2 : %f" % (r2)
                    print "RMSE: %f" % (RMSE)
                    print "RandomForest Runtime: %f" % (time.time()-full_feat_start)
                    print "Total Runtime: %f" % (time.time()-start_time)
                    paras += "_oob=%s_r2=%s_RMSE=%s_N=%s" % (oob_str, r2, RMSE, Yhat.shape[0])
                    uf.storePreds(path + 'test/', np.c_[features_key, importants], paras, full_feat_start, strs = True)

if __name__ == '__main__':
    main(sys.argv[1:])
