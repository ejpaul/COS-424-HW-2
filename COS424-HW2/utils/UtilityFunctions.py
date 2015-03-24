'''
Created on Mar 23, 2015

@author: jonathanshor
'''
import nltk, numpy
import getopt, sys, time

#BE SURE TO UPDATE TO PATH TO YOUR bed FILES
YOUR_PATH = "/Users/jonathanshor/Google Drive/Edu/2015 S/COS424 - Interacting with Data/COS424-HW2/base_files/"


TRAIN_FNAME = "intersected_final_chr1_cutoff_20_train_revised.bed"
SAMPLE_FNAME = "intersected_final_chr1_cutoff_20_sample.bed"
TEST_FNAME = "intersected_final_chr1_cutoff_20_test.bed"

def read_bed_dat(myfile):
    return numpy.loadtxt(myfile, dtype="string", delimiter="\t")
    

#Currently does not acknowledge command line args; set all paths via CONSTANTS above
def main(argv):
#     path = ''
#     
#     try:
#         opts, args = getopt.getopt(argv, "p:", ["path="])
#     except getopt.GetoptError:
#         print 'Bad args'
#         sys.exit(2)
#     for opt, arg in opts:
#         if opt in ('-p', '--path'):
#             path = arg
    path = YOUR_PATH
    
    trains = read_bed_dat(path + TRAIN_FNAME)
    sample = read_bed_dat(path + SAMPLE_FNAME)
    
    print "Sample[0] =%s" % sample[0]
    print "Trains[0] =%s" % trains[0]
    print "len(Trains) =%s" % len(trains)
if __name__ == '__main__':
    main(sys.argv[1:])
#    main("-p /Users/jonathanshor/Google Drive/Edu/2015 S/COS424 - Interacting with Data/COS424-HW2/base_files/")