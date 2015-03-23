import nltk, numpy
import getopt, sys, time

def main(argv):

    try:
        opts, args = getopt.getopt(argv, "p:", ["path="])
    except getopt.GetoptError:
        print 'python COS424-HW2-Util.py -p <path>'
        sys.exit(2)
    for opt, arg in opts:
        




if __name__ == "__main__":
    main(sys.argv[1:])