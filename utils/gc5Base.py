'''
Created on May 7, 2015

@author: jonathanshor
'''
import sys
import numpy as np
from optparse import OptionParser
import UtilityFunctions
# from __builtin__ import enumerate

DIGITS = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
OFFSETS = {'1':0, '2':560416768,'3':1154476032,'4':1636007936,'5':2098747392, '6':2536030208, '7':2946803712, '8':3695589376, '9':4043157504, '10':4335620096, '11':4653656064, '12':4971237376, '13':5286854656, '14':5517381632, '15':5729316864, '16':5924708352, '17':6110967808, '18':6294654976, '19':6669289472, '20':6470326272, '21':6884116480, '22':6800637952}
POS_MIN = 10000

# def gc5Base_intersect_to_bed(gcpathname, path, fnmask, window=400):
# gcpathname is the full path and filename of the gc5Base file 
# path is for the location of chromosonal bed files containing sites for which to pull gc data
# fnmask is a string containing 'chr' that will be used to determine bed filenames by replacing the 'chr' with 'chr1', 'chr2',...
def gc5Base_intersect_to_bed(gcpathname, sites, window=400):
    gc_prop = np.zeros(len(sites))
    with open(gcpathname) as gcf:
        cur_chrom = -1
        last_tell = 0
        for i, r in enumerate(sites):
            #Ensure no more than a chrom away in gcf
            chrom = r[0][3:]
            if chrom != cur_chrom:
                cur_chrom = chrom
                gcf.seek(OFFSETS[str(chrom)])
                gcf.readline()
                last_tell = gcf.tell()
            else:
                gcf.seek(last_tell)
            #Determine gc sites of relevance
            site_min = max(POS_MIN,int(r[1]) - window)
            site_max = int(r[2]) + window
            cur_window = 0
            #Get gc_prop val
            for line in gcf:
                if line[0] in DIGITS:
                    delim = line.find('\t')
                    cur_pos = int(line[:delim])
                    if cur_pos + 4 >= site_min:
                        if cur_pos >= site_min and cur_pos + 4 <= site_max:
                            gc_prop[i] += float(line[delim+1:-1])
                            cur_window += 5
                        elif cur_pos > site_max:
                            break
                        elif cur_pos < site_min:
                            gc_prop[i] += (float(line[delim+1:-1]) * (5-(site_min - cur_pos)))/ 5
                            cur_window += 5 - (site_min - cur_pos)
                        else:
                            gc_prop[i] += (float(line[delim+1:-1]) * (5-(site_max - cur_pos)))/ 5
                            cur_window += 5 - (site_max - cur_pos)
                    else:
                        last_tell += len(line)
                else:
                    break
            if cur_window == 0:
                gc_prop[i] = 0
            else:
                gc_prop[i] /= float(cur_window)
    gcf.close()
    GCs = np.array(gc_prop, dtype='S10')
    return GCs
        

def gc5Base_get_chrom_offsets(gcpathname):
# Return Nx2 array: col 0 = chrom ID; col 1 = tell() offset
    gc_offs = np.empty((0,2),dtype='string')
    with open(gcpathname) as gcf:
        curtell = gcf.tell()
        for line in gcf:
            if line[0] not in DIGITS:
                #line.find('chrom=chr')+len('chrom=chr') == 22
                chromID = line[22:line.find(' ',22)]
                np.append(gc_offs, [[chromID,str(curtell)]], 0)
                print "ChromID: %s, offset: %s" % (chromID, str(curtell))
            curtell = gcf.tell()
    gcf.close()
    return gc_offs

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", help='read bed data from PATH', metavar='PATH')
    parser.add_option("-g", dest="gcpathname", help='gc5Base file')
    parser.add_option("-o", dest="getoffsets", action="store_true", default=False)
    (options, _args) = parser.parse_args()     
    path = options.path
    print "PATH = " + path
    
    if options.gcpathname:
        for chrom in range(1,22):
            train = UtilityFunctions.read_bed_dat_train(path, chrom)
            GCs = np.empty((len(train),4),dtype='S10')
            GCs[:,0] = train['Chrom']
            GCs[:,1] = train['Start']
            GCs[:,2] = train['End']
            sites = GCs[:,:3]
    #         sites = np.array((train['Chrom'],train['Start'],train['End']))
            window = 400
            GCs[:,3] = gc5Base_intersect_to_bed(options.gcpathname, sites, window = window)
            np.savetxt(path + 'gc_chr%s_w=%s.bed' % (str(chrom), str(window)), GCs, fmt='%s')
    
    if options.gcpathname and options.getoffsets:
        gc_offsets = gc5Base_get_chrom_offsets(options.gcpathname)
        np.savetxt(path + 'gc_offsets.txt', gc_offsets)

if __name__ == '__main__':
    main(sys.argv[1:])