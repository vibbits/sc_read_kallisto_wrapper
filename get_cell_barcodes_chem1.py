######### v0.2 #########
import os
import sys
json_path=os.path.abspath(sys.argv[1])
if not os.path.isfile(json_path):
    print "ERROR: Please provide path to a valid config.json file..."
    print sys.argv[1]
    exit(1)

    
workpath, jsonfile = os.path.split(json_path)
os.chdir(workpath)


import numpy as np
from itertools import islice
import math
import random
from collections import Counter
#import matplotlib.pyplot as plt
import matplotlib as mpl
import gzip, time, gc
from multiprocessing import Pool
from os import listdir
from os.path import isfile, join

import json
# use argv[1] from above 
with open(json_path) as json_file:
    parameter = json.load(json_file)

# use regex from other script 

import re

barcode_filenames = []
for f in sorted(listdir(str(parameter["BASE_DIR"]))):
    pattern = r"^"+str(parameter['sample_idx'])+"_S\d+_L00\d_([IR]\d)_001.fastq.gz"
    match = re.match(pattern, f)
    if match:
        barcode_name = match.group(1)
    if isfile(join(str(parameter["BASE_DIR"]), f)) and f.startswith(str(parameter['sample_idx'])) and barcode_name == 'R2': 
        barcode_filenames.append(f)   

print "BARCODE FILES ("+str(len(barcode_filenames))+"):"+'\n'
brc_dirs=[]
for i in range(len(barcode_filenames)):
    brc_dirs+=[str(parameter["BASE_DIR"])+str(barcode_filenames[i])]
    print brc_dirs[i]

random.seed()

print "READING FILES.."

def encoding_map(ch):
    if ch=='A':return 0
    if ch=='G':return 1
    if ch=='C':return 2
    if ch=='T':return 3
    if ch=='N':return random.randint(0,3)

decoding_lst = ['A', 'G', 'C', 'T']

def encode(k):
    code = 0
    for ch in k:
        code *= 4
        code += encoding_map(ch)
    return code

def decode(code):
    ret = ''
    for _ in range(parameter["BARCODE_LENGTH"]):
        index = code & 3
        code >>= 2
        ret = decoding_lst[index] + ret
    return ret

def read_barcodes(brc_dir):
    barcodes=[]
    with gzip.open(brc_dir) as f:
        cnt=0
        for line in f:
            if cnt%4==1: barcodes+=[encode(line[:-1])]  # remove endline character
            cnt+=1;            
    return barcodes

def hamdist(s1, s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
######################################################

p=Pool()
t0 = time.time()
barcode_vec=p.map(read_barcodes, brc_dirs )
p.close()
p.join()


#all_bars=[]
#for i in range(len(brc_dirs)):
#    all_bars+=barcode_vec[i]
#barcodes=np.array(all_bars,dtype='uint32')
barcodes = np.array([item for sublist in barcode_vec for item in sublist],dtype='uint32')

del barcode_vec[:];del barcode_vec; #del all_bars
_ = gc.collect()

t1 = time.time()
print t1-t0, "sec"

print "Barcodes:\n"
for bar in barcodes[:10]:
       print decode(bar)
print "..."
NUMBER_OF_SEQUENCED_BARCODES=len(barcodes)
print "NUMBER_OF_SEQUENCED_BARCODES =",NUMBER_OF_SEQUENCED_BARCODES

print "Detecting Cells..."

counts = Counter(barcodes)

labels, values = zip(*counts.items())

# sort your values in descending order
indSort = np.argsort(values)[::-1]

# rearrange your data
labels = np.array(labels)[indSort]
values = np.array(values)[indSort]

indices = np.arange(len(labels))

NUM_OF_DISTINCT_BARCODES=len(indices)
print "NUM_OF_DISTINCT_BARCODES =",NUM_OF_DISTINCT_BARCODES

# we apply the cellranger thresholding here - 1st percentile
# the user has provided the expected cells parameter and we take the 1st percentile
EXP_CELLS=int(math.floor(parameter["EXP_CELLS"] * 0.01 - 1))
exp_values = np.where(values >= values[EXP_CELLS] / 10)

NUM_OF_BARCODES = np.shape(exp_values)[1]

mpl.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(indices, (values))
ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')
ax.set_ylabel('UMI counts', color='k')
ax.set_xlabel('barcode', color='k')
ax.set_title('UMI counts per Barcode', color='b')
ax.axvline(NUM_OF_BARCODES, color='g', linestyle='--',linewidth=2.0)
fig.savefig(str(parameter["SAVE_DIR"])+'umi_barcodes.png')

# We do not use the original implentation to set the threshold from Pachter lab
# By default we look for a number of cells in a window of 500 to 5000. 
# WINDOW = [500,5000]
#WINDOW=parameter["WINDOW"]
#print "CELL_WINDOW:", WINDOW

#from scipy.signal import savgol_filter as savgol
#valdiff=np.diff((values))
#yhat = savgol(valdiff, 151, 1)
#NUM_OF_BARCODES=np.argmax(-yhat[WINDOW[0]:WINDOW[1]])+WINDOW[0]

# We do not plot the orginal figure from the Pachter script
#fig, ax = plt.subplots(figsize=(10,5))
#ax2=ax
#ax2.plot(indices[:-1], (-valdiff),".b")
#ax.plot(-(yhat),'r',linewidth=2.0)
#ax.set_xlim([1,10000])
#ax.set_ylim([0.01,-np.min(valdiff)])
#ax.set_yscale("log")
#ax.set_xscale("log")
#ax.set_ylabel('minus umi count diff', color='k')
#ax.set_xlabel('barcode', color='k')
#ax.set_title('minus umi count diffenrence per barcode', color='k')
#ax.axvline(NUM_OF_BARCODES, color='g', linestyle='--',linewidth=2.0)
#fig.savefig(str(parameter["SAVE_DIR"])+'count_diff_umi_barcodes.png')

print "Cell_barcodes_detected:",NUM_OF_BARCODES

#NUM_OF_BARCODES=np.argmax(-yhat[WINDOW[0]:WINDOW[1]])+WINDOW[0]
#print "Cell_barcodes_detected:",NUM_OF_BARCODES

NUM_OF_READS_in_CELL_BARCODES = sum(values[:NUM_OF_BARCODES])
print "NUM_OF_READS_in_CELL_BARCODES =",NUM_OF_READS_in_CELL_BARCODES

codewords=labels[:NUM_OF_BARCODES]

print "Calculating d_min..."

Ham_dist=np.zeros([len(codewords),len(codewords)])
for i in range(len(codewords)):
    codi=decode(codewords[i])
    for j in range(i+1,len(codewords)):
        Ham_dist[i,j]=hamdist(codi,decode(codewords[j]))
        Ham_dist[j,i]=Ham_dist[i,j]
dmin=(Ham_dist+parameter["BARCODE_LENGTH"]*np.identity(len(codewords))).min(axis=1)
### to be on the safe side we suggest to correct only barcodes that have d_min>=4
d=parameter['dmin']
brc_idx_to_correct=np.arange(len(codewords))[dmin>=d]
print "number of cell barcodes to error-correct:", len(brc_idx_to_correct), "( dmin >=", d,")"

## CLEANUP
#del indices; del labels; del values; del counts; del valdiff; del indSort; del Ham_dist; del dmin
del indices; del labels; del values; del counts; del indSort; del Ham_dist; del dmin
_ = gc.collect()

print "Writing output..."

import pickle

save_dir=str(parameter["SAVE_DIR"])
#create output directory 

import os
if not os.path.isdir(save_dir):
    try:
        os.mkdir(save_dir)
    except OSError as e:
        print "OSError({0}): {1}".format(e.errno, e.strerror)

with open(save_dir+"barcodes.dat", 'wb') as f:
    pickle.dump(barcodes,f)
with open(save_dir+"codewords.dat", 'wb') as f:
    pickle.dump(codewords,f)
with open(save_dir+"brc_idx_to_correct.dat", 'wb') as f:
    pickle.dump(brc_idx_to_correct,f)
printer=""
printer+="NUMBER_OF_SEQUENCED_BARCODES: %s\n" % NUMBER_OF_SEQUENCED_BARCODES
printer+="NUM_OF_DISTINCT_BARCODES: %s\n" % NUM_OF_DISTINCT_BARCODES   
printer+="Cell_barcodes_detected: %s\n" % NUM_OF_BARCODES
printer+="NUM_OF_READS_in_CELL_BARCODES: %s\n" % NUM_OF_READS_in_CELL_BARCODES
printer+="NUM_OF_CELL_BARCODES_to_CORRECT %s (dmin >=%s)\n" % (len(brc_idx_to_correct), d)    
with open(save_dir+"run.info", 'wb') as f:
    f.write(printer)
print '\n'
print printer
print "DONE."
