"""

This file loads the data generated from the simulations and 
generates the engram-related figures.

usage:
python make_figs.py


"""
import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
import sys
import scipy.stats as stats

import matplotlib as mpl
mpl.rcParams["errorbar.capsize"] = 2
mpl.rcParams["lines.linewidth"] = 1
mpl.rcParams['pdf.fonttype'] = 42


np.set_printoptions( threshold=999999999999999)


TESTSIZE=10000
# Estimate overlap chance level by random shuffling two populations
def shuffletest(asz, bsz):
	a = np.zeros(TESTSIZE)
	b = np.zeros(TESTSIZE)

	a[0: int(asz*TESTSIZE)] = 1
	b[0: int(bsz*TESTSIZE)] = 1

	s =0
	for i in xrange(1000):
		np.random.shuffle(a)
		np.random.shuffle(b)
		s += np.sum(np.multiply(a,b)) / float(TESTSIZE)

	return s/1000.




def trevrolls(frates):
	r2=0.
	rs =0.
	n = float(frates.size)
	for i in range(frates.size): # np.nditer(frates):
		r2 += (frates[i]**2)/n
		rs += frates[i]/n
	return 1. - ((rs**2)/r2)
	

def loadspikesdat(filename, tduration):

	ff = open(filename, 'r') 
	fdata = ff.readlines()
	sx = len(fdata)
	sy = tduration;
	raster = np.zeros( (sx, sy) );
	nid=0
	for l in fdata:
		ar = np.fromstring(l, sep=' ' , dtype=int)
		raster[nid, ar] = 1
		raster[nid,0] =0 # XXX bug
		nid += 1

	return raster


def printpairstats(stat, name):
	print( name)
	print( stats.f_oneway( stat[0][0] , stat[1][0] ) ) 
	print( stats.f_oneway( stat[0][1] , stat[1][1] ) )
	print( stats.f_oneway( stat[0][0] , stat[0][1] ) )
	print( stats.f_oneway( stat[1][0] , stat[1][1] ) )

def exportcsv(data, name):
	mycsv = np.array( [data[0,0] , data[1,0], data[0, 1], data[1,1] ]);
	df = pd.DataFrame(mycsv.T)
	df.to_csv(name+".txt")

def label_diff(ax, i,j,text,X,Y):
    x = (X[i]+X[j])/2
    y = 1.1*max(Y[i], Y[j])
    dx = abs(X[i]-X[j])
    props = {'connectionstyle':'bar','arrowstyle':'-',\
                 'shrinkA':20,'shrinkB':20,'linewidth':1}
    ax.annotate(text, xy=(X[i],y-7), zorder=10, transform=ax.transData)
    #ax.text(.5, .5, "text")
    ax.annotate('', xy=(X[i],y), xytext=(X[j],y), arrowprops=props)


def makebps(pops):

	bp = plt.boxplot(  (pops[0], pops[1], pops[2]), notch=True, patch_artist=True)
	for box,color in zip( bp['boxes'], box_colors):
		box.set(color=color)
	plt.xticks( [1, 2, 3], ['Linear', 'Nonlinear\nDispersed', 'Nonlinear\nIn Branch'])

	print("[Mean , STD]")
	for i in range(3):
		print np.mean(pops[i]), np.std(pops[i])



NPYRS = 400
NRUNS=20


prp='G'
box_colors = [[189/256.,34/256.,28/256.] , [165/256.,186/256.,224/256.], [227/256.,195/256.,221/256.] ];

XLEN = 3
trs = np.zeros(( XLEN,NRUNS))
ffs = np.zeros(( XLEN,NRUNS));
pops = np.zeros(( XLEN,NRUNS));

RES_SET='single'

plot =0


plt.figure(); 
for somType in  [0, 2, 3]:

	for idx in  [0,1,2]:
	
		if (idx ==0):  #linear bar
			clustered = 0
			pvType = 2
		elif (idx ==1): #  bimodal dispersed bar
			clustered=0
			pvType=3
		elif (idx ==2): # bimodal clustered
			clustered=1
			pvType =3


		branch_hist = [];
		for run in range(NRUNS):

			spikes = np.loadtxt( '../data/%s_%d_%d_%d_%s_%d/spikesperpattern.dat'%(RES_SET, pvType, somType, clustered, prp, run), dtype=float)
			spikes = spikes[ 0:NPYRS]
			tr=  trevrolls(spikes)

			trs[idx, run] = tr

			ff = np.mean(spikes)
			ffs[idx, run] = ff/4.

			pops[idx, run] = 100.* sum(spikes>=40.) / float(NPYRS)


	plot += 1
	plt.subplot(3,3,plot)

	makebps(pops)
	plt.ylim(ymin=5, ymax=50)
	#exportcsv(pops, 'supl_pops')
	plt.ylabel('Engram Size (%) ');


	plot += 1
	plt.subplot(3,3,plot)
	makebps(ffs)
	#exportcsv(ffs, 'supl_firing')
	plt.ylabel('Mean Firing Rate (Hz) ');


	plot += 1
	plt.subplot(3,3,plot)
	makebps(trs)
	#printpairstats(trs, 'Sparsity')

	plt.ylim(ymin=0.4, ymax=1.)
	plt.ylabel('Sparsity');

	#plt.savefig("sup_%s.pdf"%(SOMCONDITION))



plt.show();



