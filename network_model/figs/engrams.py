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
	bp = plt.boxplot(  (pops[0][0], pops[0][1], pops[1][1]), notch=True, patch_artist=True)
	for box,color in zip( bp['boxes'], box_colors):
		box.set(color=color)
	plt.xticks( [1, 2, 3], ['Linear', 'Nonlinear\nDispersed', 'Nonlinear\nIn Branch'])

	print("[Mean , STD]")
	print np.mean(pops[0][0]), np.std(pops[0][0]) 
	print np.mean(pops[0][1]), np.std(pops[0][1]) 
	print np.mean(pops[1][1]), np.std(pops[1][1]) 


def plotCase(case, title):
	NRUNS=10
	for prp in ['G' ]:
		print "Case=",case

		dend_ids = [0,1]
		dend_conds = [2,3]
		dend_ticks = ['linear', 'nonlinear'];

		XLEN = len(dend_conds);

		trs = np.zeros((2, XLEN,NRUNS))


		for CLUSTERED in [0,1]:
			for did in dend_ids:
				dend_cond = dend_conds[did]

				branch_hist = [];
				for run in range(NRUNS):
					spikes = np.loadtxt( '../data/%s_%d_%d_%s_%d/spikesperpattern.dat'%(case, dend_cond, CLUSTERED, prp, run), dtype=float)
					spikes = spikes[:, 0:NPYRS]

					overlap = 100.*np.sum(np.logical_and((spikes[0,:] >40.),  ( spikes[1,:] >40.) )) / NPYRS

					trs[CLUSTERED, did, run] = overlap;
					"""
					syns = np.loadtxt('./data/two_%d_%d_%s_%d/syn-post.txt'%(dend_cond, CLUSTERED, prp, run), dtype=float)
					syns = syns[syns[:,4]>0.7]
					cols = ['input_id', 'group_id', 'branch_id','nid','weight']
					table = pd.DataFrame(syns,columns=cols)
					totals = table.groupby(['branch_id'])['weight'].count().values
					#branch_hist.extend(totals.tolist())

					ratio_clustered[did, run] = sum(totals>2) / float(sum(totals>0))
					"""

				#plt.figure()
				#plt.title('Synapses per branch cond=%d'%(dend_cond));
				#df = pd.DataFrame(np.array(branch_hist), columns=['cnt'])
				#df['cnt'].value_counts().plot(kind='bar')
				#plt.ylim(ymin=0)
				#plt.xticks(np.arange(4), ['supra', 'sub', 'linear', 'mixed']);	





		xr = np.arange(XLEN)

		#plt.title('case=%s'%(case));
		#trs = trs[2:3, :];
		#print( stats.f_oneway( trs[0,:] , trs[1,:] ) )

		tr_means = np.mean(trs, 1)
		tr_std =  np.std(trs,1)
		printpairstats(trs, title);
		exportcsv(trs, case)

		#plt.bar(xr-.2, np.mean(trs[0], 1), yerr=np.std(trs[0],1), width=.3 )
		#plt.bar(xr+.2, np.mean(trs[1], 1), yerr=np.std(trs[1],1), width=.3 )

		makebps(trs)

		plt.ylim(ymin=0, ymax=50);
		#plt.xticks(np.arange(XLEN), dend_ticks);
		plt.ylabel('Overlap (%)');
		#plt.title(title)

NPYRS = 400
NRUNS=20
CLUSTERED=0


prp='G'
dend_ids = [0,1]
dend_conds = [2,3]
dend_ticks = ['Linear', 'Nonlinear'];
box_colors = [[189/256.,34/256.,28/256.] , [165/256.,186/256.,224/256.], [227/256.,195/256.,221/256.] ];

XLEN = len(dend_conds);

trs = np.zeros((2, XLEN,NRUNS))
ffs = np.zeros((2, XLEN,NRUNS));
pops = np.zeros((2, XLEN,NRUNS));
ratio_clustered = np.zeros((2, XLEN,NRUNS));

RESULTSET='single'
SOMTWO='two'

for CLUSTERED in [0,1]:
	print "PRP=",prp


	for did in dend_ids:
		dend_cond = dend_conds[did]

		branch_hist = [];
		for run in range(NRUNS):
			spikes = np.loadtxt( '../data/%s_%d_%d_%d_%s_%d/spikesperpattern.dat'%(RESULTSET, dend_cond, 1, CLUSTERED, prp, run), dtype=float)
			spikes = spikes[ 0:NPYRS]
			tr=  trevrolls(spikes)
			trs[CLUSTERED, did, run] = tr;
			ff = np.mean(spikes);
			ffs[CLUSTERED, did, run] = ff/4.;
			pops[CLUSTERED, did, run] = 100.* sum(spikes>=40.) / float(NPYRS)

			syns = np.loadtxt('../data/%s_%d_%d_%d_%s_%d/syn-post.txt'%(RESULTSET, dend_cond, 1, CLUSTERED, prp, run), dtype=float)
			syns = syns[syns[:,4]>0.7]
			cols = ['input_id', 'group_id', 'branch_id','nid','weight']
			table = pd.DataFrame(syns,columns=cols)
			totals = table.groupby(['branch_id'])['weight'].count().values
			branch_hist.extend(totals.tolist())

			ratio_clustered[CLUSTERED, did, run] = sum(totals>2) / float(sum(totals>0))

		#plt.figure()
		#plt.title('Synapses per branch cond=%d'%(dend_cond));
		#df = pd.DataFrame(np.array(branch_hist), columns=['cnt'])
		#df['cnt'].value_counts().plot(kind='bar')
		#plt.ylim(ymin=0)
		#plt.xticks(np.arange(4), ['supra', 'sub', 'linear', 'mixed']);	




#plt.figure(); #figsize=(8,6))
#plt.subplot(3,2,1)
#ratio_clustered *= 100.
##plt.bar(xr-.2, np.mean(ratio_clustered[0], 1), yerr=np.std(ratio_clustered[0],1), width=.3 )
##plt.bar(xr+.2, np.mean(ratio_clustered[1], 1), yerr=np.std(ratio_clustered[1],1), width=.3 )
#makebps(ratio_clustered)
#printpairstats(ratio_clustered, 'Clustered')
#exportcsv(ratio_clustered, 'clustered')
#plt.ylim(ymin=10)
#plt.ylabel('Clustered Engram Synapses (%)');



plt.subplot(2,3,1)
#xr = np.arange(XLEN)
printpairstats(pops, 'Populations')
makebps(pops)
plt.ylim(ymin=5, ymax=60)
#plt.legend(['Dispersed', 'Clustered'])
exportcsv(pops, 'pops')
plt.ylabel('Engram Size (%) ');

popsize =  np.mean(pops[0][0])/100.
print("Random  overlap linear for %f  = %f"%(popsize, shuffletest(popsize, popsize)))
print( shuffletest(popsize, popsize))

popsize =  np.mean(pops[0][1])/100.
print("Random  overlap nonl/disp for %f  = %f"%(popsize, shuffletest(popsize, popsize)))

popsize =  np.mean(pops[1][1])/100.
print("Random  overlap nonl/inbr for %f  = %f"%(popsize, shuffletest(popsize, popsize)))




print np.mean(pops[0][1]), np.std(pops[0][1]) 
print np.mean(pops[1][1]), np.std(pops[1][1]) 


plt.subplot(2,3,2)

#ffs = ffs[2:3,:];
#plt.bar(xr-.2, np.mean(ffs[0], 1), yerr=np.std(ffs[0],1), width=.3 )
#plt.bar(xr+.2, np.mean(ffs[1], 1), yerr=np.std(ffs[1],1), width=.3 )

makebps(ffs)
printpairstats(ffs, 'Firing')
exportcsv(ffs, 'firing')

plt.ylim(ymin=5)
plt.ylabel('Mean Firing Rate (Hz) ');





plt.subplot(2,3,3)

#plt.bar(xr-.2, np.mean(trs[0], 1), yerr=np.std(trs[0],1), width=.3 )
#plt.bar(xr+.2, np.mean(trs[1], 1), yerr=np.std(trs[1],1), width=.3 )
makebps(trs)
printpairstats(trs, 'Sparsity')
exportcsv(trs, 'sparsity')
plt.ylim(ymin=0., ymax=1.)
plt.ylabel('Sparsity');



plt.subplot(2,3,4)
plotCase(SOMTWO, "1 hour separation")

plt.subplot(2,3,5)
plotCase(SOMTWO+"1440", "24h separation")

plt.savefig("bp_%s.pdf"%(RESULTSET))

plt.show();

