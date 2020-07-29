import numpy
#from torchvision import datasets, transforms
import matplotlib.pyplot as plt
import scipy.stats
import seaborn as sns


pearsons = {};
ttps = {};


acts=  ['linear','mixed']

for ncell in range(1,9):

    for seed in [1]:

        actual  = numpy.loadtxt('./data/actual-%d-%d.txt'%(ncell,seed ),  dtype=float)

        for act in acts:

            preds   = numpy.loadtxt('./data/predictions-%d-%d-%s.txt'%(ncell,seed, act ),  dtype=float)
            pearsonr = scipy.stats.pearsonr(actual, preds)
            pearsonr=numpy.power(pearsonr,2) #r squared
            tt, pp = scipy.stats.ttest_ind(actual, preds)

            if (act in pearsons):
                pearsons[act].append(pearsonr[0])
                ttps[act].append(pp)
            else:
                pearsons[act] = [pearsonr[0]] 
                ttps[act] = [pp] 

            print("predictions-%d-%d-%s : r=%f, p=%E"%(ncell,seed, act ,  pearsonr[0], pearsonr[1]))


#bars = []

#errs = []
#data = []
for act in acts:
    print ("Avg from all cells R^2 %s : %f +- %f"%( act, numpy.mean(pearsons[act]), numpy.std(pearsons[act])))
    print ("TTEST Avg  %s : %f +- %f"%( act, numpy.mean(ttps[act]), numpy.std(ttps[act])))
    #bars.append( numpy.mean(pearsons[act]))
    #errs.append( numpy.std(pearsons[act]))
    #data.append( pearsons[act])

print(pearsons)

# plt.figure()
# plt.bar(acts, bars, yerr= errs )
# plt.ylim( (min(bars) -0.2 , 1. ))

# plt.figure()
# plt.boxplot(data)
# plt.xticks(  range(1,len(acts)+1) , acts)

# plt.show()

#new nice figure 6.

a = [i for i in pearsons['linear']]
b = [i for i in pearsons['mixed']]

ma = numpy.mean(a)
mb = numpy.mean(b)
sa = numpy.std(a)
sb = numpy.std(b)
# sns.set()
fig, ax = plt.subplots(figsize=[8,8])
for i,j in zip(a,b):
    ax.plot([-1,1],[i,j], '-o',mfc='black', c='black', antialiased=True)
ax.plot([-1.5],[ma],'-bo')
ax.plot([1.5], [mb],'-ro')
plt.xticks([-1,1],['linear', 'nonlinear'], fontsize=14)
plt.xlim(-2,2)
plt.ylabel('R^2 all synapses', fontsize=14)

#plt.show() 

y1=numpy.loadtxt('./data/sixtylin.txt')
y2=numpy.loadtxt('./data/sixtynonl.txt')
my1=numpy.mean(y1)
my2=numpy.mean(y2)
fig, axs = plt.subplots(figsize=[8,8])
for l,m in zip(y1,y2):
  axs.plot([-1,1],[l,m], '-o',mfc='black', c='black', antialiased=True)
  axs.plot([-1.5], [my1],'-bo')
  axs.plot([1.5], [my2],'-ro')
  plt.xticks([-1,1],['linear', 'nonlinear'], fontsize=14)
  plt.xlim(-2,2)
  plt.ylabel('R^2 60 synapses', fontsize=14)


for cell in [2,7]:
    for fig in ['linear','mixed']:
      predslala   = numpy.loadtxt('./data/predictions-%d-1-%s.txt'%(cell, fig))
      actuallala  = numpy.loadtxt('./data/actual-%d-1.txt'%(cell))
      plt.figure()
      plt.scatter(actuallala,predslala, marker='.')
      pearsonrlala = scipy.stats.pearsonr(predslala, actuallala)
      pearsonrlala=numpy.power(pearsonrlala,2)
      ttlala, pplala = scipy.stats.ttest_ind(actuallala,predslala)
      print("%s: R^2=%f p=%f"%(fig, pearsonrlala[0],pplala))
      plt.title("%s ,  R^2=%f p=%f"%( fig,  pearsonrlala[0],pplala))
plt.show()