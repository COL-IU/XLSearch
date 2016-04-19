"""
Small demonstration of the hlines and vlines plots.
"""

import matplotlib.pyplot as plt
import matplotlib.text
import numpy as np
from os.path import basename
from io import StringIO
from sys import argv

file = argv[1]
base = basename(file)
seqs = base.split('_')
mz = {'alpha':[],'beta':[],'q':[]}
intensity = {'alpha':[],'beta':[],'q':[]}
annotation = {'alpha':[],'beta':[],'q':[]}
xmax = 0
ymax = 0
sequence = ''
with open(file) as f:
	#sequence = f.readline()
	for line in f:
		line = line.replace(' ','')
		lines = line.split('\t')
		mz1 = float(lines[0])
		intensity1 = float(lines[1])
		if mz1 > xmax:
			xmax = mz1
		if intensity1 > ymax:
			ymax = intensity1
		if 'Alpha' in lines[2]:
			mz['alpha'].append(mz1)
			intensity['alpha'].append(intensity1)
			annotation['alpha'].append(lines[2])
		if 'Beta' in lines[2]:
			mz['beta'].append(mz1)
			intensity['beta'].append(intensity1)
			annotation['beta'].append(lines[2])
		if '?' in lines[2]:
			mz['q'].append(mz1)
                        intensity['q'].append(intensity1)
                        annotation['q'].append(lines[2])	
xmax=1000
ymax = 100*1.5
fig = plt.figure(figsize=(30,10))
axes = fig.add_subplot(1,1,1)
plt.xlim((200,xmax))
plt.ylim((0,ymax))
axes.plot(0)
axes.vlines(mz['alpha'], [0], intensity['alpha'],colors='green',linewidths=2)
for i,txt in enumerate(annotation['alpha']) :
	if mz['alpha'][i] >= 1000:
		continue
	axes.text(mz['alpha'][i],intensity['alpha'][i],txt,color='green',fontsize=10,rotation='horizontal',family='monospace')
axes.vlines(mz['beta'],[0],intensity['beta'],colors='red',linewidths=2)
for i,txt in enumerate(annotation['beta']) :
        if mz['beta'][i] >= 1000:
                continue
        axes.text(mz['beta'][i],intensity['beta'][i],txt,color='red',fontsize=10,rotation='horizontal',family='monospace')
axes.vlines(mz['q'],[0],intensity['q'],colors='black',linewidths=2)
a = axes.text(210,90,seqs[0]+'-'+seqs[2],fontsize=20,family='monospace',color='green')
b = axes.text(210,80,seqs[1]+'-'+seqs[3],fontsize=20,family='monospace',color='red')
axes.set_xlabel('m/z')
axes.set_ylabel('Intensity')
fig.savefig(file+'.png',bbox_inches='tight')
