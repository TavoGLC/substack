#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 23:00:31 2023

@author: tavo
"""

import numpy as np
import pandas as pd 
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.colors import Normalize

from mpl_toolkits.basemap import Basemap

from scipy import interpolate
from itertools import product

###############################################################################
# Loading packages 
###############################################################################

def BottomStyle(Axes): 
    """
    Parameters
    ----------
    Axes : Matplotlib axes object
        Applies a general style to the matplotlib object

    Returns
    -------
    None.
    """    
    Axes.spines['top'].set_visible(False)
    Axes.spines['bottom'].set_visible(True)
    Axes.spines['left'].set_visible(False)
    Axes.spines['right'].set_visible(False)
    Axes.xaxis.set_tick_params(labelsize=8)
    Axes.yaxis.set_ticks([])
    
###############################################################################
# Loading packages 
###############################################################################

def MakeSingleStringSeq(xxc,yyc,data,clr,axes,swap=False):
    
    y = np.array(data)
    y = (y - y.min())/ (y.max()-y.min())
    y = list(y)
    x = np.linspace(0,1,num=int(len(y)/2)+1)
    x = list(x[1::][::-1])+list(x[1::])

    tck, u = interpolate.splprep([x + x[:1], y + y[:1]], s=0, per=True)
    unew = np.linspace(0, 1, 100)
    
    basic_form = interpolate.splev(unew, tck)
    xm,ym = np.mean(basic_form[0]),np.mean(basic_form[1])
    
    xcoords,ycoords = scale*(basic_form[0]-xm)+xxc, scale*(basic_form[1]-ym)+yyc
    if swap:
        xcr = ycoords
        ycr = xcoords
    else:
        xcr = xcoords
        ycr = ycoords
        
    axes.plot(xcr, ycr, color=clr, lw=1)
    axes.fill(xcr, ycr, color=clr, alpha=0.55)
    
    return xcoords,ycoords

###############################################################################
# Loading packages 
###############################################################################

Alphabet = ['A','C','T','G']
Blocks = []

maxSize = 5
for k in range(1,maxSize):
    
    Blocks.append([''.join(i) for i in product(Alphabet, repeat = k)])

###############################################################################
# Loading packages 
###############################################################################

MetaData = pd.read_csv('/media/tavo/storage/biologicalSequences/covid/datasets/MetaData.csv')
MetaData = MetaData[MetaData['correctdata']==1]
MetaData = MetaData[MetaData['qry']!='lat==0.0 & long==0.0']
MetaData = MetaData.set_index('id')
MetaData['qut'] = pd.qcut(MetaData['lengthofday'],500)

dataspots = pd.read_csv('/media/tavo/storage/sunspots/sunspots.csv')
dataspots['date'] = pd.to_datetime(dataspots['date'])
rollingavgspots = dataspots.groupby('date')['dailysunspots'].mean()
MetaData['spots'] = np.array(rollingavgspots.loc[MetaData['date']])


KmerData = pd.read_csv('/media/tavo/storage/biologicalSequences/covid/datasets/KmerDataUpd.csv')
KmerData = KmerData.set_index('id')
KmerData = KmerData.loc[MetaData.index]

data = pd.concat([MetaData,KmerData],axis=1)

#map01 date map02 qut
col = 'date'
grouped = data.groupby([col,'qry'])[['Length','lat','long']+Blocks[0]+Blocks[1]+Blocks[2]+Blocks[3]]

groupedmn = grouped.mean().droplevel(1)

mn = groupedmn['Length'].mean()
st = groupedmn['Length'].std()
groupedmn['zscore'] = [np.abs(val-mn)/st for val in groupedmn['Length']]

groupedmn = groupedmn[groupedmn['zscore']<15]

n_inx = np.unique(groupedmn.index)

cmap = cm.get_cmap('viridis')
normalizer = Normalize(groupedmn['Length'].min(),groupedmn['Length'].max())
im = cm.ScalarMappable(norm=normalizer)
shapeColumns = Blocks[0]
scale = 7.5

bins = pd.interval_range(start=-65, end=80,freq=15)

for k,val in enumerate(n_inx):
    
    plt.figure(figsize=(20,10))
    axmp = plt.gca()
    
    daydata = groupedmn.loc[val]
    
    x = np.array(daydata['long']).reshape(-1)
    y = np.array(daydata['lat']).reshape(-1)
    clrs = np.array(daydata['Length']).reshape(-1)
    
    m = Basemap(projection='cyl',llcrnrlat=-65, urcrnrlat=80,
                llcrnrlon=-180, urcrnrlon=180,ax=axmp)
    m.drawcoastlines(color='gray')
    m.fillcontinents(color='gainsboro')
    m.drawcountries(color='gray')
    
    xd,yd = m(-170,15)
    axmp.text(xd,yd,val[0:10],fontsize=16)
    #axmp.text(xd,yd,'Sunspots = ' + str(val),fontsize=16)
    
    msize = np.mean(clrs)
    xs,ys = m(-170,5)
    axmp.text(xs,ys,'Shape equal to mean composition by latitude',fontsize=10)
    
    shapeData = daydata[shapeColumns]
    disc = shapeData.values.shape
    
    if disc==(len(shapeColumns),):
        shapeData = shapeData.values.reshape(1,-1)
    else:
        shapeData = shapeData.values
    
    for ii,(xx,yy) in enumerate(zip(x,y)):
        
        xcrd,ycrds = MakeSingleStringSeq(xx,yy,
                                         shapeData[ii],
                                         cmap(normalizer(clrs[ii])),
                                         axmp)
    
    staxes = axmp.inset_axes([0.015,0.25,0.25,0.20])
    
    if daydata.ndim==2:
        
        latsdata = daydata.groupby(pd.cut(daydata['lat'],bins))[['Length']+shapeColumns].mean()
        mid = [(val.left+val.right)/2 for val in latsdata.index]
        
        for kl,ltr in enumerate(latsdata.index):
            
            latdta = latsdata[shapeColumns].loc[ltr]
            latclr = latsdata['Length'].loc[ltr]
            
            if np.isnan(latclr)==False:
                
                xcrd,ycrds = MakeSingleStringSeq(0,mid[kl],latdta,
                                                 cmap(normalizer(latclr)),staxes,swap=True)
                
    staxes.set_xlim([-70,85])
    BottomStyle(staxes)
        
    dateData = MetaData[MetaData[col]==val]['Length'] 
    
    a,b = np.histogram(np.array(dateData),bins=100)
    b = [(b[j]+b[j+1])/2 for j in range(len(b)-1)]
    a = (a - a.min())/(a.max() - a.min())
    
    haxes = axmp.inset_axes([0.015,0.15,0.25,0.05])
    haxes.plot(b,a,color='black')
    haxes.set_xlim([groupedmn['Length'].min(),groupedmn['Length'].max()])
    haxes.axis('off')
    
    cbaxes = axmp.inset_axes([0.015,0.1,0.25,0.025])
    cbar = plt.colorbar(im,cax=cbaxes, orientation='horizontal')
    cbar.set_label('Genome Length',fontsize=12)
    
    axmp.axis('off')
    plt.savefig('/media/tavo/storage/maps/image_map_'+str(k) +'.png',dpi=150,bbox_inches='tight')
    plt.close()
