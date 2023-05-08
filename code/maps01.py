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

MetaData = pd.read_csv('/media/tavo/storage/biologicalSequences/covid/datasets/MetaData.csv')
MetaData = MetaData[MetaData['correctdata']==1]
MetaData = MetaData[MetaData['qry']!='lat==0.0 & long==0.0']
MetaData = MetaData.set_index('id')
MetaData['qut'] = pd.qcut(MetaData['lengthofday'],500)

#map01
#grouped = MetaData.groupby(['date','qry'])[['Length','lat','long']]
#map02
grouped = MetaData.groupby(['qut','qry'])[['Length','lat','long']]

groupedmn = grouped.mean().droplevel(1)

mn = groupedmn['Length'].mean()
st = groupedmn['Length'].std()
groupedmn['zscore'] = [np.abs(val-mn)/st for val in groupedmn['Length']]

groupedmn = groupedmn[groupedmn['zscore']<15]

n_inx = np.unique(groupedmn.index)

cmap = cm.get_cmap('viridis')
normalizer = Normalize(groupedmn['Length'].min(),groupedmn['Length'].max())
im = cm.ScalarMappable(norm=normalizer)

for k,val in enumerate(n_inx):
    
    plt.figure(figsize=(20,10))
    axmp = plt.gca()
    
    x = groupedmn.loc[val]['long']
    y = groupedmn.loc[val]['lat']
    clrs = groupedmn.loc[val]['Length']
    
    m = Basemap(projection='cyl',llcrnrlat=-65, urcrnrlat=80,
                llcrnrlon=-180, urcrnrlon=180,ax=axmp)
    m.drawcoastlines(color='gray')
    m.fillcontinents(color='gainsboro')
    m.drawcountries(color='gray')
    
    sctr = m.scatter(x,y,s=100,c=clrs,cmap=cmap,norm=normalizer)
    x,y = m(-170,-20)
    #axmp.text(x,y,val[0:10],fontsize=16)
    
    msize = np.mean(clrs)
    x,y = m(-170,-30)
    axmp.text(x,y,'Mean Length = ' +str(int(msize)),fontsize=16)
    
    #map01
    #dateData = MetaData[MetaData['date']==val]['Length'] 
    #map02
    dateData = MetaData[MetaData['qut']==val]['Length']
    
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
    plt.savefig('/media/tavo/storage/maps/image_map_'+str(k) +'.png',dpi=75,bbox_inches='tight')
    plt.close()

