#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 16:13:57 2023

@author: tavo
"""

import sbmlxdf
import numpy as np

import networkx as nx
import matplotlib.pyplot as plt

from scipy.spatial import distance as ds

###############################################################################
# Loading packages 
###############################################################################

def DrawGraph(G,pos,datared,datanavy,ax):
    
    redlist,redsizes = datared
    navylist,navysizes = datanavy
    
    nx.draw_networkx_edges(G, pos,alpha=0.1,ax=ax)
    nx.draw_networkx_nodes(G, pos,nodelist=redlist,
                           node_color='red',
                           node_size=redsizes,alpha=0.5,ax=ax)
    nx.draw_networkx_nodes(G, pos,nodelist=navylist,
                           node_color='navy',
                           node_size=navysizes,alpha=0.5,ax=ax)

###############################################################################
# Loading packages 
###############################################################################

metabolism = sbmlxdf.Model('/media/tavo/storage/storage/descargas/R-HSA-1430728.xml')
metabolismdf = metabolism.to_df()

nucleotides = sbmlxdf.Model('/media/tavo/storage/storage/descargas/R-HSA-15869.xml')
nucleotidesdf = nucleotides.to_df()

###############################################################################
# Loading packages 
###############################################################################

speciesNucleotides = nucleotidesdf['species'].index
selection = [True if len(set(val.split()).intersection(set(speciesNucleotides)))>0  else False for val in metabolismdf['reactions']['reactionString']]
metabolismsubset = metabolismdf['reactions'][selection].copy()

###############################################################################
# Loading packages 
###############################################################################

adjacencylist = []

for val in metabolismsubset['reactionString']:
    splitted = val.split()
    arrowloc = [k for k,val in enumerate(splitted) if val=='=>'][0]
    
    left = [sal for sal in splitted[0:arrowloc] if len(sal)>4]
    right = [sal for sal in splitted[arrowloc::] if len(sal)>4]
    
    for xo in left:
        for yo in right:
            adjacencylist.append([xo,yo])
    
###############################################################################
# Loading packages 
###############################################################################

G = nx.Graph()
G.add_edges_from(adjacencylist)

nodesNuc = list(set(nucleotidesdf['species'].index).intersection(set(G.nodes)))
nodesMet = [val for val in G.nodes if val not in nodesNuc]

pos = nx.spring_layout(G)

plt.figure(figsize=(20,20))
ax = plt.gca()
DrawGraph(G,pos,[nodesNuc,100],[nodesMet,50],ax)

###############################################################################
# Loading packages 
###############################################################################

toRemove = []
for val in G.nodes:
    if G.degree[val]<3:# or G.degree[val]>200 :
        toRemove.append(val)
        
G.remove_nodes_from(toRemove)

nodesNuc = list(set(nucleotidesdf['species'].index).intersection(set(G.nodes)))
nodesMet = [val for val in G.nodes if val not in nodesNuc]

pos = nx.spring_layout(G)

plt.figure(figsize=(20,20))
ax = plt.gca()
DrawGraph(G,pos,[nodesNuc,100],[nodesMet,50],ax)
    
###############################################################################
# Loading packages 
###############################################################################

container = []
for val in nx.connected_components(G):
    if len(val)<100:
        container = container + list(val)

G.remove_nodes_from(container)

nodesNuc = list(set(nucleotidesdf['species'].index).intersection(set(G.nodes)))
nodesMet = [val for val in G.nodes if val not in nodesNuc]

pos = nx.spring_layout(G)

plt.figure(figsize=(20,20))
ax = plt.gca()
DrawGraph(G,pos,[nodesNuc,100],[nodesMet,50],ax)

###############################################################################
# Loading packages 
###############################################################################

centrality = nx.eigenvector_centrality_numpy(G)

sizesNuc = [1000*centrality[val] for val in nodesNuc]
sizesMet = [1000*centrality[val] for val in nodesMet]

plt.figure(figsize=(20,20))
ax = plt.gca()
DrawGraph(G,pos,[nodesNuc,sizesNuc],[nodesMet,sizesMet],ax)

threshold = 50

namesNuc = [val for val in nodesNuc if 1000*centrality[val]>threshold]
scoreNuc = [1000*centrality[val] for val in nodesNuc if 1000*centrality[val]>threshold]
orderNuc = np.argsort(scoreNuc)[::-1]
namesNuc = np.array(namesNuc)[orderNuc]
#namesNuc = namesNuc[2::]

###############################################################################
# Loading packages 
###############################################################################

G0 = G.copy()
G0.remove_nodes_from(namesNuc)

nodesNuc = list(set(nucleotidesdf['species'].index).intersection(set(G0.nodes)))
nodesMet = [val for val in G0.nodes if val not in nodesNuc]

pos = nx.spring_layout(G0)

plt.figure(figsize=(20,20))
ax = plt.gca()
DrawGraph(G0,pos,[nodesNuc,100],[nodesMet,50],ax)

###############################################################################
# Loading packages 
###############################################################################

centrality = nx.eigenvector_centrality_numpy(G)
kcentrality = nx.katz_centrality_numpy(G)
dcentrality = nx.degree_centrality(G)

centrality0 = nx.eigenvector_centrality_numpy(G0)
kcentrality0 = nx.katz_centrality_numpy(G0)
dcentrality0 = nx.degree_centrality(G0)

data = []
for val in nodesMet:
    original = [centrality[val],kcentrality[val],dcentrality[val]]
    modified = [centrality0[val],kcentrality0[val],dcentrality0[val]]
    
    data.append(ds.euclidean(original, modified))
    
data = np.array(data)
plt.figure()
plt.hist(data,bins=50)

sortorder = np.argsort(data)

selectedMet = [nodesMet[val] for val in sortorder[-50::]][::-1]
