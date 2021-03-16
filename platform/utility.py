# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 18:16:45 2019

@author: tslin
"""

import networkx as nx
        
def errorMsg(rawStr,pos,_type,msg,n_prev=10,n_after=10,skipPos=False):
    if pos-n_prev <= 0:
        prefix = ''
        start = 0
    else:
        prefix = '...'
        start = pos-n_prev
    
    if pos+n_after >= len(rawStr)-1:
        suffix = ''
        end = len(rawStr)
    else:
        suffix = '...'
        end = pos+n_after+1
    
    if skipPos:
        posTxt = ''
    else:
        posTxt = ': at (' + str(pos) + '): '
    
    msgOut = _type + posTxt + msg + ': ' \
          + '\n\t' + prefix + rawStr[start:end] + suffix \
          + '\n\t' + ' '*len(prefix) + ' '*(pos-start) + '^'
          
    print(msgOut)
    return msgOut   

def flatten_list(l,f=lambda x:x):
    flat_list = [item for sublist in l for item in sublist]
    return flat_list

def getIdStr(idNum):
    if idNum < 10:
        return str(idNum)
    else:
        return '%'+str(idNum)

def disjoint_union(G1,G2):
    def copy_graph(G,H,n):
        for node in list(G.nodes()):
            H.add_node(n+node)
            d = G.nodes[node]
            for key in d:
                if key == 'neighList':
                    H.nodes[n+node][key] = [i+n for i in G.nodes[node][key]]
                else:
                    H.nodes[n+node][key] = G.nodes[node][key]
        
        for edge in list(G.edges()):
            
            H.add_edge(edge[0]+n,edge[1]+n)
            d = G.edges[edge]
            #newedge = ,edge[1]+n
            for key in d:
                
                if key =='direction':
                    H.edges[edge[0]+n,edge[1]+n][key] = tuple([i+n for i in G.edges[edge][key]])
                else:
                    H.edges[edge[0]+n,edge[1]+n][key] = G.edges[edge][key]
    
    
    H = nx.Graph()
    copy_graph(G1,H,0)
    copy_graph(G2,H,max(G1.nodes()))
    
    
    
    return H