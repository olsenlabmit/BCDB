# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 12:15:46 2019

@author: tslin, nrebello
"""

from SMILES import SMILES
from BigSmilesPattern import BigSmilesPattern
from BigSMILES_Bond import BigSMILES_Bond
from BigSMILES_StoObj import BigSMILES_StoObj
from utility import errorMsg, flatten_list
from error import BigSMILESError,BigSMILES_BondInconsistencyError, BigSMILES_StoObjMissingTerminalError, BigSMILES_StoObjError
import networkx as nx
from collections import deque

class BigSMILES(SMILES):
    _count = 0
    def __init__(self,inStr="",pos=0,UpperBond_Dict=None,index=list()):
        BigSMILES._count = BigSMILES._count + 1
        self.rawStr = inStr.split()[0]
        self.index = index
        self.dummyStoObjStr = None
        self.pos = pos
        self.parsed = False
        self.noWriteBondDesc=False
        pattern = BigSmilesPattern._BigSmilesElement
        #self.Bond_Dict = Bond_Dict
        # UpperBond_Dict is the dictionary that stores the bonding descriptors that had already been declared within a stochastic object
        # it is used to check if inconsistency occur between distinct repeat units or end groups
        if UpperBond_Dict:
            self.Bond_Dict = dict()
            for key in UpperBond_Dict.keys():
                self.Bond_Dict[key] = UpperBond_Dict[key].copy()
        else:
            self.Bond_Dict = dict() # dictionary for storing bond
        self.StoObj_List = list() # list for storing Stochastic Obj
        
        super(BigSMILES,self).__init__(inStr=self.rawStr,pattern=pattern)
        # TODO! uncomment this and test if recursive parse works!!!
        self.parse()
        
    def __str__(self):
        s = self.writeStandard()
        return s
    
    def __len__(self):
        return len(self.StoObj_List)
    
    def __getitem__(self,key):
        if key <= self.__len__() and key >= 0:
#            return self.StoObj_List[key-1]
            return self.StoObj_List[key]
        else:
            return None
    
    def __iter__(self):
        self.n = 0
        return self
    
    def __next__(self):
        if self.n < len(self.StoObj_List):
            result = self.StoObj_List[self.n]
            self.n += 1
            return result
        else:
            raise StopIteration
    
    
    def addBigSmilesBondAtom(self,res,prevAtom,pos,prevBond):
        if 'BigSMILES_Bond' in res.keys():
            _type = 'BigSMILES_Bond'
#        else:
#            _type = 'BigSMILES_ladderBond'
 
        # add the BigSMILES bonding descriptor as if it is an atom
        self.atomCount = self.atomCount + 1
        nodeId = self.atomCount
        self.G.add_node(nodeId)
        self.G.nodes[nodeId]['rawStr'] = res.rawStr
        self.G.nodes[nodeId]['_type'] = _type
        self.G.nodes[nodeId]['neighList'] = list()
        self.G.nodes[nodeId]['pos'] = pos
        self.G.nodes[nodeId]['atom'] = _type
        self.G.nodes[nodeId]['chiral'] = ''
        
        isFirst = False
        if prevAtom == None:
            isFirst = True
        
        
        if isFirst:
            S_bond = 'u'
        else:
            if prevBond == None:
                prevBond = '-'
            S_bond = prevBond
        try:    
            self.G.nodes[nodeId]['BigSMILES_Bond'] = BigSMILES_Bond(res,S_bond,Bond_Dict=self.Bond_Dict,item=nodeId)
        except:
            errorMsg(self.rawStr,pos,'Error','Inconsistency between bonding descriptor')
            raise BigSMILES_BondInconsistencyError
        
        currentAtom = nodeId
        if prevAtom != None:
            self.createEdge(prevAtom,currentAtom,prevBond)
        prevBond = None
        return currentAtom,prevBond


            
    def addBigSmilesStoObjAtom(self,res,prevAtom,pos,prevBond):
        # add the BigSMILES dummy stochastic object descriptor as if it is an atom
        _type = 'BigSMILES_StoObj'
        self.atomCount = self.atomCount + 1
        nodeId = self.atomCount
        self.G.add_node(nodeId)
        self.G.nodes[nodeId]['rawStr'] = res.rawStr
        self.G.nodes[nodeId]['_type'] = _type
        self.G.nodes[nodeId]['neighList'] = list()
        self.G.nodes[nodeId]['pos'] = pos
        self.G.nodes[nodeId]['atom'] = _type
        self.G.nodes[nodeId]['chiral'] = ''
        
        
        # TODO! parse stochastic obj
        #print(pos)
        StoObjIdx = len(self.StoObj_List)+1
        try:
            self.G.nodes[nodeId]['BigSMILES_StoObj'] = BigSMILES_StoObj(res.rawStr,pos=pos,index=self.index+[StoObjIdx],prevBond=prevBond)
        except:
            #errorMsg(self.rawStr,pos,'Error','Inconsistency between bonding descriptor')
            raise BigSMILES_StoObjError  
        else:
            self.StoObj_List.append(self.G.nodes[nodeId]['BigSMILES_StoObj'])
            
        self.G.nodes[nodeId]['StoObjId'] = len(self.StoObj_List)
        
        left = self.G.nodes[nodeId]['BigSMILES_StoObj'].getBond(end='left')
        right = self.G.nodes[nodeId]['BigSMILES_StoObj'].getBond(end='right')
        #left = '='
        #right = None
        
        # check if the left terminal bond is specified
        if left == None and prevAtom != None:
            raise BigSMILES_StoObjMissingTerminalError
        
        if prevBond != None: 
            # raise error if the two terminal bonding descriptors found on adjacent BigSMILES stochastic objects are inconsistent
            if left != prevBond:
                raise BigSMILES_BondInconsistencyError
        else:
            prevBond = left
        

        currentAtom = nodeId
        if prevAtom == None: # open bonding site, add an explicit H atom
            if prevBond == None:
                pass
            else:
                if prevBond != '' and prevBond != '-':
                    tmpHAtom = self.createNode('[H]','organic','H',pos,isotope='1')
                else:
                    tmpHAtom = self.createNode('[H]','organic','H',pos)
                self.createEdge(tmpHAtom,currentAtom,prevBond)
        else:
            self.createEdge(prevAtom,currentAtom,prevBond)
                    
        # return the newly created atom
        prevBond = None
        return currentAtom,prevBond
        
    
    def parseOne(self,res,level,base_pos,pos,prevBond,prevAtom):
        if prevAtom != None :
            # each BigSMILES Bond ATOM should connect to at most one atom
            if self.G.nodes[prevAtom]['_type']=='BigSMILES_Bond':
                if not not self.G.nodes[prevAtom]['neighList']:
                    if 'dot' in res.keys():
                        pass
                    else:
                        errorMsg(self.rawStr,base_pos+pos,'Error','A BigSMILES bond should not connect to more than one atom.')
                        raise BigSMILESError

                # check if next bond is consistent with the right terminal
                if 'branch' in res.keys():
                    errorMsg(self.rawStr,base_pos+pos,'Error','The bonding descriptor must not be followed by a branch.')
                    raise BigSMILESError
                elif 'dot' in res.keys():
                    pass
                else:
                    if 'bond' not in res.keys():
                        trailing_bond = ''
                        if 'ringbond' in res.keys(): 
                            trailing_bond = res.ringbondtype
                        else:
                            if prevBond == None:
                                trailing_bond = '-'
                                prevBond = trailing_bond
                            else:
                                trailing_bond = prevBond
           
                            
                        #print(prevBond)
                        #print(self.G.nodes[prevAtom]['BigSMILES_Bond'].getS_bond())
                        # 
                        if self.G.nodes[prevAtom]['BigSMILES_Bond'].getS_bond() == 'u':
                            self.G.nodes[prevAtom]['BigSMILES_Bond'].setS_bond(trailing_bond)
                            if self.G.nodes[prevAtom]['BigSMILES_Bond'].checkConsistency(self.Bond_Dict,prevAtom) == -1:
                                errorMsg(self.rawStr,base_pos+pos,'Error','Inconsistent bond trailing the bonding descriptor.')
                                raise BigSMILESError
                            
                        if self.G.nodes[prevAtom]['BigSMILES_Bond'].getBondOrder() != self.G.nodes[prevAtom]['BigSMILES_Bond'].getBondOrder(trailing_bond):
                            errorMsg(self.rawStr,base_pos+pos,'Error','Inconsistent bond trailing the bonding descriptor.')
                            raise BigSMILESError

                    else:
                        pass
                        

                
            if self.G.nodes[prevAtom]['_type']=='BigSMILES_StoObj':
                # Parsing a BigSMILES stochastic object will result in an atom and a prevBond=None.
                if self.G.nodes[prevAtom]['BigSMILES_StoObj'].rightEnd == None: 
                    errorMsg(self.rawStr,base_pos+pos-1,'Error','Missing right terminal bonding descriptor for BigSMILES stochastic object connected to other atoms.')
                    raise BigSMILESError
                # Parsing a BigSMILES stochastic object will result in an atom and a trailing bond.
                if self.G.nodes[prevAtom]['BigSMILES_StoObj'].rightEnd.S_bond == 'u': 
                    errorMsg(self.rawStr,base_pos+pos-1,'Error','Undefined right terminal bonding descriptor for BigSMILES stochastic object.')
                    raise BigSMILESError
                # check if next bond is consistent with the right terminal
                if 'branch' in res.keys():
                    errorMsg(self.rawStr,base_pos+pos,'Error','The stochastic object must not be followed by a branch.')
                    raise BigSMILESError
                elif 'dot' in res.keys():
                    pass
                else:
                    if 'bond' not in res.keys(): 
                        trailing_bond = ''
                        if 'ringbond' in res.keys(): 
                            trailing_bond = res.ringbondtype
                        else:
                            #print(prevBond)
                            if prevBond == None:
                                trailing_bond = '-'
                                prevBond = trailing_bond
                            else:
                                trailing_bond = prevBond
           
                        #print(prevBond)
                        #print(self.G.nodes[prevAtom]['BigSMILES_Bond'].getS_bond())
                        # 
                        if self.G.nodes[prevAtom]['BigSMILES_StoObj'].rightEnd.getBondOrder() != self.G.nodes[prevAtom]['BigSMILES_StoObj'].rightEnd.getBondOrder(trailing_bond):
                            errorMsg(self.rawStr,base_pos+pos,'Error','Inconsistent bond trailing the stochastic object.')
                            raise BigSMILESError

                    else:
                        pass
                    
 
            

        
#        if 'BigSMILES_ladderBond' in res.keys() or 'BigSMILES_Bond' in res.keys():
        if 'BigSMILES_Bond' in res.keys():
            try:
                prevAtom,prevBond = self.addBigSmilesBondAtom(res,prevAtom,base_pos+pos,prevBond)
            except BigSMILES_BondInconsistencyError:
                raise BigSMILESError
#            prevBond = None
            return prevBond,prevAtom
        elif 'BigSMILES_StoObj' in res.keys():
            try:
                prevAtom,prevBond = self.addBigSmilesStoObjAtom(res,prevAtom,base_pos+pos,prevBond)
            except BigSMILES_BondInconsistencyError:
                errorMsg(self.rawStr,base_pos+pos,'Error','Inconsistent bond leading to stochastic object.')
                raise BigSMILESError    
            except BigSMILES_StoObjMissingTerminalError:
                errorMsg(self.rawStr,base_pos+pos,'Error','Missing left terminal bonding descriptor for BigSMILES stochastic object connected to other atoms.')
                raise BigSMILESError   
            except BigSMILES_StoObjError:
                errorMsg(self.rawStr,base_pos+pos,'Error','In parsing BigSMILES stochastic object ['+str(len(self.StoObj_List)+1)+']')
                raise BigSMILESError   
            #prevBond = None
            return prevBond,prevAtom
        else:
            return SMILES.parseOne(self,res,level,base_pos,pos,prevBond,prevAtom)
    
    def writeAtom(self,prevAtom,thisAtom,thisBond,rotCount,swapCount,hcount):
#        print(self.G.nodes[thisAtom]['_type'])
#        print(thisBond)
#        if self.G.nodes[thisAtom]['_type'] == 'BigSMILES_Bond' or self.G.nodes[thisAtom]['_type'] == 'BigSMILES_ladderBond':
        if self.G.nodes[thisAtom]['_type'] == 'BigSMILES_Bond':
            #smilesStr = '[' + 'BS' + ':' + self.G.nodes[thisAtom]['BigSMILES_Bond']._bondtype + self.G.nodes[thisAtom]['BigSMILES_Bond']._id + ']'            
            smilesStr = '[' + self.G.nodes[thisAtom]['BigSMILES_Bond'].getCompKey() + ']'
            if self.noWriteBondDesc==False:
                if thisBond:
                    return thisBond + smilesStr
                else:
                    return smilesStr
            else:
                return ''    
#            if thisBond == None:
#                return smilesStr + self.G.nodes[thisAtom]['BigSMILES_Bond'].getS_bond(isFirst=True)
#            else:
#                return self.G.nodes[thisAtom]['BigSMILES_Bond'].getS_bond(isFirst=False) + smilesStr
        elif self.G.nodes[thisAtom]['_type'] == 'BigSMILES_StoObj':
            if thisBond:
                smilesStr = thisBond + '{Sobj[' + str(self.G.nodes[thisAtom]['StoObjId']) + ']}'
            else:
                smilesStr = '{Sobj[' + str(self.G.nodes[thisAtom]['StoObjId']) + ']}'
            if len(self.G.nodes[thisAtom]['neighList']) > 1:
                nextAtom = self.G.nodes[thisAtom]['neighList'][1]
                if self.G.nodes[nextAtom]['atom'] == 'H' and self.G.nodes[nextAtom]['_type'] == 'bracket_atom':
                    smilesStr = smilesStr + self.G.nodes[thisAtom]['BigSMILES_StoObj'].getBond(end='right') + '[H]'
            return smilesStr
        else:
            if self.noWriteBondDesc==True and prevAtom!=None:
                if self.G.nodes[prevAtom]['_type'] == 'BigSMILES_Bond':
                    thisBond = ''
            return SMILES.writeAtom(self,prevAtom,thisAtom,thisBond,rotCount,swapCount,hcount)
    
    
    def writeStandard(self,noBondDesc=False, forward=True):
        # write the string with backbone as the main chain
        # ends 
        if noBondDesc==True:
            self.noWriteBondDesc=True;
        ## get all the binding sites
        # get lists of binding sites corresponding to different BigSMILES bonding descriptors
        bonding_sites = flatten_list([x[1:] for x in self.Bond_Dict.values()])
        # exit if there are less than two bonding sites
        #if len(bonding_sites) < 2:
        #    print('Error in writing standardized repeat unit: expected at least 2 bonding sites but '+str(len(bonding_sites))+' found')
        #    return None
        # choose the first two as the ends
        
        if len(bonding_sites) == 0:
            smilesStr = self.write()
            self.noWriteBondDesc=False
            return smilesStr
        elif len(bonding_sites) == 1:
            source = bonding_sites[0]
            smilesStr = self.write(source)
            self.noWriteBondDesc=False
            return smilesStr
        else:
            if (forward):
                source = bonding_sites[0]
                target = bonding_sites[1]
            else:
                source = bonding_sites[1]
                target = bonding_sites[0]
##            source = bonding_sites[0]
##            target = bonding_sites[1]
        
        # get the backbone (defined as the shortest path between the two ends)
        path = nx.shortest_path(self.G,source=source,target=target)
        #print(len(path))
        #print(path)
        #G_copy = self.G.copy() # don't bother to copy
        
        for i in range(1,len(path)-1):
            atom = path[i]
            next_atom = path[i+1]
            prev_atom = path[i-1]
            L = self.G.nodes[atom]['neighList']
            
            # first rotate the list L so that the previous atom is at the front
            rotCount = 0
            while L[0] != prev_atom:
                L = deque(L)
                L.rotate(-1)
                rotCount += 1
            L = list(L)
            
            # then check if the next_atom need to be swapped to the end
            swapCount = 0
            if L.index(next_atom) == len(L)-1: 
                #print('no swapp on atom '+str(atom))
                # no need to change things around as the next_atom is at the end of the list
                # i.e. on main chain and not on branch
                pass
            else:
                #print('swapped bonds on atom '+str(atom))
                # need to swap the next_atom with the last atom
                L[L.index(next_atom)] = L[-1]
                L[-1] = next_atom
                self.G.nodes[atom]['neighList'] = L
                swapCount += 1
                #print(L)
            
            # change chirality accordingly
            if self.G.nodes[atom]['chiral'] == '':
                # no chirality specified, whew, nothing to change
                pass
            else:
                nchiral = len(self.G.nodes[atom]['chiral']) -1
                nchiral = ((nchiral + rotCount + swapCount) % 2) + 1
                self.G.nodes[atom]['chiral'] = '@'*nchiral
        
        # done changing G_copy, start from source, and writeLinear will give chain according to path
            
        # now deal with the loops
        # get bfs tree from source. this would ensure all loops are cut on chains other than the main desired one
        self.T = nx.bfs_tree(self.G,source,reverse=True).to_undirected()
        
        self.ringDict = dict()
        self.usedRingID = [False]*100
        for edge in self.G.edges():
            if not tuple(edge) in self.T.edges():
                self.ringDict[edge] = -1
        
        #print(self.ringDict)
        #print(self.T.edges())
        #smilesStr = self.writeLinear((None,source))
        
        
        smilesStr = self.writeComponents(source)
        self.noWriteBondDesc=False
        
        return smilesStr
    
    
    def getFunctionality(self):
        bonding_sites = [x for x in self.Bond_Dict.values()]
        bonding_sites_flat = flatten_list([x[1:] for x in self.Bond_Dict.values()])
        totF = len(bonding_sites_flat)
        noTypeBond = len(bonding_sites)
        return totF,noTypeBond,[x[0].getCompleteSymbol() for x in bonding_sites],[len(x)-1 for x in bonding_sites]
    
    
    def numNeighAll(self):
        valency = {'C': 4, 'c':3, 'H': 1, 'O': 2, 'o':1, 'N': 3, 'n': 2, 'S': 2, 's':1, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1, 'Si': 4}
        typeToNum = {'': 1, '/': 1,'\\': 1, '-': 1, '=': 2,'#': 3}  
        numNeighAll = {}
        for node in self.G.nodes:
            if self.G.nodes[node]['atom'] == "BigSMILES_Bond":
                numNeighAll[node] = "BigSMILES_Bond"
            else:
                atomType = self.G.nodes[node]['atom']
                numNeigh = valency.get(atomType)
                for neigh in self.G.nodes[node]['neighList']:
                    bondType = self.G.edges[node, neigh]['type']
                    numNeigh -= typeToNum.get(bondType) - 1
                atomBond = atomType.upper() + str(numNeigh)
                numNeighAll[node] = atomBond
        return numNeighAll
    
    
    def path(self, start, length, visited= None):
        paths = []
        if visited is None:
            visited = []
        visited.append(start) 
        if length == 1:
            paths.append(visited)
        else:
            for node in self.G.adj[start]:
                visitedCopy = visited[:]
                if node not in visitedCopy:
                    paths.extend(self.path(node, length-1, visitedCopy))
        return paths
    
    
    def multiplet(self, size):
        numNeigh = self.numNeighAll()
        multiplets = {}
        for node in self.G.nodes:
            paths = self.path(node, size)
            for path in paths:
                multiplet = []
                for n in path:
                    multiplet.append(numNeigh.get(n))
                if 'BigSMILES_Bond' not in multiplet:
                    multiplet = tuple(multiplet)
                    if str(multiplet) in multiplets:
                        multiplets[str(multiplet)]+=1/2
                    elif str(multiplet[::-1]) in multiplets:
                        multiplets[str(multiplet[::-1])]+=1/2
                    else:
                        multiplets[str(multiplet)]=1/2 
        return multiplets

#######################################################################################################################    

import rdkit, collections, itertools
from rdkit import Chem
import pandas as pd
import os
def chemistry_table():
    chemistries = pd.read_csv(os.path.join(os.path.dirname( __file__ ), os.pardir, "names_dict.csv"))
    chemistries = pd.DataFrame(data = chemistries)
    chemistries = chemistries.fillna("")
    
    exact_hits = []
    isCycle = []
    index = 0
    for s in chemistries['Object']:
        Polymer = BigSMILES(s)
        G = Polymer.G
        repeat_list = []
        if chemistries['Cycle'].iloc[index] !=  "":
            isCycle.append(True)
            repeat_list = chemistries['Cycle'].iloc[index].split(",")
        else:
            isCycle.append(False)
            for obj in Polymer:
                for rep in obj:
                    repeat_list.append(rep.writeStandard(noBondDesc = True, forward = True))
        index += 1
        exact_hits.append(repeat_list)

    chemistries.insert(1, 'Exact_Hit', exact_hits)
    chemistries.insert(1, 'isCycle', isCycle)

    def join(block, n = 3):
        combos = [p for p in itertools.product(block, repeat = n)]
        unique_combos = set()
        for c in combos:
            unique_combos.add(''.join(str(elem) for elem in c))
        return list(unique_combos)

    substructure_hits = []
    for s in chemistries['Object']:
        Polymer = BigSMILES(s)
        G = Polymer.G
        nbond = []
        for obj in Polymer:
            for rep in obj:
                nbond.append(rep.writeStandard(noBondDesc = True, forward = True))
                if '$' in rep.writeStandard(noBondDesc = False, forward = True):
                    nbond.append(rep.writeStandard(noBondDesc = True, forward = False))
        substructure_hits.append(join(nbond, 3))
    chemistries.insert(1, 'Substructure_Hit', substructure_hits)
    # print(chemistries)
    return chemistries

#######################################################################################################################

if __name__ == "__main__":
    print("BigSMILES is the future of polymer informatics")
    from BigSMILES_BigSmilesObj import BigSMILES
    plga = '{[][<]C(=O)C(C)O[>],[<]C(=O)CO[>][]}'
    Polymer = BigSMILES(plga)
    polymerGraph = Polymer.G
    repeat_list = []
    for obj in Polymer:
        for rep in obj:
            repeat_list.append(rep.writeStandard(noBondDesc = True))
    print(repeat_list)
