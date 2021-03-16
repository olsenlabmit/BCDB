# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 10:28:21 2019

@author: tslin
"""
# TODO: tran-cis and extended allene-like systems
# tran-cis would require directional bond

from SmilesPattern import SmilesPattern
from utility import errorMsg, disjoint_union, getIdStr
import networkx as nx
import matplotlib.pyplot as plt
from collections import deque
from error import SmilesError, NonmatchingSymbolError, NonmatchingTransCisError, SelfLoopError, \
                  MultiBondError, HcountChiralError, AtomicWeightError


class SMILES:
    _count = 0
    def __init__(self,inStr="",inGraph=None,pattern=SmilesPattern._smilesElement):
        self.parsed = False
        self.pattern = pattern
        if inGraph == None:
            pass
        else:
            self.rawStr = ""
            self.G = inGraph
            self.T = nx.Graph()
            return
        
        try:
            self.rawStr = inStr.split()[0]
        except:
            self.rawStr = ""
        self.atomCount = 0
        self.prevAtomStack = list()
        self.position = 0
        self.G = nx.Graph()
        self.T = nx.Graph()
        self.ringDict = dict()
        self.atomPos = list()
        self.usedRingID = list()
        SMILES._count = SMILES._count + 1
        self.MolID = SMILES._count
        self.components = list()
        
    def parseBranch(self,res,prevAtom,level,pos):
        self.parse(root=prevAtom, level=level+1, inStr=res.rawStr[1:-1], base_pos=pos+1)
        #print(res.rawStr)
        
    def addBond(self,res,prevAtom):
        return res.bond
        #print(res.rawStr)
        
    def addRingBond(self,res,prevAtom):
        parsedLoopId = int(res.ringid)
        parsedBond = res.ringbondtype
        
        # case 1: ring id already exist, make ring connection
        if parsedLoopId in self.ringDict:    
            loopAtom = self.ringDict[parsedLoopId][0]
            loopBond = self.ringDict[parsedLoopId][1]
            
            if loopBond == '':
                pass
            elif parsedBond == '':
                # flip bond type for cis-trans consistency
                if loopBond == '\\':
                    parsedBond = '/'
                elif loopBond == '/':
                    parsedBond = '\\'
                else:
                    parsedBond = loopBond
            
            else:
                # check if the two occurrences of the ring bond have matching symbols
                # for usual non-directional -=# bonds
                if parsedBond != '\\' and  parsedBond != '/':
                    if parsedBond == loopBond:
                        pass
                    else:
                        raise NonmatchingSymbolError
                    #return
                # for trans-cis bonds
                else:
                    if parsedBond == '\\':
                        if loopBond == '/':
                            pass
                        else:
                            raise NonmatchingTransCisError
                    else:
                        if loopBond == '\\':
                            pass
                        else:
                            raise NonmatchingTransCisError
            
            # check if illegal ring formation occurs
            if loopAtom == prevAtom: # self-loop
                raise SelfLoopError
                return
            if self.G.has_edge(loopAtom,prevAtom): # multiple edges between a pair of atoms
                raise MultiBondError
                return
            
            #self.G.add_edge(loopAtom,prevAtom)
            #self.G.edges[loopAtom,prevAtom]['type'] = parsedBond
            self.G.add_edge(prevAtom,loopAtom)
            self.G.edges[prevAtom,loopAtom]['type'] = parsedBond
            # find in neighList and replace to maintain order
            idx = self.G.nodes[loopAtom]['neighList'].index(-1*parsedLoopId)
            self.G.nodes[loopAtom]['neighList'][idx] = prevAtom
            self.G.nodes[prevAtom]['neighList'].append(loopAtom)
            self.G.edges[prevAtom,loopAtom]['direction'] = (prevAtom,loopAtom)
            # remove the ring id from dictionary
            self.ringDict.pop(parsedLoopId)
            
        # case 2: new ring id
        else:
            self.ringDict[parsedLoopId] = [prevAtom,parsedBond]
            #G.add_edge(prevAtom,phantomLoopNodeIdx)
            self.G.nodes[prevAtom]['neighList'].append(-1*parsedLoopId)
        
    
    def createNode(self,rawStr,_type,symbol,pos,chiral='',isotope='',charge='',_class=''):
        self.atomCount = self.atomCount + 1
        nodeId = self.atomCount
        self.G.add_node((nodeId))
        self.G.nodes[nodeId]['rawStr'] = rawStr
        self.G.nodes[nodeId]['_type'] = _type
        self.G.nodes[nodeId]['atom'] = symbol
        self.G.nodes[nodeId]['neighList'] = list()
        self.G.nodes[nodeId]['pos'] = pos
        self.G.nodes[nodeId]['chiral'] = chiral
        
        self.G.nodes[nodeId]['isotope'] = isotope
        self.G.nodes[nodeId]['charge'] = charge
        self.G.nodes[nodeId]['_class'] = _class
        return nodeId

    def createEdge(self,node1,node2,bond):
        self.G.add_edge(node1,node2)
        self.G.nodes[node1]['neighList'].append(node2)
        self.G.nodes[node2]['neighList'].append(node1)
        if bond == None:
            self.G.edges[node1,node2]['type'] = ''
        else:
            self.G.edges[node1,node2]['type'] = bond
        # add directionality of bond, specifically needed for reconstructing trans-cis isomers
        self.G.edges[node1,node2]['direction'] = (node1,node2)

    def addAtom(self,res,prevAtom,prevBond,pos):
        # add the new atom into the graph
        if 'orgainc' in res.keys():
            _type = 'orgainc'
        elif 'bracket_atom' in res.keys():
            _type = 'bracket_atom'
        else:
            _type = ''
        
        currentAtom = self.createNode(res.rawStr,_type,res.symbol,pos,''.join(res.chiral),''.join(res.isotope),''.join(res.charge),''.join(res._class))
        numH = 0
        
        if res.symbol == 'H':
            if 'isotope' in res.keys():
                W = int(res.isotope)
                if W > 3 or W < 1:
                    raise AtomicWeightError     
        else:
            if 'isotope' in res.keys():
                W = int(res.isotope)
                if W < 1:
                    raise AtomicWeightError
                    
        # make connection with previous atom in stack
        if prevAtom == None:
            pass
        else:
            self.createEdge(prevAtom,currentAtom,prevBond)
            
        # deal with h-count
        if 'hcount' in res.keys():
            if res.nH == '':
                numH = 1
            else:
                numH = int(res.nH)
        for i in range(numH):
            tmpHAtom = self.createNode('','organic','H',pos)
            self.createEdge(currentAtom,tmpHAtom,'-')
        
        # deal with sp3 chirality
        if 'chiral' in res.keys():
            if numH > 1:
                raise HcountChiralError
            self.G.nodes[currentAtom]['chiral'] = res.chiral
        else:
            self.G.nodes[currentAtom]['chiral'] = ''
        
        
        
        
        # return the newly created atom
        return currentAtom           

    def finalCheck(self,level,base_pos,pos,prevBond,prevAtom):
        # entire chain/branch parsed, doing final checks and additional treatment
        if prevBond != None: # the chain ends with an open trailing bond, append explicit hydrogen to that
            tmpHAtom = self.createNode('[H]','bracket_atom','H',pos)
            self.createEdge(prevAtom,tmpHAtom,prevBond)
            # deprecated, not considered as error
            #errorMsg(self.rawStr,base_pos+pos,'Error','missing ATOM, expected ATOM trailing a bond')
            #raise SmilesError
            
        # final check after the main chain is completely parsed
        if level == 0:
            # check for any open ring
            if not not self.ringDict:
                for i in self.ringDict.keys():
                    errorMsg(self.rawStr,base_pos+pos,'Error','missing ring-bond, ring ID = '+str(i)+' not closed')
                raise SmilesError
            
            # check if app / \ bonds are adjacent to double bonds
#            for edge in self.G.edges():
#                if self.G.edges[edge]['type'] != '\\' and self.G.edges[edge]['type'] != '/':
#                    continue
#                else:
#                    doubleBond = False
#                    for node1 in edge:
#                        neighbors = self.G.nodes[node1]['neighList']
#                        for node2 in neighbors:
#                            if node2 in edge:
#                                continue
#                            else:
#                                if self.G.edges[node1,node2]['type'] == '=':
#                                    doubleBond = True
#                                    break
#                        if doubleBond:
#                            break
#                    if not doubleBond:
#                        errorMsg(self.rawStr,self.G.nodes[node1]['pos'],'Error',r'\ or / should always be adjacent to double bonds.')
#                        raise SmilesError
            
            # fill inferrable trans-cis bonds and check consistency
            for edge in self.G.edges():
                if self.G.edges[edge]['type'] != '=':
                    continue
                else:
                    for node1 in edge:
                        u = 0
                        d = 0
                        neighbors = self.G.nodes[node1]['neighList']
                        unspecified = list()
                        skipflag = False
                        for node2 in neighbors:
                            if node2 in edge:
                                continue
                            else:
                                if self.G.edges[node1,node2]['type'] == '=':
                                    skipflag = True
                                    break
                                elif self.G.edges[node1,node2]['type'] == '\\':
                                    if self.G.edges[node1,node2]['direction'] == (node1,node2):
                                        d += 1
                                    else:
                                        u += 1
                                elif self.G.edges[node1,node2]['type'] == '/':
                                    if self.G.edges[node1,node2]['direction'] == (node1,node2):
                                        u += 1
                                    else:
                                        d += 1
                                else:
                                    pass
                                if self.G.edges[node1,node2]['type'] != '\\' and self.G.edges[node1,node2]['type'] != '/':
                                    unspecified.append(node2)
                        # skip to next node if both connection of the atom are double bonds
                        if skipflag:
                            continue
                        # otherwise check consistency
                        if u > 1 or d > 1:
                            errorMsg(self.rawStr,self.G.nodes[node1]['pos'],'Error','errorneous isomeric specification; cannot have two bonds on an allenal carbon pointing in the same relative "up" or "down" direction.')
                            raise SmilesError
                        # complete inferrable up/down specification
                        # case 1: up and down completely agnostic
                        if u+d == 0:
                            continue
                        # case 2: up and down fully specified
                        elif u+d == 2:
                            continue
                        # case 3: only one of the two bonds are specified
                        else:
                            if len(unspecified) > 1:
                                errorMsg(self.rawStr,self.G.nodes[node1]['pos'],'Error','errorneous isomeric specification; an allenal carbon with valency larger than 4 is found.')
                                raise SmilesError
                            elif len(unspecified) < 1:
                                pass
                            else:
                                node2 = unspecified[0]
                                if u == 0:
                                    self.G.edges[node1,node2]['type'] = '/'
                                else:
                                    self.G.edges[node1,node2]['type'] = '\\'
                                self.G.edges[node1,node2]['direction'] = (node1,node2)
            
            # check if all chiral centers have 4 bonds
            flag = False
            for node in self.G.nodes():
                if self.G.nodes[node]['_type'] == 'organic' or self.G.nodes[node]['_type'] == 'bracket_atom':
                    if self.G.nodes[node]['chiral'] != '':
                        # TODO: allene chiral centers
                        if len(self.G.nodes[node]['neighList']) != 4:
                            errorMsg(self.rawStr,self.G.nodes[node]['pos'],'Error','errorneous chiral center specification, chiral centers should always have 4 distinct neighbors')
                            flag = True
                        else:
                            numH = [0,0,0]
                            for atom in self.G.nodes[node]['neighList']:
                                if self.G.nodes[atom]['atom'] == 'H':
                                    if self.G.nodes[atom]['isotope'] == '' or self.G.nodes[atom]['isotope'] == '1':
                                        numH[0] += 1
                                    elif self.G.nodes[atom]['isotope'] == '2':
                                        numH[1] += 1
                                    elif self.G.nodes[atom]['isotope'] == '3':
                                        numH[2] += 1
                            if numH[0] > 1 or numH[1] > 1 or numH[2] > 1:
                                errorMsg(self.rawStr,self.G.nodes[node]['pos'],'Error','the ATOM cannot be chiral with chirality '+self.G.nodes[atom]['chiral']+'with more than one H/D/T attached')
                                flag = True
        
                    else:
                        continue
            if flag:
                raise SmilesError

    def parseError(self,tmpStr,base_pos,pos):
        if tmpStr[0] == '(':
            #raise Exception(errorMsg(self.rawStr,base_pos+pos,'unbalanced parenthesis'))
            errorMsg(self.rawStr,base_pos+pos,'Error','in parsing branch, unbalanced parenthesis or illegal expression within branch')
            raise SmilesError
        elif tmpStr[0] == ')':
            errorMsg(self.rawStr,base_pos+pos,'Error','unbalanced parenthesis')
            raise SmilesError
        else:
            #raise Exception(errorMsg(self.rawStr,base_pos+pos,'cannot recognize element'))
            errorMsg(self.rawStr,base_pos+pos,'Error','cannot recognize element')
            raise SmilesError

    def parseOne(self,res,level,base_pos,pos,prevBond,prevAtom):
        # the case next element is a dot
        if 'dot' in res.keys():
            if prevBond != None: # the chain ends with an open trailing bond, append explicit hydrogen to that
                tmpHAtom = self.createNode('[H]','bracket_atom','H',pos)
                self.createEdge(prevAtom,tmpHAtom,prevBond)
            # check if dot is within branches # support this unstandard notation
            #if level > 0:
            #    errorMsg(self.rawStr,base_pos+pos,'dot can only be on base level and cannot appear within parenthesizes')
            #    raise SmilesError
            # clear prevAtom
            prevAtom = None
            
        # the case next element is a branch
        elif 'branch' in res.keys():
            if len(res.rawStr) == 2:
                errorMsg(self.rawStr,base_pos+pos,'Warning','missing expressions within parenthesizes')
                #raise Warning
            if prevAtom == None:
                errorMsg(self.rawStr,base_pos+pos,'Error','expected ATOM or BOND but got BRANCH, a branch should always follow an atom')
                raise SmilesError
            if pos == 0:
                #if level == 0:
                #    errorMsg(self.rawStr,base_pos+pos,'Error','expected ATOM or BOND but got BRANCH, a branch should always follow an atom')
                #    raise SmilesError
                #else:
                errorMsg(self.rawStr,base_pos+pos,'Warning','ignoring extra parentheses, expected ATOM or BOND but got BRANCH, a branch should always follow an atom')
                    #raise Warning
            if prevBond != None:
                errorMsg(self.rawStr,base_pos+pos,'Warning','potentially misleading expression, a branch should follow immediately an ATOM and the BOND on the main chain should come directly before the corresponding ATOM.')
                #expected ATOM but got BRANCH, a branch should always follow an atom'))
                #raise SmilesError
            self.parseBranch(res,prevAtom,level,base_pos+pos)
            
        # next element is a bond
        elif 'bond' in res.keys():
            if prevAtom == None:
                errorMsg(self.rawStr,base_pos+pos,'Error','expected ATOM but got BOND, a BOND should always follow an ATOM')
                raise SmilesError
                #raise Exception(errorMsg(self.rawStr,base_pos+pos,'expected ATOM but got BOND, a bond should always follow an atom'))
            if prevBond != None:
                if prevBond == self.addBond(res,prevAtom):
                    #errorMsg(self.rawStr,base_pos+pos,'Warning','ignoring repeated BOND symbol')
                    # raise Warning
                    pass
                else:
                    errorMsg(self.rawStr,base_pos+pos,'Error','expected ATOM but got BOND, a BOND should always be followed by an ATOM')
                    raise SmilesError
            
            prevBond = self.addBond(res,prevAtom)
            
        # next element is a ring-bond
        elif 'ringbond' in res.keys():
            if prevAtom == None:
                errorMsg(self.rawStr,base_pos+pos,'Error','expected ATOM but got RING-BOND, a RING-BOND should always follow an ATOM')
                raise SmilesError                    
                #raise Exception(errorMsg(self.rawStr,base_pos+pos,'expected ATOM but got RINGBOND, a ring-bond should always follow an atom'))
                
            if pos == 0:
                errorMsg(self.rawStr,base_pos+pos,'Warning','potentially ambiguous expression, a RING-BOND should directly follow an ATOM')
                #raise SmilesError
                #raise Exception(errorMsg(self.rawStr,base_pos+pos,'expected ATOM or BOND but got RINGBOND, a ring-bond should always follow an atom'))
            
            if prevBond != None:
                errorMsg(self.rawStr,base_pos+pos,'Warning','potentially misleading expression, a RING-BOND should follow immediately an ATOM and the BOND on the main chain should come directly before the corresponding ATOM.')
                #expected ATOM but got RINGBOND, a ring-bond should always follow an atom'))
                #raise SmilesError
                
            
            try:
                self.addRingBond(res,prevAtom)
            except NonmatchingSymbolError:
                errorMsg(self.rawStr,base_pos+pos,'Error','bond symbol does not match for RING-BONDS with the same ID')
                raise SmilesError
                #raise Exception(errorMsg(self.rawStr,base_pos+pos,'bond symbol does not match for ring bonds with the same ID'))
            except NonmatchingTransCisError:
                errorMsg(self.rawStr,base_pos+pos,'Error','cis-trans in RING-BOND is not consistent')
                raise SmilesError
            except SelfLoopError:
                errorMsg(self.rawStr,base_pos+pos,'Error','self-loops on a single ATOM is not allowed')
                raise SmilesError
            except MultiBondError:
                errorMsg(self.rawStr,base_pos+pos,'Error','more than one BOND between a pair of atoms is found')
                raise SmilesError
            
        # otherwise, the next element is an atom     
        else:
            try: 
                prevAtom = self.addAtom(res,prevAtom,prevBond,base_pos+pos)
            except HcountChiralError:
                errorMsg(self.rawStr,base_pos+pos,'Error','the ATOM cannot be chiral with chirality '+res.chiral+'with more than one Hydrogen attached')
                raise SmilesError
            except AtomicWeightError:
                errorMsg(self.rawStr,base_pos+pos,'Error','erroneous isotope specified.')
                raise SmilesError
            prevBond = None
        
        return prevBond,prevAtom


    def parse(self,root=None,level=0,inStr=None,base_pos=0):
        if self.parsed == True:
            return self.G
        if root == None:
            inStr = self.rawStr
        
        pos = 0
        tmpStr = inStr
        prevAtom = root
        prevBond = None
        
        # iteratively read from left to right until the entire string is read
        while pos < len(inStr):
            tmpStr = inStr[pos:]
            try:
                res = self.pattern.parseString(tmpStr)
                #print(tmpStr)
                #print(res.asDict())
            except:
                self.parseError(tmpStr,base_pos,pos)
            prevBond,prevAtom = self.parseOne(res,level,base_pos,pos,prevBond,prevAtom)
            
            # iterate over the input string
            pos = pos + len(res.rawStr)
        
        # entire chain/branch parsed, doing final checks and additional treatment
        self.finalCheck(level,base_pos,pos,prevBond,prevAtom)
          # reduce allenic chiral centers
            # TODO: allene chiral centers
        
        if level == 0:
            for g in nx.connected_components(self.G):
                self.components.append(g)
            
            self.parsed = True    
        return self.G
        
    
    def swapBranch(self,swapStart,thisAtom,swapList,condition):
        swapIdx = swapStart
        swapCount = 0
        
        while swapIdx < len(swapList):
            atom = swapList[swapIdx]
            if condition(atom):
                swapIdx += 1
                continue
            else:
                flag = False
                for i in range(swapIdx+1,len(swapList)):
                    atom2 = swapList[i]
                    if condition(atom2):
                        swapList[swapIdx] = atom2
                        swapList[i] = atom
                        swapCount += 1
                        swapIdx += 1
                        flag = True
                        break
                    else:
                        continue
                if not flag:
                    break
        
        return swapCount,swapIdx-swapStart
    
    
    def getBondType(self,prevAtom,thisAtom):
        thisBond = self.G.edges[prevAtom,thisAtom]['type']
        if thisBond == '\\' or thisBond == '/':
            direction = self.G.edges[prevAtom,thisAtom]['direction']
            if direction == (prevAtom,thisAtom):
                pass
            else:
                if thisBond == '\\':
                    thisBond = '/'
                else:
                    thisBond = '\\'
        return thisBond
    
    def writeAtom(self,prevAtom,thisAtom,thisBond,rotCount,swapCount,hcount):
        smilesStr = ''
        # obtain correct chiral symbol
        if self.G.nodes[thisAtom]['chiral'] =='':
            chiralSymbol = ''
        else:
            nchiral = len(self.G.nodes[thisAtom]['chiral']) -1
            nchiral = ((nchiral + rotCount + swapCount) % 2) + 1
            chiralSymbol = '@'*nchiral
        # obtain correct hcount
        if hcount == 0:
            hcountSymbol = ''
        elif hcount == 1:
            hcountSymbol = 'H'
        else:
            hcountSymbol = 'H'+str(hcount)

        if thisBond == None:
            thisBond = ''
        # write bracket atom
        #print(self.G.nodes[thisAtom])
        if self.G.nodes[thisAtom]['isotope'] != '' \
           or self.G.nodes[thisAtom]['chiral'] != '' \
           or self.G.nodes[thisAtom]['charge'] != '' \
           or self.G.nodes[thisAtom]['_class'] != '':
            smilesStr = smilesStr + thisBond + '[' \
                        + self.G.nodes[thisAtom]['isotope'] \
                        + self.G.nodes[thisAtom]['atom'] \
                        + chiralSymbol \
                        + hcountSymbol \
                        + self.G.nodes[thisAtom]['charge'] \
                        + self.G.nodes[thisAtom]['_class'] \
                        + ']'
        
        # or normal atom
        else:
            if self.G.nodes[thisAtom]['atom'] != 'H':
                ATOM = self.G.nodes[thisAtom]['atom']
                if ATOM!='Cl' and ATOM!='Br' and ATOM!='B' and ATOM!='C' and ATOM!='N'and \
                   ATOM!='O' and ATOM!='S' and ATOM!='P' and ATOM!='P' and ATOM!='F' and ATOM!='I' and \
                   ATOM!='c' and ATOM!='b' and ATOM!='n' and ATOM!='o' and ATOM!='s' and ATOM!='p':
                    smilesStr = smilesStr + thisBond + '[' + self.G.nodes[thisAtom]['atom'] + ']'
                else:
                    smilesStr = smilesStr + thisBond + self.G.nodes[thisAtom]['atom']
            else:
                if self.G.nodes[thisAtom]['_type'] == 'bracket_atom':
                    smilesStr = smilesStr + thisBond + '[H]' # TODO
                
        
        return smilesStr
    
    
    def writeLinear(self,base):
        smilesStr = ''
        
        
        prevAtom = base[0]
        thisAtom = base[1]
        
        
        while thisAtom != None:
            nextAtoms = self.G.nodes[thisAtom]['neighList']
            #nextAtomsTree = list(self.T[thisAtom])
            
            rotCount = 0
            if prevAtom != None:
                while nextAtoms[0] != prevAtom:
                    nextAtoms = deque(nextAtoms)
                    nextAtoms.rotate(-1)
                    rotCount += 1
                nextAtoms = list(nextAtoms)
                
            swapCount = 0
            if prevAtom != None:
                start = 1
            else:
                start = 0
                
            hswap,hcount = self.swapBranch(start,thisAtom,nextAtoms,
                                           lambda atom:self.G.nodes[atom]['atom'] == 'H' and \
                                                       self.G.nodes[atom]['isotope'] == '' and \
                                                       self.G.nodes[atom]['charge'] == '' and \
                                                       self.G.nodes[atom]['_class'] == '' )
            ringswap,ringcount = self.swapBranch(start,thisAtom,nextAtoms,
                                                 lambda atom:tuple(sorted((thisAtom,atom))) in self.ringDict.keys())
            swapCount = hswap + ringswap
            hcount + ringcount
            
            # write bond
            if prevAtom != None:
                thisBond = self.getBondType(prevAtom,thisAtom)
            else:
                thisBond = None
            # defer writing bond to check if next atom is H
            #smilesStr = smilesStr + thisBond
            
            # write atom
            smilesStr = smilesStr + self.writeAtom(prevAtom,thisAtom,thisBond,rotCount,swapCount,hcount)
        
            # write ring branches
            for nextAtom in nextAtoms[start+hcount:start+hcount+ringcount]:
                ringbond = self.getBondType(thisAtom,nextAtom)
                edge = tuple(sorted((thisAtom,nextAtom)))
                if self.ringDict[edge] == -1:
                    ringID = self.usedRingID.index(False)
                    self.ringDict[edge] = ringID
                    self.usedRingID[ringID] = True
                else:
                    ringID = self.ringDict[edge]
                    self.ringDict[edge] = -2
                    self.usedRingID[ringID] = False
                
                ringIDstr = getIdStr(ringID+1)
                
                smilesStr = smilesStr + ringbond + ringIDstr
                
            # write branches
            for nextAtom in nextAtoms[start+hcount+ringcount:-1]:
                smilesStr = smilesStr + '(' + self.writeLinear((thisAtom,nextAtom)) +')'
            
            prevAtom = thisAtom
            if start+hcount+ringcount == len(nextAtoms):
                thisAtom = None
            else:
                thisAtom = nextAtoms[-1]
                continue
        
        return smilesStr
        

    
    def writeComponents(self,source):
        tmpComponents = self.components.copy()
        smilesStrs = list()
        smilesStrs.append(self.writeLinear((None,source)))
        
        for item in tmpComponents:
            if source in item:
                tmpComponents.remove(item)
                break
        
        for item in tmpComponents:
            if len(item) > 0:
                base = next(iter(item))
                smilesStrs.append(self.writeLinear((None,base)))
        
        return '.'.join(smilesStrs)
         
    
    
    def write(self,base=1):

        
        self.T = nx.minimum_spanning_tree(self.G)
        
        self.ringDict = dict()
        self.usedRingID = [False]*100
        for edge in self.G.edges():
            if not tuple(edge) in self.T.edges():
                self.ringDict[edge] = -1
        
        #print(self.ringDict)
        #print(self.T.edges())
               
        
        smilesStr = self.writeComponents(base)
        #return smilesStr,self.G,self.T
        return smilesStr
        
            
            
    
    
    def draw(self):
        colorDict = { # CPK color code
                'H' : 'w', 
                'C' : 'k', 
                'c' : 'k',
                'N' : 'b', 
                'n' : 'b', 
                'O' : 'r', 
                'o' : 'r', 
                'F' : 'g', 
                'Cl': 'g', 
                'Br': 'brown',  
                'I' : 'purple', 
                'S' : 'y'  }        
        color_map = []
        for node in self.G:
            atomType = self.G.nodes[node]['atom']
            if atomType in colorDict:
                color_map.append(colorDict[atomType])
            else:
                color_map.append('pink') # should be pink
        
        edgeColorDict = { 
                ''  : 'k', 
                '-' : 'k', 
                '=' : 'r', 
                '#' : 'purple',
                '/' : 'b', 
                '\\' : 'g', 
                ':' : 'k' }        
        edge_colors = []
        for edge in list(self.G.edges()):
            edgeType = self.G.edges[edge]['type']
            if atomType in colorDict:
                edge_colors.append(edgeColorDict[edgeType])
            else:
                edge_colors.append('w') # should be pink
        atoms = nx.get_node_attributes(self.G, 'atom')
        chirals = nx.get_node_attributes(self.G, 'chiral')
        labels = dict()
        for key in atoms.keys():
            labels[key] = str(key)+'\n'+atoms[key]+chirals[key]
        labelsH = dict()
        for key in atoms.keys():
            if atoms[key] == 'H':
                labelsH[key] = str(key)+'\n'+atoms[key]+chirals[key]
            else:
                labelsH[key] = ''
        # https://networkx.github.io/documentation/stable/reference/drawing.html
        #nx.drawing.nx_pylab.draw_kamada_kawai(self.G, 
        #nx.drawing.nx_pylab.draw_spring(self.G,
        plt.figure()
        position=nx.spring_layout(self.G)
        nx.drawing.nx_pylab.draw_networkx(self.G,
                                          pos = position,
                                          with_labels = True,
                                          labels = labels,
                                          font_color = 'w',
                                          font_weight = 'bold',
                                          node_size = 1000,
                                          node_color = color_map,
                                          edge_color=edge_colors)
        nx.draw_networkx_labels(self.G,
                                pos = position,
                                font_weight = 'bold',
                                font_color = 'k',
                                labels = labelsH )
        #nx.draw_networkx_nodes(self.G,position,node_color = color_map,node_shape='o',node_size = 1000)

if __name__ == "__main__":
    testStr = 'C1(O[2H:1])=CC=CC=C1I'
    testStr = 'C[C@H](Cl)I'
    #testStr = 'C[C@](Cl)(I)[H]'
    testStr = 'c1ccccc1'
    #testStr = 'C1CCCCC1'
    #testStr = 'C(=2)[CH2](C)(C)CC2'
    #testStr = 'OCC(CCC)C(C(C)C)CCC'
    #testStr = 'N1CC2CCCC2CC1'
    #testStr = 'C12(CCCCC1)CCCCC2'
    #testStr = '[H]C([2H])([H])[H]'
    #testStr = 'Oc1ccccc1.NCCO'
    #testStr = 'c1c2c3c4cc1.Br2.Cl3.Cl4'
    #testStr = 'C0CCCCC0'
    #testStr = 'C.CC(C)C.C' 
    #testStr = 'C1CCCCC%01'
    #testStr = 'C%012(CCCCC1)CCCCC2'
    #testStr = 'C[C@H]1CCC(CC1)C'
    #testStr = '[F]/C=C/F'
    #testStr = 'C(\F)=C/F'
    #testStr = r'C/1=C/C=C\C=C/C=C\1'
    #testStr = 'C\\1=CCCCCCC/1'
    
    # NONSTANDARD Examples compatible to OpenSMILES standard
    #testStr = 'C1.C12.C2' # ring-bonds used to connect dot-separated structures
    #testStr = 'C%00CC%00' # multi-digit single digit ring ID
    
    
    
    # NONSTANDARD Examples that are still parsed
    #testStr = 'C(C(C.C)C)C' # dot within branch, parsed as if the part after the dot in the branch is a distinct molecule
    #testStr = 'C(C)1CC1' # ring-bond appearing not directly after atom
    #testStr = 'C(.C)C' # dot withing branch
    #testStr = 'C()' # empty branch
    #testStr = 'C(C())C' # empty branch
    #testStr = 'C((C))C' # branch not following atom, additional parentheses
    #testStr = 'C((C)C)C' # branch not following atom
    #testStr = '.C' # starting with dot
    #testStr = 'C..C' # multiple dots
    #testStr = 'C.' # trailing dots
    #testStr = 'C=(O)C' # separated bond and atom. The atom that appears first is counted first when chirality is considered
    #testStr = 'C=#1CCC1' # ring not following atom. The atom that appears first is counted first when chirality is considered 
    #testStr = 'C(1CC)CC1' # ring bond in branch
    #testStr = 'C1CC(=1)' # ring bond in branch
    #testStr = 'C1CC(1)' # ring bond in branch
    #testStr = 'C(C.)(C)C' # trailing dot in branch
    #testStr = 'C==CC' # repeated bond symbols are ignored
    #testStr = 'C(1CC1)' # ring bond in branch
    #testStr = 'C(1)CC1' # ring bond in branch
    
    
    # ILLEGAL Examples
    # ring errors
    #testStr = 'C=1CCCC-1' # nonmatching symbols
    #testStr = 'C12CCCCC12' # multiple bond between two atoms
    #testStr = 'C12C2CCC1' # multiple bond between atoms
    #testStr = 'C11' # self-loop
    #testStr = '1CC' # ring not following atom
    #testStr = 'C1CC2CC' # unclosed rings
    #testStr = 'C.1CCCCC.1' # dot cannot appear in front of ring bond
    #testStr = r'C/1=C/C=C\C=C/C=C/1' # inconsistent cis-trans ring-bond
    
    # dot errors
    #testStr = 'C.1CCCCC.1' # dot cannot appear in front of ring bond
    
    # bond error
    #testStr = '=CCCC' # bond not following an atom
    #testStr = 'C(CCCC=)' # bond not followed by an atom
    #testStr = 'C=-CC' # unmatching repeated bond
    
    # branch error
    #testStr = '(CO)=O' # branch not following atom, ambiguous connection
    #testStr = 'C.(C)' # branch not following atom
    #testStr = 'C(=(C))C' # branch not following atom
    #testStr = 'C(C(C)))C' # unbalanced parenthesizes
    #testStr = 'C(C(C)C' # unbalanced parenthesizes
    
    # chiral center error
    #testStr = 'C[C@H2]1CCC(CC1)C' # too many H
    #testStr = 'C[C@H]([H])CCC(CC)C' # too many H
    #testStr = 'C[C@]([2H])([2H])CCC(CC)C' # too many D
    #testStr = 'C[C@](C)CCC(CC)C' # too few neighbors
    #testStr = 'C[C@](C)(CC)([H])C' # too many neighbors
    
    # cis-trans error
    #testStr = r'C\C(/C)=C' # conflicting up's 
    #testStr = r'C\C(\C)=C(/C)/C' # conflicting down's 
    #testStr = r'C\C(\C)=C/C\C' # \ or / not adjacent to double bonds
    #testStr = r'C\C(\C)=C(/C)(C)C' # too many neighbors
    
    # other error
    #testStr = 'C?C' # undefined characters
    #testStr = 'C\C=C(/Br)\C(O)=2C(Cl)(I)O2' # valence error
    #testStr = 'CC([4H])1CCC(CC1)C' # erroneous isotope
    # valency only checked for carbons with specified chirality and allenic carbons 

    
    #s = testStr.split()
    testSMILES = SMILES(testStr)
    G = testSMILES.parse()
    s,G,T = testSMILES.write()
    print(s)
    #H = disjoint_union(G,G)
    #testSMILES.draw()

    