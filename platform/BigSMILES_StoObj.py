# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 15:43:33 2019

@author: tslin
"""

from BigSmilesPattern import BigSmilesPattern
#from BigSmiles import BigSMILES
from utility import errorMsg
from BigSMILES_Bond import BigSMILES_Bond
from error import BigSMILES_StoObjError,BigSMILES_BondInconsistencyError
#import networkx as nx


class BigSMILES_StoObj:
    _count = 0

    
    def __init__(self,inStr,pos=0,index=list(),prevBond=None):
        # get unique id associated with the object
        self._uid = BigSMILES_StoObj._count
        self.index = index
        self.tmpHolder = 0
        self.pos = pos
        BigSMILES_StoObj._count = BigSMILES_StoObj._count + 1
        # Save the string to raw string
        self.rawStr = inStr
        # the repUnit and endGrp are containers for the repeating units and end groups in the stochastic object
        # they are lists of BigSMILESChain 
        self.repUnit = []
        self.endGrp = []
        # leftEnd and rightEnd are BigSMILESBond objects that holds the bonding descriptor object 
        # of the two ends of the stochastic object (the [$1] or [<=1] descriptors seen in bigsmiles strings)
        self.leftEnd = None
        self.rightEnd = None
        # the Bond_Dict dictionary holds the BigSMILES bonds that had already appeared. 
        # each key correspond to the [$<>][id] key
        # each value is a list with exactly one element, which is a copy of the first BigSMILES_Bond object that was declared
        # this dictionary needs to be updated after every repeat unit or end group is parsed
        self.Bond_Dict = dict()
        
        #print(inStr)
        # parse the left terminal bonding descriptor
        peeledStr = self.rawStr[1:-1]
        try:
            peeledStr,openerLen = self.parseLeftTerminal(peeledStr,prevBond)
        except:
            errorMsg(self.rawStr,1,'Error','unrecognized element encountered during parsing the left bonding descriptor')
            raise BigSMILES_StoObjError
        # parse the right terminal bonding descriptor
        try:
            peeledStr = self.parseRightTerminal(peeledStr)
        except BigSMILES_BondInconsistencyError:
            errorMsg(self.rawStr,len(self.rawStr)-self.tmpHolder,'Error','Inconsistency in the right-end bonding descriptor')
            raise BigSMILES_StoObjError
        except:
            errorMsg(self.rawStr,len(self.rawStr)-1,'Error','unrecognized element encountered during parsing the right bonding descriptor')
            raise BigSMILES_StoObjError
        
        
        
        # Parse the peeledStr to get the two lists
        try:
            res = BigSmilesPattern._StoObjLists.parseString(peeledStr)
        except:
            errorMsg(inStr,openerLen+1,'Error','unrecognized element encountered during parsing of Stochastic Object')
            raise BigSMILES_StoObjError
        
        from BigSMILES_BigSmilesObj import BigSMILES
        
        localpos = openerLen+1
        for i in range(len(res.repUnit)):
            #print(repUnit)
            repUnit = res.repUnit[i]
            try:
                self.repUnit.append(BigSMILES(inStr=repUnit,pos=localpos+pos,UpperBond_Dict=self.Bond_Dict,index=self.index+[i+1]))
            except:
                errorMsg(inStr,localpos,'Error','error in parsing repeat unit ['+str(len(self.repUnit)+1)+']' )
                raise BigSMILES_StoObjError 
            else:
                # update Bond_Dict
                for key in self.repUnit[-1].Bond_Dict.keys():
                    if key in self.Bond_Dict:
                        pass
                    else:
                        self.Bond_Dict[key] = list(self.repUnit[-1].Bond_Dict[key][0:1])
                localpos += len(self.repUnit[-1].rawStr)+1
            #print(self.repUnit)
            #print(repUnit)
            
        
        #print('here2')
        #print(res.endGrp)
        for i in range(len(res.endGrp)):
            endGrp = res.endGrp[i]
            try:
                self.endGrp.append(BigSMILES(inStr=endGrp,pos=localpos+pos,UpperBond_Dict=self.Bond_Dict,index=self.index+[len(res.repUnit)+i+1]))
            except:
                errorMsg(inStr,localpos,'Error','error in parsing end group ['+str(len(self.endGrp)+1) +']' )
                raise BigSMILES_StoObjError 
            else:
                # update Bond_Dict
                for key in self.repUnit[-1].Bond_Dict.keys():
                    if key in self.Bond_Dict:
                        pass
                    else:
                        self.Bond_Dict[key] = list(self.repUnit[-1].Bond_Dict[key][0:1])
                localpos += len(self.endGrp[-1].rawStr)+1
        #print('here3')
        if len(res.rawStr) < len(peeledStr):
            errorMsg(inStr,openerLen+1+len(res.rawStr),'Error','unrecognized element encountered during parsing of Stochastic Object')
            raise BigSMILES_StoObjError
        #print('here4')
        
        # update the right terminal bonding descriptor after entire object is parsed
        if self.rightEnd:
            self.rightEnd.checkConsistency(self.Bond_Dict,item=None)
        
        
        return None

    def __str__(self):
        #return self.rawStr
        #return '{' + ','.join([str(x) for x in self.repUnit]) + ';' + ','.join([str(x) for x in self.endGrp]) + '}'
        string = '{' 
        if self.leftEnd == None:
            string = string + '[]'
        else:
            string = string + '[' + self.leftEnd.getCompKey() + ']'
        for obj in self.repUnit:
            string = string + str(obj) + ','
        
        string = string[:-1]+';'
        
        for obj in self.endGrp:
            string = string + str(obj) + ','
        
        string = string[:-1] 
        if self.rightEnd == None:
            string = string + '[]'
        else:
            string = string + '[' + self.rightEnd.getCompKey() + ']'
        string = string + '}'
        return string
        
        
        
    def __len__(self):
        return len(self.repUnit)+len(self.endGrp)
    
    def __getitem__(self,key):
        if key <= len(self.repUnit) and key >= 0:
#            return self.repUnit[key-1]
            return self.repUnit[key]
        elif key <= self.__len__():
#            return self.endGrp[key-1-len(self.repUnit)]
            return self.endGrp[key-len(self.repUnit)]
        else:
            return None
    
    def __iter__(self):
        self.n = 0
        return self
    
    def __next__(self):
        if self.n < len(self.repUnit):
            result = self.repUnit[self.n]
            self.n += 1
            return result
        elif self.n < len(self.repUnit)+len(self.endGrp):
            result = self.endGrp[self.n-len(self.repUnit)]
            self.n += 1
            return result
        else:
            raise StopIteration
        
        
    def parseLeftTerminal(self,peeledStr,prevBond):
        # peeledStr is the string to be parsed (the raw string with the outer {} stripped)
        
        # parse the left terminal
        try:
            res = BigSmilesPattern._opener.parseString(peeledStr)
        except:
            #errorMsg(self.rawStr,1,'Error','unrecognized element encountered during parsing the left bonding descriptor')
            raise BigSMILES_StoObjError        
        
        if 'BigSMILES_Bondtype' not in res.keys():
            self.leftEnd = None
        else:
            self.leftEnd = BigSMILES_Bond(res,S_bond=prevBond,Bond_Dict=self.Bond_Dict,item=None)
        
        peeledStr = peeledStr[len(res.rawStr):]
        
        return peeledStr,len(res.rawStr)
    
    def parseRightTerminal(self,peeledStr):
        # parse the right terminal
        try:
            res,start,end = BigSmilesPattern._closer.scanString(peeledStr).__next__()
        except:
            errorMsg(self.rawStr,len(self.rawStr)-1,'Error','unrecognized element encountered during parsing the right bonding descriptor')
            raise BigSMILES_StoObjError
        if 'BigSMILES_Bondtype' not in res.keys():
            self.rightEnd = None
        else:
            try:
                self.rightEnd = BigSMILES_Bond(res,S_bond='u',Bond_Dict=self.Bond_Dict,item=None)
            except:
                self.tmpHolder = len(res.rawStr)
                raise BigSMILES_BondInconsistencyError
                
        if len(res.rawStr) > 0:
            peeledStr = peeledStr[:-len(res.rawStr)]
        else:
            peeledStr = peeledStr[:]

        return peeledStr
    
    def getBond(self,end='left'):
        if end == 'left':
            tmpBond = self.leftEnd
        elif end == 'right':
            tmpBond = self.rightEnd
        else:
            return None
        
        if tmpBond == None:
            return None
        else:
            return tmpBond.getS_bond()
        

        
    
if __name__ == "__main__":
        
    #print(SmilesPattern._bond)
    #print(bigSmilesPattern._bondDesc)
    
    testStr = '$[$-1]1'
    testStr = '$=1'
    testStr = 'C<CO>{[$1]$C{[>]<CCO>[<]}CC$,$CC$;$C[$]}CC(C)C{[>]<CCO>[<]}{CC(Cl);CC}C'
    testStr = '{[$]$CC$[>]}'
    
    P = BigSMILES_StoObj(testStr,pos=0)
    #P.setBigSMILESstring(testStr)
    #P.parse()    