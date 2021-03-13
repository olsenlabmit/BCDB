# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 14:08:58 2019

@author: tslin
"""
from error import BigSMILES_BondInconsistencyError

class BigSMILES_Bond:
    
    def __init__(self,res,S_bond,Bond_Dict=None,item=None):
        if 'BigSMILES_Bond' in res.keys():
#            self._isLadder = False
            self.res = res
            self._rawStr = res.rawStr
            self._bondtype = res.BigSMILES_Bondtype
            
            if 'BigSMILES_Bondid' not in res.keys():
                self._id = '1'
            else:
                if res.BigSMILES_Bondid == '':
                    self._id = '1'
                else:                
                    self._id = res.BigSMILES_Bondid
            
            self.S_bond = ''
            
            self.setS_bond(S_bond)
            flag = self.checkConsistency(Bond_Dict,item)
            if flag == -1:
                raise BigSMILES_BondInconsistencyError
        else:
            self.noBond = True


    def setS_bond(self,S_bond):
        if S_bond == None:
            self.S_bond = '-'
        else:
            if S_bond == '\\' or S_bond == '/' or S_bond == '-':
                self.S_bond = '-'
            elif S_bond == '=':
                self.S_bond = '='
            elif S_bond == ':':
                self.S_bond = ':'
            elif S_bond == '#':
                self.S_bond = '#'
            elif S_bond == '$':
                self.S_bond = '$'
            elif S_bond == 'u':
                self.S_bond = 'u'
            else:
                self.S_bond = None
            
    def getS_bond(self):
        return self.S_bond
    
    def checkConsistency(self,Bond_Dict,item):
        # check consistency with Bond_Dict entries
        if Bond_Dict != None:
            #print('checking bond consistency')
            key = self.getCompKey()
            #print(key)
            if key in Bond_Dict:
                if self.S_bond == 'u':
                    if Bond_Dict[key][0].S_bond != 'u':
                        self.S_bond = Bond_Dict[key][0].S_bond
                    else:
                        pass
                else:
                    if Bond_Dict[key][0].S_bond == 'u':
                        Bond_Dict[key][0].S_bond = self.S_bond
                    else:
                        if not self.compare(Bond_Dict[key][0]):
        #                    raise BigSMILES_BondInconsistencyError
                            return -1
                        else:
                            pass
            else:
                Bond_Dict[key] = list()                  
                Bond_Dict[key].append(BigSMILES_Bond(self.res,self.S_bond))
            
            if item != None:
                if item not in Bond_Dict[key]:
                    Bond_Dict[key].append(item)
            
            return 1
        else:
            return 0

    def getCompKey(self):
        return self._bondtype + str(self._id)
        
    # check if the bigSmiles bond is consistent with B
    # if consistent, return True
    # if not, return False
    # (TODO?) in addition, update self if  
    def compare(self,B):
        selfOrder = self.getBondOrder()
        BOrder = B.getBondOrder()
        
        if selfOrder == 0:
            return True
        elif selfOrder == BOrder:
            return True
        elif BOrder ==0 and selfOrder == 1:
            return True
        else:
            return False
            
        
        
    def getBondOrder(self,S_bond=None):
        if S_bond == None:
            S_bond = self.S_bond
        if S_bond == '-' or S_bond == '/' or S_bond == '\\':
            order = 1
        elif S_bond == ':':
            order = 1.5
        elif S_bond == '=':
            order = 2
        elif S_bond == '#':
            order = 3
        elif S_bond == 'u':
            order = -1
        else:
            order = 0
        return order
    
        
    def getCompleteSymbol(self):
#        return self._bondtype + self.S_bond + self._id
        return self._bondtype + self._id
    
    def __str__(self):
        return self._rawStr