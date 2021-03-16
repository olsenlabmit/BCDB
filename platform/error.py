# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 18:24:02 2019

@author: tslin
"""

class Error(Exception):
    pass

class SmilesError(Error):
    def __init__(self, message=''):
        self.message = message

class NonmatchingSymbolError(SmilesError):
    def __init__(self):
        pass

class NonmatchingTransCisError(SmilesError):
    def __init__(self):
        pass

class SelfLoopError(SmilesError):
    def __init__(self):
        pass

class MultiBondError(SmilesError):
    def __init__(self):
        pass

class HcountChiralError(SmilesError):
    def __init__(self):
        pass

class AtomicWeightError(SmilesError):
    def __init__(self):
        pass
    
class bigSmilesError(Error):
    def __init__(self, message=''):
        self.message = message
        
class Augmented_SMILESError(bigSmilesError):
    def __init__(self):
        pass
    
class BigSMILES_StoObjError(bigSmilesError):
    def __init__(self):
        pass

class BigSMILES_BondInconsistencyError(bigSmilesError):
    def __init__(self):
        pass
    
class BigSMILES_StoObjMissingTerminalError(bigSmilesError):
    def __init__(self):
        pass
        
class BigSMILESError(bigSmilesError):
    def __init__(self):
        pass
