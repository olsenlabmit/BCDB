# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:02:36 2019

@author: tslin
"""

from SmilesPattern import SmilesPattern
from pyparsing import Word, Literal, nestedExpr, nums, StringStart, StringEnd, printables,alphanums
from utility import flatten_list
  

class BigSmilesPattern(SmilesPattern):
    
#### DEFINITIONS of patterns involved in BigSMILES_Bond ####
    _BigSmilesBondChar = "$<>"
    _BondDesc =  Word(_BigSmilesBondChar,exact=1).setResultsName('BigSMILES_Bondtype') + \
                ( (Word(nums,exact=1).setResultsName('BigSMILES_Bondid') | \
                  Literal('%')+Word(nums,exact=2).setResultsName('BigSMILES_Bondid') ) )*(0,1) 
#    _ladderBondDesc = Word(_bigsmilesBondChar,exact=1).setResultsName('BigSMILES_outerBondtype') + \
#                      '[' + _bondDesc + ']' + \
#                      (Word(nums,exact=1).setResultsName('BigSMILES_outerbondid') | \
#                      Literal('%')+Word(nums,exact=2).setResultsName('BigSMILES_outerbondid') ) 
                      
#    _bigsmilesBond = _ladderBondDesc.setResultsName('BigSMILES_ladderBond') | _bondDesc.setResultsName('BigSMILES_Bond')
    _BigSmilesBond = (Literal('[') + _BondDesc.setResultsName('BigSMILES_Bond') + Literal(']') )

#### DEFINITIONS of patterns involved in Augmented_SMILES ####
    
    # redefinition for the elements used in parsing of Augmented SMILES strings 
    _AugmentedSmilesChar = SmilesPattern._smilesChar | _BigSmilesBond
    
    _AugmentedBranchContent = _AugmentedSmilesChar*(1,None)
    _AugmentedBranchContent.setParseAction(lambda toks: ''.join(toks))
    
    _AugmentedBranch = nestedExpr('(',')',content=_AugmentedBranchContent)
    #_AugmentedBranch.setParseAction(lambda toks: '('+''.join([str(item) for sublist in toks for item in sublist])+')')
    _AugmentedBranch.setParseAction(lambda toks: '('+''.join(flatten_list(toks,str))+')')
    
    # _AugmentedSmilesElement explicitly used in Augmented_SMILES()
    _AugmentedSmilesElement = _AugmentedSmilesChar | _AugmentedBranch.setResultsName('branch')
    _AugmentedSmilesElement.addParseAction(SmilesPattern.addRawStr)

    
    
#### DEFINITIONS of stochastic object
    _StoObjSepChar = ",;"
    #_BracketedBond = (Literal('[') + _BigSmilesBond + Literal(']')).setResultsName('BigSMILES_bracketedBond')
    _TerminalBond = (Literal('[') + (_BondDesc*(0,1)).setResultsName('BigSMILES_Bond') + Literal(']')).setResultsName('BigSMILES_terminalBond')
    _opener = StringStart() + _TerminalBond
    _opener.setParseAction(SmilesPattern.addRawStr)
    _closer = _TerminalBond + StringEnd()
    _closer.setParseAction(SmilesPattern.addRawStr)
    _StoObjSep = Word(',;',exact=1)
    
    printableExceptCurly = printables.replace('{', '').replace('}', '')
    _StoObjContent = Word(printableExceptCurly)#.setResultsName('StoObjCont')
    _StoObjContent.setParseAction(lambda toks: ''.join(toks))
    _StoObj = nestedExpr('{','}',content=_StoObjContent)
    _StoObj.setParseAction(lambda toks: '{'+''.join(flatten_list(toks,str))+'}')
    #_StoObjDummy = (Literal('{') + Word(nums).setResultsName('StoObjId') + Literal('}')).setResultsName('BigSMILES_StoObj')
    
    
    def separateList(toks):
        L = [x for x in toks if x != ',']
        if not 'endGrp' in toks.keys():
            toks['repUnit'] = L
            toks['endGrp'] = list()
        else:
            n = L.index(';')
            toks['repUnit'] = L[:n]
            toks['endGrp'] = L[n+1:]
        toks['rawStr'] = ''.join(toks)
        return toks
    
    printableExceptSemicolon = printables.replace(';', '')
    printableExceptCommaSemicolon = printableExceptSemicolon.replace(',', '')
    _StoObjUnit = Word(printableExceptCommaSemicolon)
    _StoObjList = _StoObjUnit + ("," + _StoObjUnit)*(0,None)
    _StoObjLists = _StoObjList.setResultsName('repUnit') + \
                   (Literal(';') + _StoObjList.setResultsName('endGrp'))*(0,1)
    _StoObjLists.setParseAction(separateList)
    
    
    
#### DEFINITIONS of patterns involved in BigSMILES() ####
    _BigSmilesChar = SmilesPattern._smilesChar | _BigSmilesBond | _StoObj.setResultsName('BigSMILES_StoObj')
    
    _BigSmilesBranchContent = _BigSmilesChar*(1,None)
    _BigSmilesBranchContent.setParseAction(lambda toks: ''.join(toks))
    
    _BigSmilesBranch = nestedExpr('(',')',content=_BigSmilesBranchContent)
    _BigSmilesBranch.setParseAction(lambda toks: '('+''.join(flatten_list(toks,str))+')')
    
    _BigSmilesElement = _BigSmilesChar | _BigSmilesBranch.setResultsName('branch')
    _BigSmilesElement.addParseAction(SmilesPattern.addRawStr)
    
    

    # additional definition of augmented SMILES for parsing the entire augmented SMILES segments (BigSMILES Chain Objects)
#    _augBranchContent = (SmilesPattern._smilesChar | _bigsmilesBond)*(1,None)
#    _augBranch = nestedExpr('(',')',content=_augBranchContent)
#    _augBranch.setParseAction(lambda toks: '('+''.join([str(item) for sublist in toks for item in sublist])+')')
#    _augSmilesElement = SmilesPattern._smilesChar | _augBranch 
#    _bigsmilesChar = _augSmilesElement | _bigsmilesBond
                     
#    _bigsmileschainObj = _bigsmilesChar*(1,None)
#    _bigsmileschainObj.setParseAction(lambda toks: ''.join(toks))

    
    # bracketed bonds and definition for starting/ending patterns for stochastic objects
    # _opener, _closer, _stoObjExactContent explicitly used in BigSMILES_StoObj()
#    _bracketedBond = (Literal('[') + _bigsmilesBond + Literal(']')).setResultsName('BigSMILES_bracketedBond')
#    _opener = StringStart()+(_bracketedBond)*(0,1)
#    _opener.setParseAction(SmilesPattern.addRawStr)
#    _closer = (_bracketedBond)*(0,1)+StringEnd()
#    _closer.setParseAction(SmilesPattern.addRawStr)
#    _stoObjSep = Word(',;',exact=1)
    
#    _stoObjContent = (_AugmentedSmilesChar | _stoObjSep | _bracketedBond)*(1,None)
#    _stoObjContent.setParseAction(lambda toks: ''.join(toks))
    
#    _stoObj = nestedExpr('{','}',content=_stoObjContent)
    #_stoObj.setParseAction(lambda toks: '{'+''.join([str(item) for sublist in toks for item in sublist])+'}')
#    _stoObj.setParseAction(lambda toks: '{'+''.join(flatten_list(toks,str))+'}')
    
#    _bigsmilesElement = _bigsmileschainObj.setResultsName('Augmented_SMILES') | \
#                        _stoObj.setResultsName('BigSMILES_StoObj')
#    _bigsmilesElement.setParseAction(SmilesPattern.addRawStr)
    
    def separateList(toks):
        L = [x for x in toks if x != ',']
        if not 'endGrp' in toks.keys():
            toks['repUnit'] = L
            toks['endGrp'] = list()
        else:
            n = L.index(';')
            toks['repUnit'] = L[:n]
            toks['endGrp'] = L[n+1:]
        toks['rawStr'] = ''.join(toks)
        return toks
    
#    _stoObjUnit = (_stoObj | _bigsmileschainObj)*(1,None)
#    _stoObjUnit.setParseAction(''.join)
#    _stoObjList = _stoObjUnit + ("," + _stoObjUnit)*(0,None)
    #_stoObjList.setParseAction(removeDelimiter)
    #_stoObjList.setParseAction(lambda toks: [x for x in toks if x != ','])
#    _stoObjExactContent = _stoObjList.setResultsName('repUnit') + \
#                          (Literal(';') + _stoObjList.setResultsName('endGrp'))*(0,1)
#    _stoObjExactContent.setParseAction(separateList)
    
    #_smilesElement = _smilesChar | _branch.setResultsName('branch')
    #_smilesElement.setParseAction(addRawStr)
    
#    _branch = nestedExpr('(',')',content=_branchContent)
#    _branch.setParseAction(lambda toks: '('+''.join([str(item) for sublist in toks for item in sublist])+')')
    
    
if __name__ == "__main__":
        
    #print(SmilesPattern._bond)
    #print(bigSmilesPattern._bondDesc)
    
    testStr = '$[$-1]1'
    testStr = '$=1'
    testStr = 'C(C<)O>{[$1]$C{[>]<CCO>[<]}CC$,$CC$;$C[$]}{[>]<CCO>[<]}'
    testStr = '[$1]$CC$[$2]'
    testStr = 'CC{1{2}2}{2}C{3}CC'
    testStr = '<C,CCC,C[1]{123},OCC;CCCC,CC'
    testStr = '[$]X'
    #testStr = '$CC$[$[<=1]2]'
    #testStr = '$CC$'
    #testStr = '$CC,{[>]C(C<)CO>[<]}$,$C(C)C(C)C$,$CC;$C,$O'
    
    pattern = BigSmilesPattern._BigSmilesBond
    #results = pattern.scanString(testStr)
    #for res,start,end in results:
    #    print(start)
    #    print(end)
    #    print(res)
    #res,start,end = pattern.scanString(testStr).__next__()
    #print(res.rawStr)
    res = pattern.parseString(testStr)