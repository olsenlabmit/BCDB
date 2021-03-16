# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 16:42:29 2019

@author: tslin
"""


from pyparsing import Literal,Word, White, alphas, nestedExpr, quotedString, cStyleComment, alphanums, nums, StringStart, StringEnd

    
class SmilesPattern:
    def __init__(self):
        pass
    
    def addRawStr(toks):
        if 'branch' in toks:
            toks['rawStr']=toks['branch']
        else:
            toks['rawStr']=''.join(toks[:])
        return toks
    
    #whitespace = " \t\n"
    whitespace = White().leaveWhitespace()
    ### ATOM SECTION ### 
    # Organic Subset section
    _aliphatic_organic = ( Literal('Cl').setResultsName('symbol') \
                         | Literal('Br').setResultsName('symbol') \
                         | Word('BCNOSPFI',exact=1).setResultsName('symbol') ).setResultsName('organic')
    _aromatic_organic = ( Literal('c').setResultsName('symbol') \
                         | Word('bnosp',exact=1).setResultsName('symbol') ).setResultsName('organic')
    #_aliphatic_organic.setResultsName('organic')
    #_aromatic_organic.setResultsName('organic')
    
    # Bracketed Atoms section
    _isotope = Word(nums,min=1)
    _element_symbols =Literal('He') | Literal('Li') | Literal('Be') | Literal('Ne') | Literal('Na') | Literal('Mg') \
                    | Literal('Al') | Literal('Si') | Literal('Cl') | Literal('Ar') | Literal('Ca') | Literal('Sc') \
                    | Literal('Ti') | Literal('Cr') | Literal('Mn') | Literal('Fe') | Literal('Co') | Literal('Ni') \
                    | Literal('Cu') | Literal('Zn') | Literal('Ga') | Literal('Ge') | Literal('As') | Literal('Se') \
                    | Literal('Br') | Literal('Kr') | Literal('Rb') | Literal('Sr') | Literal('Zr') | Literal('Nb') \
                    | Literal('Mo') | Literal('Tc') | Literal('Ru') | Literal('Rh') | Literal('Pd') | Literal('Ag') \
                    | Literal('Cd') | Literal('In') | Literal('Sn') | Literal('Sb') | Literal('Te') | Literal('Xe') \
                    | Literal('Cs') | Literal('Ba') | Literal('Hf') | Literal('Ta') | Literal('Re') | Literal('Os') \
                    | Literal('Ir') | Literal('Pt') | Literal('Au') | Literal('Hg') | Literal('Tl') | Literal('Pb') \
                    | Literal('Bi') | Literal('Po') | Literal('At') | Literal('Rn') | Literal('Fr') | Literal('Ra') \
                    | Literal('Rf') | Literal('Db') | Literal('Sg') | Literal('Bh') | Literal('Hs') | Literal('Mt') \
                    | Literal('Ds') | Literal('Rg') | Literal('La') | Literal('Ce') | Literal('Pr') | Literal('Nd') \
                    | Literal('Pm') | Literal('Sm') | Literal('Eu') | Literal('Gd') | Literal('Tb') | Literal('Dy') \
                    | Literal('Ho') | Literal('Er') | Literal('Tm') | Literal('Yb') | Literal('Lu') | Literal('Ac') \
                    | Literal('Th') | Literal('Pa') | Literal('Np') | Literal('Pu') | Literal('Am') | Literal('Cm') \
                    | Literal('Bk') | Literal('Cf') | Literal('Es') | Literal('Fm') | Literal('Md') | Literal('No') \
                    | Literal('Lr') \
                    | Literal('H') | Literal('B') | Literal('C') | Literal('N') | Literal('O') | Literal('F') | Literal('P')  \
                    | Literal('S') | Literal('K') | Literal('V') | Literal('Y') | Literal('I') | Literal('W') | Literal('U')
    _aromatic_symbols = Literal('se') | Literal('as') | Word('cnops',exact=1) 
    _symbol = _element_symbols | _aromatic_symbols | Literal('*')
    
    # Chirality section
    _chiral = Literal('@@') | Literal('@')  #|  Literal('@TH1') | Literal('@TH2') \
            #| Literal('@SP1') | Literal('@SP2') | Literal('@SP3') \
            #| Literal('@AL1') | Literal('@AL2') | '@TB'+Word(nums,min=1,max=2) | '@OH'+Word(nums,min=1,max=2) 
    _chiral.setParseAction(''.join)
    
    # Hydrogens section
    _hcount = Literal('H') + (Word('123456789',exact=1)*(0,1)).setResultsName('nH')
    #_hcount.setParseAction(''.join)
    
    # Charge section
    _charge = ('-' + Word('123456789',exact=1)*(0,1)) | ('+' + Word('123456789',exact=1)*(0,1)) | Literal('--') | Literal('++')
    #_charge.setParseAction(''.join)
    
    # Atom Class section
    _class = ':' + Word(nums,min=1)
    
    # Bracketed Atom definition
    _bracket_atom = '[' + _isotope.setResultsName('isotope')*(0,1)  \
                        + _symbol.setResultsName('symbol')          \
                        + _chiral.setResultsName('chiral')*(0,1)    \
                        + _hcount.setResultsName('hcount')*(0,1)    \
                        + _charge.setResultsName('charge')*(0,1)    \
                        + _class.setResultsName('_class')*(0,1)      \
                        + ']'
    #_bracket_atom.setResultsName('bracket_atom')
    
    # Atom definition
    #_atom = _aliphatic_organic | _aromatic_organic | _bracket_atom | Literal('*').setResultsName('symbol')
    _atom = _aliphatic_organic \
            | _aromatic_organic \
            | _bracket_atom.setResultsName('bracket_atom') \
            | Literal('*').setResultsName('symbol')
    #def addRawStr(toks):
    #    toks['rawStr']=''.join(toks)
    #    return toks
    #_atom.setParseAction(addRawStr)
    _atom.leaveWhitespace()
    #_atom.setParseAction(''.join)
    #_atom.setParseAction(lambda locn,tokens: (locn,''.join(tokens[:])))
    
    ### BOND SECTION ###
    _bond = Word('-=#:\/',exact=1)
    _bond.leaveWhitespace()
    #_bond.setParseAction(addRawStr)
    
    #_ringbond = _bond*(0,1) + \
    #            (Word(nums,exact=1).setParseAction(lambda tok:[''.join(tok)] ) | \
    #            (Literal('%')+Word(nums,exact=2).setResultsName('ringid')).setParseAction(lambda tok:[''.join(tok[:])] ) )
    _ringbond = (_bond*(0,1)).setResultsName('ringbondtype') + \
                (Word(nums,exact=1).setResultsName('ringid') | \
                 Literal('%')+Word(nums,exact=2).setResultsName('ringid') )
    
    _ringbond.leaveWhitespace()
    #_ringbond.setParseAction(addRawStr)
    
    _dot = Literal('.')
    #_dot.setParseAction(addRawStr)
    
    _smilesChar = _ringbond.setResultsName('ringbond') | _bond.setResultsName('bond') \
                | _atom.setResultsName('atom') | _dot.setResultsName('dot')
    
    _branchContent = _smilesChar*(1,None)
    _branchContent.setParseAction(lambda toks: ''.join(toks))
    
    
    _branch = nestedExpr('(',')',content=_branchContent)
    _branch.setParseAction(lambda toks: '('+''.join([str(item) for sublist in toks for item in sublist])+')')
    
    
    _smilesElement = _smilesChar | _branch.setResultsName('branch')
    _smilesElement.setParseAction(addRawStr)
    
    #_ringbond.setParseAction(lambda tok:tok[0],tok[1])
    
    
    #_bonds = (_bond | _ringbond)*(0,1)
    #_bonds.setParseAction(''.join)

if __name__ == "__main__":
        
    testStr = '[12C]/[OH](C(CC(C)C)C)[CH2](C)(C)CC'
    testStr = 'F/C=C/F'
    
    pattern = SmilesPattern._smilesElement
    res = pattern.parseString(testStr)
    if 'dot' in res.keys():
        print('dot')
    elif 'branch' in res.keys():
        print('branch')
    elif 'bond' in res.keys():
        print('bond')
    elif 'ringbond' in res.keys():
        print('ringbond')
    else:
        print('atom')
    print(res.rawStr)