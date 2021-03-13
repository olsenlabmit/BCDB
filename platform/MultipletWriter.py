#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 16:36:49 2019

@author: Cshi
"""
import pandas as pd
from BigSMILES_BigSmilesObj import BigSMILES

df = pd.read_csv("/Users/Cshi/Desktop/BigSMILES/BigSMILESCodes/new_bond_desc/Tg_dataset_newSyntax.csv")
tg = df["Tg"]
index = 0
for name in df["BigSMILES"]:
    polymer = BigSMILES(name)
    dic = polymer[0][0].multiplet(3) #change this parameter to get different multiplets
    df = pd.DataFrame.from_dict(dic, orient='index', columns = [name])
    if index == 0:
        df1 = df
    else:
        df1= pd.concat([df1,df],sort=False, axis =1)
    index += 1
df2 = df1.transpose()
df2.insert(0,"Tg", list(tg))
df2.to_csv('triplets.csv')
    





