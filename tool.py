import sqlite3
import pandas as pd
import numpy as np
import copy
from BigSMILES_BigSmilesObj import chemistry_table
from plot import visualize
import re

import tkinter
from tkinter import *
from tkinter import colorchooser

root = Tk()
root.title("BCDB")
root.geometry("900x500")

global default, search
default = (192,192,192)
search = []
chemistries = chemistry_table()    

def subgr_search(queries, bigsmiles, names):
    """
    Execute all queries on the same bigsmiles
    """
    
    import re
    ends = re.split("{|}",bigsmiles)

    import rdkit
    from rdkit import Chem
    list_of_negatives = []
    for q in range(len(queries)):
        subgraph = queries[q][0]
        if subgraph == "":
            continue
        for index in range(len(names)):
            match = False
            if " rep" in subgraph:
                search = subgraph[0:-4]
                exact = chemistries[chemistries.Name==names[index]]['Exact_Hit']
                exact = list(exact)[0]
                isCycle = list(chemistries[chemistries.Name==names[index]]['isCycle'])[0]
                if isCycle:
                    search = search + search
                else:
                    search = search
                for e in exact:
                    target = Chem.MolFromSmiles(e)
                    a = target.HasSubstructMatch(Chem.MolFromSmarts(search))
                    b = target.GetNumAtoms() == Chem.MolFromSmarts(search).GetNumAtoms()
                    if a and b:
                        match = True
                        break
                if not match:
                    list_of_negatives.append([q,index])

            else:
                substr = chemistries[chemistries.Name==names[index]]['Substructure_Hit']
                substr = list(substr)[0]
                for s in substr:
                    string = ends[2 * index] + s + ends[2 * index + 2]
                    target = Chem.MolFromSmiles(string)
                    if target.HasSubstructMatch(Chem.MolFromSmarts(subgraph)):
                        match = True
                        break
                if not match:
                    list_of_negatives.append([q,index])
    return list_of_negatives

def match_criterion(matrix):
    match = False
    for j in range(len(matrix[0])):
        start = matrix[0][j]
        counter = set()
        if start == 1 and len(matrix) == 1:
            return True
        elif start == 1:
            counter.add(j)
            for k in range(1, len(matrix)):
                for l in range(len(matrix[k])):
                    if matrix[k][l] == 1 and l not in counter:
                        counter.add(l)
                        break
            if len(counter) == len(matrix):
               match = True
    return match

def replace(sql, letter):
    sql = sql.replace(" name "," name" + letter + " ")
    sql = sql.replace(" Mn "," Mn" + letter + " ")
    sql = sql.replace(" Mw "," Mw" + letter + " ")
    sql = sql.replace(" D "," D" + letter + " ")
    sql = sql.replace(" N "," N" + letter + " ")
    sql = sql.replace(" f "," f" + letter + " ")
    sql = sql.replace(" w "," w" + letter + " ")
    sql = sql.replace(" p "," p" + letter + " ")
    return sql

def execute_search(all_color, all_size, all_shape, all_overall_sql, all_query, all_target_blocks):
    connection = sqlite3.connect('BCDB.db')
    cursor = connection.cursor()
    matches = []
    for c in range(len(all_color)):
        color =         all_color[c]
        size =          all_size[c]
        shape =         all_shape[c]
        overall_sql =   all_overall_sql[c]
        query =         all_query[c]
        if len(query) == 0:
            query.append(["",""])
        table =         all_target_blocks[c]
    
        matches_table = []
        for target_blocks in table:
            try:
                table_name = ""
                if target_blocks == 2:
                    table_name = "diblock"
                elif target_blocks == 3:
                    table_name = "triblock"
                elif target_blocks == 4:
                    table_name = "tetrablock"
                elif target_blocks == 6:
                    table_name = "hexablock"
                else:
                    table_name = "undecablock"

                subset = cursor.execute("select ind,BigSMILES from " + table_name)
                hits = []
                for s in subset:
                    hits.append([s[0], s[1], 0, np.zeros((len(query), target_blocks)), [None]*target_blocks])
                
                if overall_sql != "":
                    subset = cursor.execute("select ind from " + table_name + " where " + overall_sql)
                    for s in subset:
                        hits[s[0]][2] = 1
                else:
                    for i in range(len(hits)):
                        hits[i][2] = 1

                for i in range(len(query)):
                    for j in range(target_blocks):
                        sql = query[i][1]
                        subset = cursor.execute("select ind,name" + str(j+1) + " from " + table_name)
                        for s in subset:
                            hits[s[0]][4][j] = s[1]
                        if sql == "":
                            subset = cursor.execute("select ind from " + table_name)
                            for s in subset:
                                hits[s[0]][3][i][j] = 1
                        else:
                            sql = replace(" " + sql, str(j+1))
                            subset = cursor.execute("select ind from " + table_name + " where " + sql)
                            for s in subset:
                                hits[s[0]][3][i][j] = 1
                prune = []
                for i in range(len(hits)):
                    if hits[i][2] == 1 and np.count_nonzero(hits[i][3]) != 0 and match_criterion(hits[i][3]):
                        prune.append(hits[i])
                print("SEARCHING TABLE " + table_name)
                print("DONE WITH SQL")
                print("Number of Hits after SQL " + str(len(prune)))
                match = []
                bigsmiles_prev = ""
                list_of_negatives = []
                for i in range(len(prune)):
                    if prune[i][1] != bigsmiles_prev:
                        bigsmiles_prev = prune[i][1]
                        list_of_negatives = subgr_search(query,prune[i][1],prune[i][4])
                            
                    for j in range(len(list_of_negatives)):
                        q = list_of_negatives[j][0]
                        index = list_of_negatives[j][1]
                        prune[i][3][q][index] = 0
                        
                    if match_criterion(prune[i][3]):
                        match.append(int(prune[i][0]))
                print("Number of Hits after Subgraph Search: " + str(len(match)) + "\n")
                matches_table.append([table_name, match, color, size, shape])
            except:
                print('SEARCH FAILED')
                matches_table.append([table_name, [], color, size, shape])
        matches.append(matches_table)
    connection.close()
    return matches    

def svd(download_data = False):
    all_color = []
    all_size = []
    all_shape = []
    all_overall_sql = []
    all_query = []
    all_target_blocks = []
    for c in range(len(search)):
        if search[c][0] == 0:
            try:
                a_c = search[c][1]
                a_size = int(search[c][2].get())
                a_shape = search[c][3].get()
                a_o_sql = search[c][4].get()
                a_q = []
                for i in range(len(search[c][5])):
                    a_q.append([search[c][5][i][0].get(), search[c][5][i][1].get()])
                nblocks = re.split(",| ",search[c][6].get())
                nblocks = [int(i) for i in nblocks]
            except:
                continue
            all_color.append(a_c)
            all_size.append(a_size)
            all_shape.append(a_shape)
            all_overall_sql.append(a_o_sql)
            all_query.append(a_q)
            all_target_blocks.append(nblocks)
    print(all_color)
    print(all_size)
    print(all_shape)
    print(all_overall_sql)
    print(all_query)
    print(all_target_blocks)
    if len(all_color) == len(all_size) == len(all_shape) == len(all_overall_sql) == len(all_query) == len(all_target_blocks):
        matches = execute_search(all_color, all_size, all_shape, all_overall_sql, all_query, all_target_blocks)
    else:
        matches = [[[""]]]
    for c in range(len(search)):
        if search[c][0] == 1:
            try:
                temp_entry = re.split(', | |,', search[c][7].get())
                temp_entry = [float(i) for i in temp_entry]
                fA_entry = [float(search[c][5].get())] * len(temp_entry)
                Mn_entry = [float(search[c][6].get())] * len(temp_entry)
                matches.append([["benchmark", temp_entry, fA_entry, Mn_entry, search[c][1], int(search[c][2].get()), search[c][3].get()]])
            except:
                continue
    if download_data:
        connection = sqlite3.connect('BCDB.db')
        cursor = connection.cursor()
        import csv
        for s in range(len(matches)):
            with open("Search"+str(s+1)+".csv", "w", encoding="UTF-8", newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                for t in range(len(matches[s])):
                    if matches[s][t][0] == "benchmark":
                        continue
                    subset = list(cursor.execute("select * from " + matches[s][t][0]))
                    columns = list(map(lambda x: x[0], cursor.description))
                    csvwriter.writerow(["Search Information"])
                    csvwriter.writerow(["System/Overall Polymer SQL:", all_overall_sql[s]])
                    csvwriter.writerow(["Individual Block Search"])
                    query = all_query[s]
                    for i in range(len(query)):
                        csvwriter.writerow(["Block #" + str(i + 1), "Subgraph:", query[i][0], "SQL:", query[i][1]])
                    csvwriter.writerow([""])
                    csvwriter.writerow(["BCDB HITS IN TABLE " + matches[s][t][0]])
                    csvwriter.writerow(columns)
                    for index in matches[s][t][1]:
                        new_subset = []
                        for i in list(subset[index]):
                            new_subset.append(i)
                        csvwriter.writerow(new_subset)
                    csvwriter.writerow([""])
        print("FINISHED DOWNLOAD")
        connection.close()
    else:
        visualize(matches, default, choose_plots)
 
def add_search(input):
    index = len(search)
    
    def assign_color(f1, index):
        chosen = colorchooser.askcolor()
        search[index][1] = chosen[0]
        color_chosen = Button(f1, text = "Color", bg = chosen[1], command = lambda: assign_color(f1, index)).grid(row = 0, column = 0, sticky=NSEW)

    def delete(f, index):
        search[index][0] = -1
        f.destroy()
        
    if input == 1:
        f = Frame(frame_scrolling_existing)

        f1 = Frame(f)
        color_button = Button(f1, text = "Color", command = lambda: assign_color(f1, index)).grid(row=0, column=0, sticky=NSEW)
        clear = Button(f1, text = "Delete", command = lambda: delete(f, index)).grid(row=0, column=1, sticky=NSEW)
        size = Entry(f1, width = 3, justify=CENTER)
        size.insert(END, '10')
        size.grid(row=1, column=0, sticky=NSEW)
        marker = Entry(f1, width = 3, justify=CENTER)
        marker.insert(END, 'o')
        marker.grid(row=1, column=1, sticky=NSEW)
        f1.pack(side = TOP, fill = X)

        def add_block_query(f2, index, num_blocks):
            global search
            try:
                n = int(num_blocks.get())
            except:
                return
            for i in range(n):
                block = Frame(f2)
                qLabel = Label(block, text = "#" + str(i+1)+" Subgraph: ").grid(row = 0, column = 0, columnspan = 2, sticky=NSEW)
                subgraph = Entry(block, width = 25, justify=CENTER)
                subgraph.grid(row=0, column=2, columnspan = 3, sticky=NSEW)
                qLabel = Label(block, text = " and SQL Condition(s): ").grid(row = 0, column = 5, columnspan = 2, sticky=NSEW)
                sql = Entry(block, width = 76, justify=CENTER)
                sql.grid(row=0, column=7, columnspan = 3, sticky=NSEW)
                block.pack(side = TOP)
                search[index][5].append([subgraph, sql])

        qLabel = Label(f1, text = "System and Overall Polymer SQL Condition(s): ").grid(row = 0, column = 2, sticky=NSEW)
        overall_sql = Entry(f1, width = 80, justify=CENTER)
        overall_sql.grid(row=0, column=3, columnspan = 4, sticky=NSEW)
        qLabel = Label(f1, text = "Individual Block Search: How Many?").grid(row = 1, column = 2, sticky=NSEW)
        num_blocks = Entry(f1, justify=CENTER)
        num_blocks.grid(row = 1, column = 3, sticky=NSEW)
        
        f2 = Frame(f)
        add_button = Button(f1, text = "Add", command = lambda: add_block_query(f2, index, num_blocks))
        add_button.grid(row=1, column=4, sticky=NSEW)
        Label(f1, text = "Target Block Table(s): ").grid(row = 1, column = 5, sticky=NSEW)
        target_blocks = Entry(f1, width = 10, justify=CENTER)
        target_blocks.grid(row = 1, column = 6, sticky=NSEW)
        f2.pack(side = TOP, fill = X)

        Label(f).pack(side=TOP,fill=X)

        f.pack(side = TOP, fill = X)

        search.append([0, default, size, marker, overall_sql, [], target_blocks])

    elif input == 2:
        f = Frame(frame_scrolling_existing) 
        
        color_button = Button(f, text = "Color", command = lambda: assign_color(f, index)).grid(row = 0, column = 0, sticky=NSEW)
        clear = Button(f, text = "Delete" , command = lambda: delete(f, index)).grid(row = 0, column = 1, sticky=NSEW)
        size = Entry(f, width = 3, justify=CENTER)
        size.insert(END, '10')
        size.grid(row=1, column=0, sticky=NSEW)
        marker = Entry(f, width = 3, justify=CENTER)
        marker.insert(END, 'o')
        marker.grid(row=1, column=1, sticky=NSEW)

        qLabel = Label(f, text = "BigSMILES:").grid(row = 0, column = 2, sticky=NSEW)
        bigsmiles_entry = Entry(f, width = 111, justify=CENTER)
        bigsmiles_entry.grid(row=0, column=3, columnspan = 4, sticky=NSEW)

        qLabel = Label(f, text = "f\u1D00: ").grid(row = 1, column = 2, sticky=NSEW)
        fA_entry = Entry(f, justify=CENTER)
        fA_entry.grid(row=1,column=3, sticky=NSEW)
        qLabel = Label(f, text = "M\u2099 (g/mol): ").grid(row = 1, column = 4, sticky=NSEW)
        Mn_entry = Entry(f, justify=CENTER)
        Mn_entry.grid(row = 1, column = 5, sticky=NSEW)

        qLabel = Label(f, text = "T (\u00B0C): ").grid(row = 4, column = 2, sticky=NSEW)
        temp_entry = Entry(f, justify=CENTER)
        temp_entry.grid(row = 4, column = 3, columnspan = 3, sticky=NSEW)

        qLabel = Label(f, text = "phase1: ").grid(row = 5, column = 2, sticky=NSEW)
        phase1 = Entry(f, justify=CENTER)
        phase1.grid(row = 5, column = 3, sticky=NSEW)
        qLabel = Label(f, text = "phase2: ").grid(row = 5, column = 4, sticky=NSEW)
        phase2 = Entry(f, justify=CENTER)
        phase2.grid(row = 5, column = 5, sticky=NSEW)

        Label(f).grid(row = 6, column = 0)
        
        f.pack(side = TOP, fill = X)

        search.append([1, default, size, marker, bigsmiles_entry, fA_entry, Mn_entry, temp_entry, phase1, phase2])
        
    existing_canvas.configure(scrollregion=existing_canvas.bbox("all"))

master = Frame(root).pack(side = TOP, fill = X)
searches = Frame(master)
a = Button(searches, text = "Add Search", command=lambda: add_search(1)).grid(row = 0, column = 0)
b = Button(searches, text = "Add New Data", command=lambda: add_search(2)).grid(row = 0, column = 1)
c = Button(searches, text = "Download", command=lambda: svd(download_data = True)).grid(row = 0, column = 2)
def change_default(searches):
    chosen = colorchooser.askcolor()
    global default
    default = chosen[0]
    d = Button(searches, text = "Default Color", bg = chosen[1], command = lambda: change_default(searches)).grid(row = 0, column = 3)
d = Button(searches, text = "Default Color", command = lambda: change_default(searches)).grid(row = 0, column = 3)
e = Button(searches, text = "Visualize", command=svd).grid(row = 0, column = 4)
searches.pack(side = TOP)

choose_frame = Frame(master)
choose_plots = [IntVar(),IntVar(),IntVar(),IntVar(),IntVar(),IntVar()]
choose_plots[0].set(1)
for i in range(1, 6):
    choose_plots[i].set(0)
c = Checkbutton(choose_frame, text="T vs. f\u1D00", variable=choose_plots[0]).pack(side = LEFT)
c = Checkbutton(choose_frame, text="T vs. M\u2099", variable=choose_plots[1]).pack(side = LEFT)
c = Checkbutton(choose_frame, text="M\u2099 vs. f\u1D00", variable=choose_plots[2]).pack(side = LEFT)
c = Checkbutton(choose_frame, text="M\u2099 hist.", variable=choose_plots[3]).pack(side = LEFT)
c = Checkbutton(choose_frame, text="T hist.", variable=choose_plots[4]).pack(side = LEFT)
c = Checkbutton(choose_frame, text="f\u1D00 hist.", variable=choose_plots[5]).pack(side = LEFT)
choose_frame.pack(side = TOP)

existing_data = Frame(master).pack(side = TOP)
existing_canvas = Canvas(existing_data)
existing_canvas.pack(side = TOP, fill = BOTH, expand = 1)
existing_scroll_y = Scrollbar(existing_canvas, orient=VERTICAL, command=existing_canvas.yview)
existing_scroll_y.pack(side=RIGHT, fill=Y)
existing_canvas.configure(yscrollcommand=existing_scroll_y.set)
existing_canvas.bind('<Configure>', lambda e: existing_canvas.configure(scrollregion = existing_canvas.bbox("all")))
frame_scrolling_existing = Frame(existing_canvas)
existing_canvas.create_window((0,0), window=frame_scrolling_existing, anchor="nw")

root.mainloop()
