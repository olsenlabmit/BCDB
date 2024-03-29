import pandas as pd
from pandas import DataFrame
import sqlite3
import rdkit
from rdkit import Chem

from BigSMILES_BigSmilesObj import BigSMILES

def create_table():
    def include_uncertainty():
        system_overall = ["T","meas","anneal","Talt","Mn","Mw","D","N"]
        block = ["Mn","Mw","D","N","f","ftot","w","rho"]
        n = [2,3,4,6,11]
        tot = []
        for table in n:
            all_columns = []
            error1 = ["std","se"]
            for i in range(len(error1)):
                for j in range(len(system_overall)):
                    all_columns.append(system_overall[j] + error1[i])
                for j in range(table):
                    for k in range(len(block)):
                        all_columns.append(block[k] + str(j + 1) + error1[i])
            error2 = ["unc","desc"]
            for j in range(len(system_overall)):
                all_columns.append(system_overall[j] + error2[0])
                all_columns.append(system_overall[j] + error2[1])
            for j in range(table):
                for k in range(len(block)):
                    all_columns.append(block[k] + str(j + 1) + error2[0])
                    all_columns.append(block[k] + str(j + 1) + error2[1])
            tot.append(all_columns)
        return tot
    uncertainty = include_uncertainty()

    table = ["diblock","triblock","tetrablock","hexablock","undecablock"]
    n = [2,3,4,6,11]
    connection = sqlite3.connect('BCDB.db')
    cursor = connection.cursor()
    for t in range(len(table)):
        create = "CREATE TABLE "  + table[t]
        alter = """
        (
            ind int,
            ORCID text,
            DOI text,
            phase1 text,
            phase2 text,
            phase_method text,
            T double,
            T_meas double,
            T_anneal double,
            T_alt double,
            T_describe text,
            notes text,
            BigSMILES text,
            Mn double,
            Mn_method text,
            Mw double,
            Mw_method text,
            D double,
            D_method text,
            N double,
            N_method text
        )
        """
        create += alter
        cursor.execute(create)
        columns = ["name","Mn","Mw","D","N","f","f_tot","w","rho"]
        data_type = ["text","double","double","double","double","double","double","double","double"]
        for i in range(n[t]):
            alter = "ALTER TABLE " + table[t] + " ADD " + columns[0] + str(i+1) + " " + data_type[0]
            cursor.execute(alter)
            for j in range(1, len(columns)):
                alter = "ALTER TABLE " + table[t] + " ADD " + columns[j] + str(i+1) + " " + data_type[j]
                cursor.execute(alter)
                alter = "ALTER TABLE " + table[t] + " ADD " + columns[j] + str(i+1) + "_method " + data_type[j]
                cursor.execute(alter)
        for i in range(len(uncertainty[t])):
            if "desc" in uncertainty[t][i]:
                alter = "ALTER TABLE " + table[t] + " ADD " + uncertainty[t][i] + " text"
            else:
                alter = "ALTER TABLE " + table[t] + " ADD " + uncertainty[t][i] + " double"
            cursor.execute(alter)
        connection.commit()

        bcdb = pd.read_excel('../data/BCPs.xlsx',table[t],skiprows=1)
        columns = list(bcdb.columns) 
        bcdb = bcdb.values.tolist()
        for i in range(len(bcdb)):
            for j in range(len(bcdb[i])):
                if type(bcdb[i][j]) == str:
                    bcdb[i][j] = bcdb[i][j].strip()
        bcdb = DataFrame(bcdb, columns = columns)
        bcdb.to_sql(name=table[t],con=connection,if_exists='append',index=False)
            
    connection.close()

def phase_check():
    bcdb = pd.read_excel('../data/BCPs.xlsx',"diblock",skiprows=1)
    bcdb = bcdb.fillna("")

    # check Figure 3
    phases = [0,0,0,0,0]
    p1 = bcdb['phase1']
    p2 = bcdb['phase2']
    T = bcdb['T']
    Mn = bcdb['Mn']
    fA = bcdb['f1']
    Tplot = [[],[],[],[],[],[]]
    Mnplot = [[],[],[],[],[],[]]
    fAplot = [[],[],[],[],[],[]]
    for i in range(len(p1)):
        if p1[i] == "lamellar" and p2[i] == "":
            phases[0] += 1
            Tplot[0].append(T[i])
            Mnplot[0].append(Mn[i])
            fAplot[0].append(fA[i])
        elif p1[i] == "cylinder" and p2[i] == "":
            phases[1] += 1
            Tplot[1].append(T[i])
            Mnplot[1].append(Mn[i])
            fAplot[1].append(fA[i])
        elif p1[i] == "gyroid" and p2[i] == "":
            phases[2] += 1
            Tplot[2].append(T[i])
            Mnplot[2].append(Mn[i])
            fAplot[2].append(fA[i])
        elif p1[i] == "sphere" and p2[i] == "":
            phases[3] += 1
            Tplot[3].append(T[i])
            Mnplot[3].append(Mn[i])
            fAplot[3].append(fA[i])
        elif p1[i] == "disordered" and p2[i] == "":
            phases[4] += 1
            Tplot[4].append(T[i])
            Mnplot[4].append(Mn[i])
            fAplot[4].append(fA[i])
        else:
            Tplot[5].append(T[i])
            Mnplot[5].append(Mn[i])
            fAplot[5].append(fA[i])
    print(phases)

    #Generate Figure 3
    import matplotlib.pyplot as plt
    from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
    import matplotlib.colors as mcolors
    import matplotlib.markers
    from matplotlib.figure import Figure
    from matplotlib import rcParams, cycler, rc
    import numpy as np
    from plot import plot_style

    matplotlib.rcParams.update(plot_style())

    C = np.array([[228, 26, 28],[55,126,184],[152,78,163],[77,175,74],[0,0,0],[160, 160, 160]])
    marker = ["s", "^", "d", "o", "*", "o"]
    label = ["Lam", "Cyl", "Gyr", "Sph", "Dis", "Other"]

    fig, ax = plt.subplots()
    for i in range(len(label)-1, -1, -1):
        plt.scatter(fAplot[i], Tplot[i], marker = marker[i], s = 10, c = C[i]/256, label = label[i])
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels))
    ax.set_xlim(0, 1.0)
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(axis = 'both', which='major', length = 10, width = 1, direction='in', top = True, right=True)
    ax.tick_params(axis = 'both', which='minor', length = 5, width = 1, direction='in', top = True, right=True)
    plt.xlabel("$f_A$", fontsize = 20)
    plt.ylabel("$T$ ($^{\circ}$C)", fontsize = 20)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.tight_layout()

    fig, ax = plt.subplots()
    for i in range(len(label)-1, -1, -1):
        plt.scatter(Mnplot[i], Tplot[i], marker = marker[i], s = 10, c = C[i]/256, label = label[i])
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(reversed(handles), reversed(labels))
    ax.set_xscale('log')
    ax.set_xlim(10**3, 10**6)
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(axis = 'both', which='major', length = 10, width = 1, direction='in', top = True, right=True)
    ax.tick_params(axis = 'both', which='minor', length = 5, width = 1, direction='in', top = True, right=True)
    plt.xlabel("$M_n$ (g/mol)", fontsize = 20)
    plt.ylabel("$T$ ($^{\circ}$C)", fontsize = 20)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.tight_layout()
    plt.show()
    
    # check Figure 6
    fA = bcdb['f1']
    fB = bcdb['f2']
    T = bcdb['T']
    nameA = bcdb['name1']
    nameB = bcdb['name2']
    Mn = bcdb['Mn']
    queries = [0,0,0,0,0]
    unique_diblock = set()
    for i in range(len(fA)):
        s = [p1[i], p2[i]]
        if "lamellar" in s:
            queries[0] += 1
        if "lamellar" in s and "disordered" in s:
            queries[1] += 1
        if "lamellar" in s and "cylinder" in s:
            queries[2] += 1
        if T[i] > 100 and ("lamellar" in s or Mn[i] > 30000):
            queries[3] += 1
        unique_diblock.add(sorted([nameA[i],nameB[i]]))
    print(len(unique_diblock))
    
def database_summary():
    table = ["diblock","triblock","tetrablock","hexablock","undecablock"]

    tot = 0
    orcid = set()
    blocks = set()
    b = [2, 3, 4, 6, 11]
    for i in range(len(table)):
        bcdb = pd.read_excel('../data/BCPs.xlsx', table[i], skiprows=1)
        bcdb = bcdb.fillna("")
        
        tot += len(bcdb["DOI"])

        for o in bcdb["DOI"]:
            orcid.add(o)
        
        for j in range(1, b[i] + 1):
            n = "name" + str(j)
            for o in bcdb[n]:
                blocks.add(o) 
        
        if b[i] == 2:
            # unique diblocks
            block_dict = {}
            diblock_dict = {}
            for j in range(len(bcdb["name1"])):
                d = tuple(sorted([bcdb["name1"][j], bcdb["name2"][j]]))
                if d not in diblock_dict:
                    diblock_dict[d] = 0
                diblock_dict[d] += 1 / len(bcdb["name1"])
                if bcdb["name1"][j] not in block_dict:
                    block_dict[bcdb["name1"][j]] = 0
                block_dict[bcdb["name1"][j]] += 1/ (2 * len(bcdb["name1"]))
                if bcdb["name2"][j] not in block_dict:
                    block_dict[bcdb["name2"][j]] = 0
                block_dict[bcdb["name2"][j]] += 1/ (2 * len(bcdb["name1"]))
            block_dict = sorted(block_dict.items(), key=lambda x:x[1], reverse = True)
            diblock_dict = sorted(diblock_dict.items(), key=lambda x:x[1], reverse = True)
            print(block_dict)
            print(diblock_dict)
            print("# unique blocks in diblock page: ", len(block_dict))
            print("# unique diblocks in diblock page: ", len(diblock_dict))
            
            # single phase
            single = 0
            phase_dict = {}
            for j in range(len(bcdb["phase1"])):
                if bcdb["phase2"][j] == "":
                    single += 1
                    if bcdb["phase1"][j] not in phase_dict:
                        phase_dict[bcdb["phase1"][j]] = 0
                    phase_dict[bcdb["phase1"][j]] += 1/ (len(bcdb["phase1"]))
            phase_dict = sorted(phase_dict.items(), key=lambda x:x[1], reverse = True)
            print(phase_dict)
            print("% pure phase in diblock page: ", single/len(bcdb["phase1"]))
                
    print("total data points: ", tot)
    print("total papers:  ", len(orcid))
    print("total unique blocks: ", len(blocks))

def main_Figure5():
    bcdb = pd.read_excel('../data/BCPs.xlsx', "diblock", skiprows=1)
    bcdb = bcdb.fillna("")
    query = [0, 0, 0]
    smarts = Chem.MolFromSmarts("[R]")
    before = ""
    before_R = False
    for i in range(len(bcdb["phase1"])):
        if bcdb["phase1"][i] == "lamellar" and bcdb["phase2"][i] == "":
            a = bcdb["name1"][i]
            b = bcdb["name2"][i]
            if a == "PS" and b == "PI" or a == "PI" and b == "PS":
                query[0] += 1
        try:
            if bcdb["BigSMILES"][i] == before:
                if before_R:
                    query[1] += 1
            else:
                before = bcdb["BigSMILES"][i]
                Polymer = BigSMILES(bcdb["BigSMILES"][i])
                both_R = True
                for obj in Polymer:
                    found = False
                    for rep in obj:
                        smiles = Chem.MolFromSmiles(rep.writeStandard(noBondDesc = True))
                        if smiles.HasSubstructMatch(smarts):
                            found = True
                            break
                    if not found:
                        both_R = False
                        break
                before_R = both_R
                if both_R:
                    query[1] += 1
        except:
            x = 1
        if float(bcdb["f1"][i]) < 0.25 or float(bcdb["f2"][i]) < 0.25:
            query[2] += 1
        if i % 100 == 0:
            print(i)
    print(query)

def SI_Table4():
    bcdb = pd.read_excel('../data/BCPs.xlsx', "diblock", skiprows=1)
    bcdb = bcdb.fillna("")
    query = [0, 0, 0, 0, 0]
    for i in range(len(bcdb["phase1"])):
        if bcdb["phase1"][i] == "lamellar" and bcdb["phase2"][i] == "":
            query[0] += 1
        if bcdb["phase1"][i] == "lamellar" or bcdb["phase2"][i] == "lamellar":
            query[1] += 1
        if bcdb["phase1"][i] == "lamellar" and bcdb["phase2"][i] == "disordered" or \
            bcdb["phase1"][i] == "disordered" and bcdb["phase2"][i] == "lamellar":
            query[2] += 1
        if bcdb["phase1"][i] == "lamellar" and bcdb["phase2"][i] == "cylinder" or \
            bcdb["phase1"][i] == "cylinder" and bcdb["phase2"][i] == "lamellar":
            query[3] += 1
        if (bcdb["phase1"][i] == "lamellar" or bcdb["phase2"][i] == "lamellar" or int(bcdb["Mn"][i]) > 30000) and \
            int(bcdb["T"][i]) > 100:
            query[4] += 1
    print(query)

# create_table()
# database_summary()
# SI_Table4()
main_Figure5()

# check table columns
# connection = sqlite3.connect('BCDB.db')
# cursor = connection.cursor()
# list(cursor.execute("select * from diblock"))
# columns = list(map(lambda x: x[0], cursor.description))
# print(columns)