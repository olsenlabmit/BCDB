import pandas as pd
import sqlite3

def create_table(uncertainty_columns):
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
            T double,
            notes text,
            BigSMILES text,
            Mn double,
            Mw double,
            D double,
            N double
        )
        """
        create += alter
        cursor.execute(create)
        columns = ["name","Mn","Mw","D","N","f","ftot","w","rho"]
        data_type = ["text","double","double","double","double","double","double","double","double"]
        for i in range(n[t]):
            for j in range(len(columns)):
                alter = "ALTER TABLE " + table[t] + " ADD " + columns[j] + str(i+1) + " " + data_type[j]
                cursor.execute(alter)
        for i in range(len(uncertainty_columns[t])):
            if "desc" in uncertainty_columns[t][i]:
                alter = "ALTER TABLE " + table[t] + " ADD " + uncertainty_columns[t][i] + " text"
            else:
                alter = "ALTER TABLE " + table[t] + " ADD " + uncertainty_columns[t][i] + " double"
            cursor.execute(alter)
        connection.commit()

        bcdb = pd.read_excel('../BCPs.xlsx',table[t],skiprows=1)
        columns = list(bcdb.columns) 
        bcdb = bcdb.values.tolist()
        for i in range(len(bcdb)):
            for j in range(len(bcdb[i])):
                if type(bcdb[i][j]) == str:
                    bcdb[i][j] = bcdb[i][j].strip()
        from pandas import DataFrame
        bcdb = DataFrame(bcdb, columns = columns)
        bcdb.to_sql(name=table[t],con=connection,if_exists='append',index=False)
            
    connection.close()

def checks():
    bcdb = pd.read_excel('../BCPs.xlsx',"diblock",skiprows=1)
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
    
    fig, ax = plt.subplots()
    C = np.array([[228, 26, 28],[55,126,184],[152,78,163],[77,175,74],[0,0,0],[160, 160, 160]])
    marker = ["s", "^", "d", "o", "*", "o"]
    label = ["Lam", "Cyl", "Gyr", "Sph", "Dis", "Other"]
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
##    ax.set_xscale('log')
##    ax.set_xlim(10**3, 10**6)
##    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
##    ax.tick_params(axis = 'both', which='major', length = 10, width = 1, direction='in', top = True, right=True)
##    ax.tick_params(axis = 'both', which='minor', length = 5, width = 1, direction='in', top = True, right=True)
##    plt.xlabel("$M_n$ (g/mol)", fontsize = 20)
##    plt.ylabel("$T$ ($^{\circ}$C)", fontsize = 20)
##    plt.xticks(fontsize = 12)
##    plt.yticks(fontsize = 12)
##    plt.tight_layout()
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

def modifiedBS():  
    bcdb = pd.read_excel('../BCPs.xlsx',"diblock",skiprows=1)
    x = []
    for i in bcdb["BigSMILES"]:
        if i[0:4] in ["{[<]","{[$]","{[>]"]:
            i = "{[]" + i[4:]
        if i[-4:] in ["[<]}","[$]}","[>]}"]:
            i = i[0:-4] + "[]}"
        print(i)
        x.append(i)
    from pandas import DataFrame
    x = DataFrame(x, columns=["BigSMILES"])
    x.to_excel("modifiedBS.xlsx","Sheet1")

def orcid_doi():
     articles = pd.read_excel('../BCPs.xlsx',"Articles")
     articles_o = articles["ORCID"]
     articles_d = articles["DOI"]
     articles_docID = articles["docID"]
     bcdb = pd.read_excel('../BCPs.xlsx',"Rheology")
     bcdb_docID = bcdb["docID"]
     d = []
     print(articles_docID)
     for i in range(len(bcdb_docID)):
         for j in range(len(articles_docID)):
             if bcdb_docID[i].srip() == articles_docID[j].strip():
                 d.append([articles_o[j], articles_d[j]])
     from pandas import DataFrame
     d = DataFrame(d, columns=["ORCID","DOI"])
     d.to_excel("modifiedORCID.xlsx","Sheet1")

def include_uncertainty():
    system_overall = ["T","Mn","Mw","D","N"]
    block = ["Mn","Mw","D","N","f","ftot","w","rho"]
    n = [2,3,4,6,11]
    from pandas import DataFrame
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
            
#create_table()
#orcid_doi()  
#checks()
l = include_uncertainty()
create_table(l)

#check table columns
##connection = sqlite3.connect('BCDB.db')
##cursor = connection.cursor()
##list(cursor.execute("select * from diblock"))
##columns = list(map(lambda x: x[0], cursor.description))
##print(columns)
     
