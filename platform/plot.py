import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.colors as mcolors
import matplotlib.markers
from matplotlib.figure import Figure
from matplotlib import rcParams, cycler, rc
import numpy as np
import sqlite3

def plot_style():

    dark_gray = "k" 
    light_gray = ".8"

    style_dict = { 
                    "figure.facecolor": "white",
                    "text.color": dark_gray,
                    "axes.labelcolor": dark_gray,
                    "legend.edgecolor": "black",
                    "legend.frameon": True,
                    "legend.numpoints": 1,
                    "legend.scatterpoints": 1,
                    "legend.fancybox": False,
                    "legend.labelspacing": 0.40,
                    "legend.handlelength": 1.25,
                    "legend.handletextpad": 0.40,
                    "legend.borderaxespad": 0.75,
                    "legend.borderpad": 0.40,
                    "xtick.direction": "out",
                    "ytick.direction": "out",
                    "xtick.color": dark_gray,
                    "ytick.color": dark_gray,
                    "axes.axisbelow": True,
                    "image.cmap": "Greys",
                    "font.family": ["sans-serif"],
                    "font.size": 12, 
                    "font.sans-serif": [
                        "Arial",
                    ],  
                    "grid.linestyle": "-",
                    "axes.grid": False,
                    "lines.solid_capstyle": "round",
                    "axes.facecolor": "white",
                    "axes.edgecolor": "black",
                    "axes.linewidth": 1.0,
                    "grid.color": light_gray,
                    "xtick.major.size": 0,
                    "ytick.major.size": 0,
                    "xtick.minor.size": 0,
                    "ytick.minor.size": 0,
                    "text.usetex": True,
                    "text.latex.preamble": r'\usepackage{sansmath}'
                }
    return style_dict

def visualize(matches, default, choose_plots):
    
    connection = sqlite3.connect('BCDB.db')
    cursor = connection.cursor()
    subset = list(cursor.execute("""select * from diblock"""))

    T = [[]]
    Mn = [[]]
    fA = [[]]
    color = [default]
    size = [5]
    shape = ["o"]
    for s in subset:
        T[-1].append(float(s[6]))
        Mn[-1].append(float(s[10]))
        fA[-1].append(float(s[27]))

    for s in range(len(matches)):
        for t in range(len(matches[s])):
            if matches[s][t][0] == "diblock":
                hits = matches[s][t][1]
                color.append(matches[s][t][2])
                size.append(matches[s][t][3])
                shape.append(matches[s][t][4])
                T.append([])
                Mn.append([])
                fA.append([])
                for index in hits:
                    T[-1].append(float(subset[index][6]))
                    Mn[-1].append(float(subset[index][10]))
                    fA[-1].append(float(subset[index][27]))
            elif matches[s][t][0] == "benchmark":
                T.append(matches[s][t][1])
                fA.append(matches[s][t][2])
                Mn.append(matches[s][t][3])
                color.append(matches[s][t][4])
                size.append(matches[s][t][5])
                shape.append(matches[s][t][6])

    ### Comment out this line if latex is not installed
    # matplotlib.rcParams.update(plot_style())

    if choose_plots[0].get() == 1:
        fig, ax = plt.subplots()
        for i in range(0, len(T)):
            try:
                plt.scatter(fA[i], T[i], c = np.array([color[i]])/256.0, marker = shape[i], s = size[i])
            except:
                continue
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.set_xbound(lower = 0, upper = 1)
        ax.tick_params(axis = 'both', which='major', length = 10, width = .5, direction='in', top = True, right=True)
        ax.tick_params(axis = 'both', which='minor', length = 5, width = .5, direction='in', top = True, right=True)
        plt.xlabel("$f_A$", fontsize = 20)
        plt.ylabel("$T$ ($^{\circ}$C)", fontsize = 20)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.tight_layout()
    if choose_plots[1].get() == 1:
        fig, ax = plt.subplots()
        for i in range(0, len(T)):
            try:
                plt.scatter(Mn[i], T[i], c = np.array([color[i]])/256.0, marker = shape[i], s=size[i])
            except:
                continue
        ax.set_xscale('log')
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.tick_params(axis = 'both', which='major', length = 10, width = 1, direction='in', top = True, right=True)
        ax.tick_params(axis = 'both', which='minor', length = 5, width = 1, direction='in', top = True, right=True)
        plt.xlabel("$M_n$ (g/mol)", fontsize = 20)
        plt.ylabel("$T$ ($^{\circ}$C)", fontsize = 20)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.tight_layout()
    if choose_plots[2].get() == 1:
        fig, ax = plt.subplots()
        for i in range(len(T)):
            try:
                plt.scatter(fA[i], Mn[i], c = np.array([color[i]])/256.0, marker = shape[i], s=size[i])
            except:
                continue
        ax.set_yscale('log')
        ax.set_xbound(lower = 0, upper = 1)
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax.tick_params(axis = 'both', which='major', length = 10, width = 1, direction='in', top = True, right=True)
        ax.tick_params(axis = 'both', which='minor', length = 5, width = 1, direction='in', top = True, right=True)
        plt.xlabel("$f_A$", fontsize = 20)
        plt.ylabel("$M_n$ (g/mol)", fontsize = 20)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.tight_layout()        
    if choose_plots[3].get() == 1:
        fig, ax = plt.subplots()
        plt.hist(Mn[0], bins = 25, color=(0.8941176470588236, 0.10196078431372549, 0.10980392156862745), log = True)
        ax.ticklabel_format(axis = 'x', style = 'sci', scilimits=(0,0))
        ax.xaxis.offsetText.set_fontsize(18)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(axis = 'both', which='major', length = 10, width = 1, direction='in', top = True, right=True)
        ax.tick_params(axis = 'both', which='minor', length = 5, width = 1, direction='in', top = True, right=True)
        plt.xlabel("$M_n$ (g/mol)", fontsize = 20)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.tight_layout()
    if choose_plots[4].get() == 1:
        fig, ax = plt.subplots()
        plt.hist(T[0], bins = 100, color=(0.8941176470588236, 0.10196078431372549, 0.10980392156862745))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.tick_params(axis = 'both', which='major', length = 10, width = 1, direction='in', top = True, right=True)
        ax.tick_params(axis = 'both', which='minor', length = 5, width = 1, direction='in', top = True, right=True)
        plt.xlabel("$T$ ($^{\circ}$C)", fontsize = 20)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.tight_layout()
    if choose_plots[5].get() == 1:
        fig, ax = plt.subplots()
        plt.hist(fA[0], bins =  50, color=(0.8941176470588236, 0.10196078431372549, 0.10980392156862745))
        ax.set_xbound(lower = -.05, upper = 1.05)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.tick_params(axis = 'both', which='major', length = 10, width = 1, direction='in', top = True, right=True)
        ax.tick_params(axis = 'both', which='minor', length = 5, width = 1, direction='in', top = True, right=True)
        plt.xlabel("$f_A$", fontsize = 20)
        plt.xticks(fontsize = 12)
        plt.yticks(fontsize = 12)
        plt.tight_layout()
    plt.show()
    connection.close()
