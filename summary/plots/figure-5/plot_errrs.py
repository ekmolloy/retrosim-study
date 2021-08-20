import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
import sys


plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{helvet} \usepackage{sfmath}']  


fig = plt.figure(figsize=(7.5, 3))
gs = gridspec.GridSpec(1, 2)
ax0 = plt.subplot(gs[0,0])
ax1 = plt.subplot(gs[0,1])

axs = [ax0, ax1]

df = pandas.read_csv("../../csvs/mean-species-tree-error.csv")
df = df[df["NRET"] != 10]
df = df[df["NRET"] != 50]

for i, ax in enumerate(axs):
    if i == 0:
        modl = "26taxa"
    else:
        modl = "Palaeognathae"

    xdf = df[df["MODL"] == modl]

    pars = xdf[xdf["MTHD"] == "parsimony-strict"]
    cs = xdf[xdf["MTHD"] == "camin-sokal-strict"]
    dollo = xdf[xdf["MTHD"] == "dollo-strict"]
    astrid = xdf[xdf["MTHD"] == "astrid-bp"]
    astral = xdf[xdf["MTHD"] == "astral-bp"]
    mdc = xdf[xdf["MTHD"] == "mdc-bp-ur"]

    # black, gray

    xs = [100 - 10, 500 - 50, 1000 - 100, 5000 - 500, 10000 - 1000, 50000 - 5000, 100000 - 10000]
    ax.errorbar(xs, pars.AVG_FNR.values, yerr=pars.SE_FNR.values,
                fmt='.-', markersize=8, capthick=2, label="Unordered")
                #, color="black", )

    xs = [100 - 6, 500 - 30, 1000 - 60, 5000 - 300, 10000 - 600, 50000 - 3000, 100000 - 6000]
    ax.errorbar(xs, dollo.AVG_FNR.values, yerr=dollo.SE_FNR.values,
                fmt='.-', markersize=8, capthick=2, label="Dollo")
                #color="gray")
    
    xs = [100 - 2, 500 - 10, 1000 - 20, 5000 - 100, 10000 - 200, 50000 - 1000, 100000 - 2000]
    #xs = [100 - 1.5, 500 - 4.5, 1000 - 15, 5000 - 45, 10000 - 150, 50000 - 450, 100000 - 1500]
    ax.errorbar(xs, cs.AVG_FNR.values, yerr=cs.SE_FNR.values,
                fmt='.-', markersize=8, capthick=2, label="Camin-Sokal")
                #color="red")

    xs = [100 + 2, 500 + 10, 1000 + 20, 5000 + 100, 10000 + 200, 50000 + 1000, 100000 + 2000]
    #xs = [100 + 1.5, 500 + 4.5, 1000 + 15, 5000 + 45, 10000 + 150, 50000 + 450, 100000 + 1500]
    ax.errorbar(xs, mdc.AVG_FNR.values, yerr=mdc.SE_FNR.values,
                fmt='.-', markersize=8, capthick=2, label=r"MDC\_BP")
                #color="magenta")

    xs = [100 + 6, 500 + 30, 1000 + 60, 5000 + 300, 10000 + 600, 50000 + 3000, 100000 + 6000]
    ax.errorbar(xs, astrid.AVG_FNR.values, yerr=astrid.SE_FNR.values,
                fmt='.-', markersize=8, capthick=2, label=r"ASTRID\_BP")
                #color="green")

    xs = [100 + 10, 500 + 50, 1000 + 100, 5000 + 500, 10000 + 1000, 50000 + 5000, 100000 + 10000]
    ax.errorbar(xs, astral.AVG_FNR.values, yerr=astral.SE_FNR.values,
                fmt='.-', markersize=8, capthick=2, label=r"ASTRAL\_BP")
                #color="blue")

    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=11)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xscale('log')
    ax.set_xlabel(r'Number of RIs', fontsize=14)
    ax.set_xticks(astral.NRET.values)
    ax.set_yticks([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
    ax.set_ylim(0, 0.35)

    if i == 0:
        #ax.set_ylabel(r'FN rate', fontsize=14)
        ax.set_ylabel(r'Species Tree Error', fontsize=14)
        ax.text(85, 0.35, r"A)", fontsize=12)
        ax.set_title(r"26-taxa", fontsize=16)
        ax.legend(frameon=False, ncol=1, fontsize=10, loc="upper right")
    else:
        ax.text(85, 0.35, r"B)", fontsize=12)
        ax.set_title(r"Palaeognathae-GT-AZ", fontsize=16)


gs.tight_layout(fig, rect=[0, 0.0, 1, 1])
plt.savefig("errrs.pdf", format='pdf', dpi=300)
