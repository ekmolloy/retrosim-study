import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
import sys


plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{helvet} \usepackage{sfmath}']  


fig = plt.figure(figsize=(7.5, 5))
gs = gridspec.GridSpec(2, 2)
ax0 = plt.subplot(gs[0,0])
ax1 = plt.subplot(gs[0,1])
ax2 = plt.subplot(gs[1,0])
ax3 = plt.subplot(gs[1,1])

axs = [ax0, ax1, ax2, ax3]

df = pandas.read_csv("../../csvs/mean-species-tree-error-all.csv")
df = df[df["MTHD"] == "astral-bp"]
df = df[df["NRET"] != 10]
df = df[df["NRET"] != 50]

for i, ax in enumerate(axs):
    if i == 0:
        modl = "26taxa"
    else:
        modl = "Palaeognathae"

    xdf = df[df["MODL"] == modl]

    t00 = xdf[xdf["THRC"] == 0.0]
    t40 = xdf[xdf["THRC"] == 0.4]
    t50 = xdf[xdf["THRC"] == 0.5]
    t60 = xdf[xdf["THRC"] == 0.6]
    t70 = xdf[xdf["THRC"] == 0.7]
    t80 = xdf[xdf["THRC"] == 0.8]
    t90 = xdf[xdf["THRC"] == 0.9]

    if (i == 0) or (i == 1):
        ax.errorbar(t00.NRET.values, t00.AVG_FNR.values, yerr=t00.SE_FNR.values,
                    fmt='.-', markersize=8, capthick=2, label="0.0")
                    #color="black")
        ax.errorbar(t50.NRET.values, t50.AVG_FNR.values, yerr=t50.SE_FNR.values,
                    fmt='.-', markersize=8, capthick=2, label="0.5")
                    #color="gray")
        ax.errorbar(t60.NRET.values, t60.AVG_FNR.values, yerr=t60.SE_FNR.values,
                    fmt='.-', markersize=8, capthick=2, label="0.6")
                    #color="red")
        ax.errorbar(t70.NRET.values, t70.AVG_FNR.values, yerr=t70.SE_FNR.values,
                    fmt='.-', markersize=8, capthick=2, label=r"0.7")
                    #color="magenta")
        ax.errorbar(t80.NRET.values, t80.AVG_FNR.values, yerr=t80.SE_FNR.values,
                    fmt='.-', markersize=8, capthick=2, label=r"0.8")
                    #color="green")
        ax.errorbar(t90.NRET.values, t90.AVG_FNR.values, yerr=t90.SE_FNR.values,
                    fmt='.-', markersize=8, capthick=2, label=r"0.9")
                    #color="blue")
    else:
        ax.errorbar(t00.NRET.values, t00.AVG_FPR.values, yerr=t00.SE_FPR.values,
                    fmt='.-', markersize=8, capthick=2, label="0.0")
                    #color="black")
        ax.errorbar(t50.NRET.values, t50.AVG_FPR.values, yerr=t50.SE_FPR.values,
                    fmt='.-', markersize=8, capthick=2, label="0.5")
                    #color="gray")
        ax.errorbar(t60.NRET.values, t60.AVG_FPR.values, yerr=t60.SE_FPR.values,
                    fmt='.-', markersize=8, capthick=2, label="0.6")
                    #color="red")
        ax.errorbar(t70.NRET.values, t70.AVG_FPR.values, yerr=t70.SE_FPR.values,
                    fmt='.-', markersize=8, capthick=2, label=r"0.7")
                    #color="magenta")
        ax.errorbar(t80.NRET.values, t80.AVG_FPR.values, yerr=t80.SE_FPR.values,
                    fmt='.-', markersize=8, capthick=2, label=r"0.8")
                    #color="green")
        ax.errorbar(t90.NRET.values, t90.AVG_FPR.values, yerr=t90.SE_FPR.values,
                    fmt='.-', markersize=8, capthick=2, label=r"0.9")
                    #color="blue")

    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=11)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xscale('log')
    ax.set_xticks(t00.NRET.values)

    if (i == 0) or (i == 1):
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
        ax.set_ylim(0, 0.6)
    else:
        ax.set_yticks([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.30])
        ax.set_ylim(0, 0.30)

    if (i == 0):
        ax.set_ylabel(r'FN rate', fontsize=14)
        ax.set_title(r"26-taxon-Sim", fontsize=16)
        ax.legend(frameon=False, ncol=1, fontsize=10, loc="upper right")
    elif (i == 1):
        ax.set_title(r"Palaeognathae-GT-AZ", fontsize=16)
    elif (i == 2):
        ax.set_ylabel(r'FP rate', fontsize=14)
        ax.set_xlabel(r'Number of RIs', fontsize=14)
    else:
        ax.set_xlabel(r'Number of RIs', fontsize=14)

    ax0.text(85, 0.6, r"A)", fontsize=12)
    ax1.text(85, 0.6, r"C)", fontsize=12)
    ax2.text(85, 0.3, r"B)", fontsize=12)
    ax3.text(85, 0.3, r"D)", fontsize=12)

gs.tight_layout(fig, rect=[0, 0.0, 1, 1])
plt.savefig("astral.pdf", format='pdf', dpi=300)
