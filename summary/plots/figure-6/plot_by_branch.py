import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
import sys


plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{helvet} \usepackage{sfmath}']  


fig = plt.figure(figsize=(7.5, 6))
gs = gridspec.GridSpec(3, 2)
ax0 = plt.subplot(gs[0,0])
ax1 = plt.subplot(gs[0,1])
ax2 = plt.subplot(gs[1,0])
ax3 = plt.subplot(gs[1,1])
ax4 = plt.subplot(gs[2,0])
ax5 = plt.subplot(gs[2,1])

axs = [ax0, ax1, ax2, ax3, ax4, ax5]

for i, ax in enumerate(axs):
    if (i == 0) or (i == 2) or (i == 4):
        df = pandas.read_csv("../../csvs/26taxa-mean-simulated-branch-lengths-varied.csv")
        # print(numpy.unique(df.BRLN_TRUE.values))    
        brlns = [0.1, 0.2, 0.25, 0.3, 0.4, 0.5,
                 0.6, 0.7, 0.8, 0.9, 1.3, 1.5,
                 1.9, 2., 2.5, 2.9, 3.6, 4.1,
                 6., 7.]
        brlns = [0.1, 0.2, 0.5, 1.3, 2.]
    else:
        df = pandas.read_csv("../../csvs/Palaeognathae-mean-simulated-branch-lengths-varied.csv")
        # print(numpy.unique(df.BRLN_TRUE.values))
        brlns = [0.019359, 0.053167, 0.309131, 0.387391,
                 1.0655, 1.15143, 1.95468, 3.01202, 3.0872,
                 4.16743]
        brlns = [0.019359, 0.053167, 0.387391, 1.0655, 1.95468]

    df = df[df["NRET"] != 10]
    df = df[df["NRET"] != 50]

    if i < 2:
        for brln in brlns:
            xdf = df[df["BRLN_TRUE"] == brln]
            ax.errorbar(xdf.NRET.values, xdf.AVG_EN.values, #yerr=xdf.STD_EN.values,
                        fmt='.-', markersize=8, capthick=2,
                        label=str("%1.2f" % brln))
            ax.set_yticks([1, 100, 1000, 10000])
            ax.set_ylim(1, 20000)
            ax.set_yscale('log')
            if i == 0:
                ax.legend(frameon=False, ncol=2, fontsize=10, loc="upper left")
            else:
                ax.legend(frameon=False, ncol=2, fontsize=10, loc="lower right")
    elif (i > 1) and (i < 4):
        for brln in brlns:
            xdf = df[df["BRLN_TRUE"] == brln]
            ax.errorbar(xdf.NRET.values, (xdf.NREP.values / 25) * 100,
                        fmt='.-', markersize=8, capthick=2,
                        label=str("%1.2f" % brln))
            #ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
            #ax.set_ylim(0, 1.05)
            ax.set_yticks([0, 20, 40, 60, 80, 100])
            ax.set_ylim(0, 105)
    elif (i > 3):
        for brln in brlns:
            xdf = df[df["BRLN_TRUE"] == brln]
            ax.errorbar(xdf.NRET.values, xdf.AVG_PERR_MLE_RI.values, 
                        fmt='.-', markersize=8, capthick=2,
                        label=str("%1.2f" % brln))
            ax.set_yticks([0, 20, 40, 60, 80, 100])
            ax.set_ylim(0, 105)

    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=11)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xscale('log')

    if (i == 0):
        ax.set_ylabel(r'$EN$', fontsize=14)
        ax.set_title(r"26-taxa", fontsize=16)
    elif (i == 1):
        ax.set_title(r"Palaeognathae-GT-AZ", fontsize=16)
    elif (i == 2):
        ax.set_ylabel(r"\% replicates", fontsize=14)
    if (i == 4):
        ax.set_ylabel(r"\% error for length", fontsize=14)

    if (i > 3):
        ax.set_xlabel(r'Number of RIs', fontsize=14)

    ax0.text(85, 20000, r"A)", fontsize=12)
    ax1.text(85, 20000, r"D)", fontsize=12)
    ax2.text(85, 105, r"B)", fontsize=12)
    ax3.text(85, 105, r"E)", fontsize=12)
    ax4.text(85, 105, r"C)", fontsize=12)
    ax5.text(85, 105, r"F)", fontsize=12)

gs.tight_layout(fig, rect=[0, 0.0, 1, 1])
plt.savefig("by_branch.pdf", format='pdf', dpi=300)
