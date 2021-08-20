import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
import sys


plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{helvet} \usepackage{sfmath}']  


fig = plt.figure(figsize=(7.5, 6))
gs = gridspec.GridSpec(3,2)
ax0 = plt.subplot(gs[0,0])
ax1 = plt.subplot(gs[0,1])
ax2 = plt.subplot(gs[1,0])
ax3 = plt.subplot(gs[1,1])
ax4 = plt.subplot(gs[2,0])
ax5 = plt.subplot(gs[2,1])

axs = [ax0, ax1, ax2, ax3, ax4, ax5]


df = pandas.read_csv("../../csvs/26taxa-mean-simulated-branch-lengths.csv")
inds = numpy.argsort(df.BRLN_TRUE.values)
xs = df.BRLN_TRUE.values[inds]
keep = numpy.where(xs <= 5.0)

avg_brln_map_gt = df.AVG_BRLN_MAP_GT.values[inds][keep]
std_brln_map_gt = df.STD_BRLN_MAP_GT.values[inds][keep]
avg_aerr_map_gt = df.AVG_AERR_MAP_GT.values[inds][keep]
std_aerr_map_gt = df.STD_AERR_MAP_GT.values[inds][keep]
avg_perr_map_gt = df.AVG_PERR_MAP_GT.values[inds][keep]
std_perr_map_gt = df.STD_PERR_MAP_GT.values[inds][keep]

avg_brln_mle_gt = df.AVG_BRLN_MLE_GT.values[inds][keep]
std_brln_mle_gt = df.STD_BRLN_MLE_GT.values[inds][keep]
avg_aerr_mle_gt = df.AVG_AERR_MLE_GT.values[inds][keep]
std_aerr_mle_gt = df.STD_AERR_MLE_GT.values[inds][keep]
avg_perr_mle_gt = df.AVG_PERR_MLE_GT.values[inds][keep]
std_perr_mle_gt = df.STD_PERR_MLE_GT.values[inds][keep]

avg_brln_mle_ri = df.AVG_BRLN_MLE_RI.values[inds][keep]
std_brln_mle_ri = df.STD_BRLN_MLE_RI.values[inds][keep]
avg_aerr_mle_ri = df.AVG_AERR_MLE_RI.values[inds][keep]
std_aerr_mle_ri = df.STD_AERR_MLE_RI.values[inds][keep]
avg_perr_mle_ri = df.AVG_PERR_MLE_RI.values[inds][keep]
std_perr_mle_ri = df.STD_PERR_MLE_RI.values[inds][keep]

#ax0.errorbar(xs[keep], avg_brln_map_gt, yerr=std_brln_map_gt,
#	fmt='.', markersize=8, capthick=2, color="blue", label=r"MAP-GT")
ax0.errorbar(xs[keep], avg_brln_mle_gt, yerr=std_brln_mle_gt,
	fmt='.', markersize=8, capthick=2, color="black", label=r"MLE-GT")
ax0.errorbar(xs[keep], avg_brln_mle_ri, yerr=std_brln_mle_ri,
	fmt='.', markersize=8, color="gray", label=r"MLE-RI")
ax0.plot(xs[keep], xs[keep], '-', color="black", label=r"True")
ax0.set_yticks([0, 1, 2, 3, 4, 5, 6])
ax0.set_ylim(0, 6)
ax0.legend(frameon=False, ncol=1, fontsize=10, loc="upper left")

#ax2.errorbar(xs[keep], avg_aerr_map_gt, yerr=std_aerr_map_gt,
#	fmt='.-', markersize=8, capthick=2, color="black", label=r"MAP-GT")
ax2.errorbar(xs[keep], avg_aerr_mle_gt, yerr=std_aerr_mle_gt,
	fmt='.-', markersize=8, capthick=2, color="black", label=r"MLE-GT")
ax2.errorbar(xs[keep], avg_aerr_mle_ri, yerr=std_aerr_mle_ri,
	fmt='.-', markersize=8, color="gray", label=r"MLE-RI")
ax2.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
ax2.set_ylim(0, 2.0)
ax2.legend(frameon=False, ncol=1, fontsize=10, loc="upper left")

#ax4.errorbar(xs[keep], avg_perr_map_gt, yerr=std_perr_map_gt,
#	fmt='.-', markersize=8, capthick=2, color="black", label=r"MAP-GT")
ax4.errorbar(xs[keep], avg_perr_mle_gt, yerr=std_perr_mle_gt,
	fmt='.-', markersize=8, capthick=2, color="black", label=r"MLE-GT")
ax4.errorbar(xs[keep], avg_perr_mle_ri, yerr=std_perr_mle_ri,
	fmt='.-', markersize=8, color="gray", label=r"MLE-RI")
ax4.set_yticks([0, 10, 20, 30, 40, 50])
ax4.set_ylim(0, 50)
#ax4.legend(frameon=False, ncol=1, fontsize=10, loc="center right")

df = pandas.read_csv("../../csvs/Palaeognathae-mean-simulated-branch-lengths.csv")
inds = numpy.argsort(df.BRLN_TRUE.values)
xs = df.BRLN_TRUE.values[inds]
keep = numpy.where(xs <= 5.0)

avg_brln_map_gt = df.AVG_BRLN_MAP_GT.values[inds][keep]
std_brln_map_gt = df.STD_BRLN_MAP_GT.values[inds][keep]
avg_aerr_map_gt = df.AVG_AERR_MAP_GT.values[inds][keep]
std_aerr_map_gt = df.STD_AERR_MAP_GT.values[inds][keep]
avg_perr_map_gt = df.AVG_PERR_MAP_GT.values[inds][keep]
std_perr_map_gt = df.STD_PERR_MAP_GT.values[inds][keep]

avg_brln_mle_gt = df.AVG_BRLN_MLE_GT.values[inds][keep]
std_brln_mle_gt = df.STD_BRLN_MLE_GT.values[inds][keep]
avg_aerr_mle_gt = df.AVG_AERR_MLE_GT.values[inds][keep]
std_aerr_mle_gt = df.STD_AERR_MLE_GT.values[inds][keep]
avg_perr_mle_gt = df.AVG_PERR_MLE_GT.values[inds][keep]
std_perr_mle_gt = df.STD_PERR_MLE_GT.values[inds][keep]

avg_brln_mle_ri = df.AVG_BRLN_MLE_RI.values[inds][keep]
std_brln_mle_ri = df.STD_BRLN_MLE_RI.values[inds][keep]
avg_aerr_mle_ri = df.AVG_AERR_MLE_RI.values[inds][keep]
std_aerr_mle_ri = df.STD_AERR_MLE_RI.values[inds][keep]
avg_perr_mle_ri = df.AVG_PERR_MLE_RI.values[inds][keep]
std_perr_mle_ri = df.STD_PERR_MLE_RI.values[inds][keep]

#ax1.errorbar(xs[keep], avg_brln_map_gt, yerr=std_brln_map_gt,
#	fmt='.', markersize=8, capthick=2, color="blue", label=r"MAP-GT")
ax1.errorbar(xs[keep], avg_brln_mle_gt, yerr=std_brln_mle_gt,
	fmt='.', markersize=8, capthick=2, color="black", label=r"MLE-GT")
ax1.errorbar(xs[keep], avg_brln_mle_ri, yerr=std_brln_mle_ri,
	fmt='.', markersize=8, color="gray", label=r"MLE-RI")
ax1.plot(xs[keep], xs[keep], '-', color="black", label=r"True")
#ax1.set_ylabel(r'Branch Length', fontsize=14)
ax1.set_yticks([0, 1, 2, 3, 4, 5, 6])
ax1.set_ylim(0, 6)
#ax1.legend(frameon=False, ncol=1, fontsize=10, loc="upper left")

#ax3.errorbar(xs[keep], avg_aerr_map_gt, yerr=std_aerr_map_gt,
#	fmt='.-', markersize=8, capthick=2, color="black", label=r"MAP-GT")
ax3.errorbar(xs[keep], avg_aerr_mle_gt, yerr=std_aerr_mle_gt,
	fmt='.-', markersize=8, capthick=2, color="black", label=r"MLE-GT")
ax3.errorbar(xs[keep], avg_aerr_mle_ri, yerr=std_aerr_mle_ri,
	fmt='.-', markersize=8, color="gray", label=r"MLE-RI")
#ax3.set_ylabel(r'Error', fontsize=14)
ax3.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
ax3.set_ylim(0, 2.0)
#ax3.legend(frameon=False, ncol=1, fontsize=10, loc="upper left")

#ax5.errorbar(xs[keep], avg_perr_map_gt, yerr=std_perr_map_gt,
#	fmt='.-', markersize=8, capthick=2, color="black", label=r"MAP-GT")
ax5.errorbar(xs[keep], avg_perr_mle_gt, yerr=std_perr_mle_gt,
	fmt='.-', markersize=8, capthick=2, color="black", label=r"MLE-GT")
ax5.errorbar(xs[keep], avg_perr_mle_ri, yerr=std_perr_mle_ri,
	fmt='.-', markersize=8, color="gray", label=r"MLE-RI")
#ax5.set_ylabel(r'Percent Error', fontsize=14)
ax5.set_yticks([0, 10, 20, 30, 40, 50])
ax5.set_ylim(0, 50)
#ax5.legend(frameon=False, ncol=1, fontsize=10, loc="center right")

for i, ax in enumerate(axs):
    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=11)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])

ax0.set_title(r'26-taxa', fontsize=16)
ax1.set_title(r'Palaeognathae-GT-AZ', fontsize=16)

ax0.text(0.0, 6, r"A)", fontsize=12)
ax2.text(0.0, 2, r"B)", fontsize=12)
ax4.text(0.0, 50, r"C)", fontsize=12)
ax1.text(0.0, 6, r"D)", fontsize=12)
ax3.text(0.0, 2, r"E)", fontsize=12)
ax5.text(0.0, 50, r"F)", fontsize=12)

ax4.set_xlabel(r'True Branch Length', fontsize=14)
ax5.set_xlabel(r'True Branch Length', fontsize=14)

ax0.set_ylabel(r'Estimated Length', fontsize=14)
ax2.set_ylabel(r'Absolute Error', fontsize=14)
ax4.set_ylabel(r'Percent Error', fontsize=14)

gs.tight_layout(fig, rect=[0, 0.0, 1, 1])
plt.savefig("brlns.pdf", format='pdf', dpi=300)
