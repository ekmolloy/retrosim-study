import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from scipy import special
import sys

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{helvet} \usepackage{sfmath}']  


def prior_ri(x, l):
    y = ((2. / 3) / (1 - x)) - 1.
    a = (1 + y)**2 / (special.lambertw(y).real + 1)
    b = (special.lambertw(y).real / y)**(2 + l)
    return l * a * b

def prior_gt(x, l):
    a = (2*l - 1)
    b = ((3. / 2) * (1 - x))**a
    return l * b


xs = numpy.array([0.001, 0.0025, 0.005, 0.0075,
                 0.01, 0.025, 0.05, 0.075,
                 0.1, 0.25, 0.5, 0.75,
                 1.0, 2.5, 3.0])


fig = plt.figure(figsize=(8, 6))
gs = gridspec.GridSpec(1,2)
ax0 = plt.subplot(gs[0,0])
ax1 = plt.subplot(gs[0,1])

axs = [ax0, ax1]

xs = numpy.arange(0.33, 0.99, 0.001)

# Lemma 2 in Sayyari and Mirarab
ax0.set_title(r"Gene Trees", fontsize=20)
ax0.plot(xs, prior_gt(xs, 1.0), '-', label=r"$\lambda = 1.0$", linewidth=2)
ax0.plot(xs, prior_gt(xs, 0.75), '-', label=r"$\lambda = 0.75$", linewidth=2)
ax0.plot(xs, prior_gt(xs, 0.5), '-', label=r"$\lambda = 0.5$", linewidth=2)
ax0.plot(xs, prior_gt(xs, 0.25), '-', label=r"$\lambda = 0.25$", linewidth=2)
ax0.plot(xs, prior_gt(xs, 0.1), '-', label=r"$\lambda = 0.1$", linewidth=2)
ax0.plot(xs, prior_gt(xs, 0.05), '-', label=r"$\lambda = 0.05$", linewidth=2)
ax0.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax0.set_xlabel(r'$t$', fontsize=14)
ax0.set_ylabel(r'$f_{\theta_j(t)} = \lambda \big( \frac{3}{2}(1 - t) \big)^{2\lambda - 1}$', fontsize=14)
ax0.set_ylim(0, 1.0)

ax1.set_title(r"RIs", fontsize=20)
ax1.plot(xs, prior_ri(xs, 1.0), '-', label=r"$\lambda = 1.0$", linewidth=2)
ax1.plot(xs, prior_ri(xs, 0.75), '-', label=r"$\lambda = 0.75$", linewidth=2)
ax1.plot(xs, prior_ri(xs, 0.5), '-', label=r"$\lambda = 0.5$", linewidth=2)
ax1.plot(xs, prior_ri(xs, 0.25), '-', label=r"$\lambda = 0.25$", linewidth=2)
ax1.plot(xs, prior_ri(xs, 0.1), '-', label=r"$\lambda = 0.1$", linewidth=2)
ax1.plot(xs, prior_ri(xs, 0.05), '-', label=r"$\lambda = 0.05$", linewidth=2)
ax1.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

ax1.set_ylim(0, 1.0)
ax1.set_xlabel(r'$t$', fontsize=14)
ax1.set_ylabel(r"$f_{\theta_j}(t) = \lambda \bigg( \frac{(1 + y)^2}{W[y] + 1} \bigg) \bigg( \frac{W[y]}{y} \bigg)^{2 + \lambda}$, where $y = \frac{2}{3}(1-t)^{-1} - 1$", fontsize=14)
ax1.legend(frameon=False, ncol=1, fontsize=14)


for i, ax in enumerate(axs):
	ax.tick_params(axis='x', labelsize=11)
	ax.tick_params(axis='y', labelsize=11)
	ax.get_xaxis().tick_bottom() 
	ax.get_yaxis().tick_left() 
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)


gs.tight_layout(fig, rect=[0, 0.0, 1, 1])
plt.savefig("prior.pdf", format='pdf', dpi=300)
