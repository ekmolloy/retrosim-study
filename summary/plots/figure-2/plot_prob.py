import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from scipy import special
import sys

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{helvet} \usepackage{sfmath}']  


def ml_brln_r_inv(x):
    ex = numpy.exp(-1.0 * x)
    return (((1. / 3.) * ex) + x) / (x + ex)


def ml_brln_r(t):
    # y = (t - (1./3)) / (1. - t)
    if t <= 1. / 3:
        return 0
    elif t == 1.:
        return 10
    y = ((2. / 3) / (1. - t)) - 1
    x = special.lambertw(y, k=0, tol=1e-08)
    return x.real


def ml_brln_g_inv(x):
    return 1.0 - (2.0 / 3.0) * numpy.exp(-1.0 * x)


def ml_brln_g(t):
    if t <= 1. / 3:
        return 0
    if t == 1:
        return 10
    return -1.0 * numpy.log(1.5 * (1.0 - t))


def smallangle(x):
    return (1.0 + 2.0*x) / 3.0


xs = numpy.array([0.001, 0.0025, 0.005, 0.0075,
                 0.01, 0.025, 0.05, 0.075,
                 0.1, 0.25, 0.5, 0.75,
                 1.0, 2.5, 3.0])


fig = plt.figure(figsize=(6, 4))
gs = gridspec.GridSpec(1,1)
ax0 = plt.subplot(gs[0])

axs = [ax0]

xs = numpy.arange(0.0, 5.0, 0.001)

ax0.plot(xs, ml_brln_g_inv(xs), 'k--', label=r"Gene Trees: $p(\tau) = 1 - \frac{2}{3}e^{-\tau}$", linewidth=2)
ax0.plot(xs, ml_brln_r_inv(xs), 'k-', label=r"RIs: $p(\tau) = (\tau + \frac{1}{3}e^{-\tau}) / (\tau + e^{-\tau})$", linewidth=2)
ax0.plot(xs, smallangle(xs), 'k:', label=r"Small Angle: $p(\tau) = \frac{1}{3} + \frac{2}{3}\tau$", linewidth=2)

print("%f vs. %f" % (ml_brln_g_inv(0.1), ml_brln_r_inv(0.1)))

#ax0.axhline(0.4, color='gray', linewidth=1, linestyle='-.')
#ax0.axhline(1.0, color='gray', linewidth=1, linestyle='-.')
ax0.axvline(0.3, color='gray', linewidth=1, linestyle='-.', label=r"$\tau = 0.3$")
#ax0.axvline(4.0, color='gray', linewidth=1, linestyle='-.')
ax0.set_xlabel(r'Internal Branch Length $\tau$ in CUs', fontsize=12)
ax0.set_ylabel(r'Probability of Dominant Quartet', fontsize=12)
ax0.set_ylim(0.33, 1.0)

for i, ax in enumerate(axs):
	ax.tick_params(axis='x', labelsize=11)
	ax.tick_params(axis='y', labelsize=11)
	ax.get_xaxis().tick_bottom() 
	ax.get_yaxis().tick_left() 
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)

ax0.legend(frameon=False, ncol=1, fontsize=12)

gs.tight_layout(fig, rect=[0, 0.0, 1, 1])
plt.savefig("prob.pdf", format='pdf', dpi=300)
