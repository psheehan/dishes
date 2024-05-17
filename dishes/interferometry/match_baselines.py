from .libinterferometry import Visibilities
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

def match_baselines(vis1, vis2, fcut=0.25, plot=None, ndim=1, iepoch=0,\
        jepoch=1):
    dclose1 = []
    dclose2 = []

    for i in range(vis1.uvdist.size):
        dclose1.append((((vis1.u[i] - vis2.u)**2 + (vis1.v[i] - \
                vis2.v)**2)**0.5 / vis1.uvdist[i]).min())

    for i in range(vis2.uvdist.size):
        dclose2.append((((vis2.u[i] - vis1.u)**2 + (vis2.v[i] - \
                vis1.v)**2)**0.5 / vis2.uvdist[i]).min())

    good1 = np.array(dclose1) < 0.2
    good2 = np.array(dclose2) < 0.2

    vis1_close =  uv.Visibilities(vis1.u[good1], vis1.v[good1], vis1.freq, \
            vis1.real[good1,:], vis1.imag[good1,:], vis1.weights[good1,:])
    vis2_close =  uv.Visibilities(vis2.u[good2], vis2.v[good2], vis2.freq, \
            vis2.real[good2,:], vis2.imag[good2,:], vis2.weights[good2,:])

    # Put the data together.

    if ndim == 1:
        data1 = np.log10(vis1_close.uvdist)
        data2 = np.log10(vis2_close.uvdist)

        weights1 = vis1_close.weights.sum(axis=1)
        weights2 = vis2_close.weights.sum(axis=1)

        x_eval = np.linspace(min(np.log10(vis1_close.uvdist).min(), \
                np.log10(vis2_close.uvdist).min()), \
                max(np.log10(vis1_close.uvdist).max(), \
                np.log10(vis2_close.uvdist).max()), 500)
    elif ndim == 2:
        data1 = np.concatenate((vis1_close.u[np.newaxis,:], \
                vis1_close.v[np.newaxis,:]))
        data2 = np.concatenate((vis2_close.u[np.newaxis,:], \
                vis2_close.v[np.newaxis,:]))

        u, v = np.meshgrid(np.linspace(-vis1_close.uvdist.max(), \
                vis1_close.uvdist.max(), 100), \
                np.linspace(-vis2_close.uvdist.max(), \
                vis2_close.uvdist.max(), 100))
        x_eval = np.vstack([u.ravel(), v.ravel()])

    # Generate a KDE for both datasets.

    kernel1 = scipy.stats.gaussian_kde(data1, weights=weights1)
    kernel2 = scipy.stats.gaussian_kde(data2, weights=weights2)

    # For the purposes of plotting, evaluate on a uniform array of baselines.

    kde1 = kernel1(x_eval) / kernel1(x_eval).max()
    kde2 = kernel2(x_eval) / kernel2(x_eval).max()

    inv_kde2 = (1./kde2) / (1./kde2[kde2 > 0.2]).max()

    selection_prob = kde1 * inv_kde2
    selection_prob /= selection_prob[kde2 > 0.2].max()

    # Now do the real calculation to evaluate which baselines from vis2_close 
    # to keep.

    vis2_close_kde1 = kernel1(data2) / kernel1(data2).max()
    vis2_close_kde2 = kernel2(data2) / kernel2(data2).max()

    vis2_close_inv_kde2 = (1./vis2_close_kde2) / (1./vis2_close_kde2).max()

    vis2_close_selection_prob = vis2_close_kde1 * vis2_close_inv_kde2
    vis2_close_selection_prob /= vis2_close_selection_prob[vis2_close_kde2 > \
            0.2].max()

    good = np.random.uniform(size=vis2_close.uvdist.size) < \
            vis2_close_selection_prob

    #print(good.sum() / good.size)

    # Now that we know which baselines to keep, regenerate the kernel for 
    # vis2_close to evaluate how it compares to vis1_close.

    if ndim == 1:
        kernel2_new = scipy.stats.gaussian_kde(data2[good], \
                weights=weights2[good])
    else:
        kernel2_new = scipy.stats.gaussian_kde(data2[:,good])

    kde2_new = kernel2_new(x_eval) / kernel2_new(x_eval).max()

    if ndim == 1:
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(9,3))

        ax[0].hist(dclose1, 30, alpha=0.25, density=True)
        ax[0].hist(dclose2, 30, alpha=0.25, density=True)
        ax[0].set_xlim(0.,min(1., max(max(dclose1), max(dclose2))*1.1))

        ax[0].set_xlabel("$d_{close}$")

        ax[1].plot(x_eval-3, kde1, \
                label="Epoch {0:d}".format(iepoch)+" post $d_{close} < 0.1$")
        ax[1].plot(x_eval-3, kde2, \
                label="Epoch {0:d}".format(jepoch)+" post $d_{close} < 0.1$")
        ax[1].plot(x_eval-3, kde2_new, color=u'#d62728', \
                label="Epoch {0:d} matched to Epoch {1:d}".format(jepoch,iepoch))

        ax[1].set_ylim(0.,1.3)

        ax[1].set_xlabel("log$_{10}$ Baseline (k$\lambda$)")

        ax[1].legend(framealpha=0.25, prop={'size':'small'}, loc='upper left')

        ax[2].scatter(vis1.u/1000, vis1.v/1000, s=1, alpha=0.05, color="gray")
        ax[2].scatter(vis2.u/1000, vis2.v/1000, s=1, alpha=0.05, color="gray")
        ax[2].scatter(vis1_close.u/1000, vis1_close.v/1000, s=1, alpha=0.75, \
                label="Epoch {0:d}, final".format(iepoch), marker=".")
        ax[2].scatter(vis2_close.u[good]/1000, vis2_close.v[good]/1000, s=1, \
                alpha=1., label="Epoch {0:d}, final".format(jepoch), \
                color=u'#d62728', marker=".")

        ax[2].legend(framealpha=0.25, prop={'size':'small'}, loc='lower left')

        lim = max(vis1_close.uvdist.max(), vis2_close.uvdist.max())/1000
        ax[2].set_xlim(-lim, lim)
        ax[2].set_ylim(-lim, lim)

        ax[2].set_xlabel("u (k$\lambda$)")
        ax[2].set_ylabel("v (k$\lambda$)")

        fig.tight_layout()

        if plot == None:
            plt.show()
        else:
            fig.savefig(plot)
    elif ndim == 2:
        fig, ax = plt.subplots(nrows=1, ncols=1)

        ax.imshow(kde1.reshape(u.shape), origin="lower", \
                interpolation="nearest")
        plt.show()

    return vis1_close, uv.Visibilities(vis2_close.u[good], vis2_close.v[good], \
            vis2_close.freq, vis2_close.real[good,:], vis2_close.imag[good,:], \
            vis2_close.weights[good,:])
