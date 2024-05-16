import numpy
import h5py
import astropy
import time
import numba

class Visibilities:
    def __init__(self, u=None, v=None, freq=None, real=None, imag=None, \
            weights=None, baseline=None, array_name="CARMA"):

        if (type(u) != type(None)) and (type(v) != type(None)):
            self.u = u
            self.v = v
            self.uvdist = numpy.sqrt(u**2 + v**2)

        if (type(freq) != type(None)):
            self.freq = freq

        if (type(real) != type(None)) and (type(imag) != type(None)):
            self.real = real
            self.imag = imag
            self.amp = numpy.sqrt(real**2 + imag**2)
            self.phase = numpy.arctan2(imag, real)

            if type(weights) != type(None):
                self.weights = weights
            else:
                self.weights = numpy.ones((self.real.shape[0],self.real.shape[1]))

        self.baseline = baseline
        self.array_name = array_name

    def __reduce__(self):
        return (self.rebuild, (self.u, self.v, self.freq, self.real, self.imag,\
                self.weights, self.baseline, self.array_name))

    def rebuild(u, v, freq, real, imag, weights, baseline, array_name):
        return Visibilities(u, v, freq, real, imag, weights, baseline, \
                array_name)

    def get_baselines(self, num):
        incl = self.baseline == num

        return Visibilities(self.u[incl], self.v[incl], self.freq, \
                self.real[incl,:], self.imag[incl,:], \
                self.weights[incl,:], baseline=self.baseline[incl])

    def set_header(self, header):
        self.header = header

    def asFITS(self):
        hdulist = astropy.fits.HDUList([])

        nvis = self.u.size
        hdu = astropy.fits.PrimaryHDU(numpy.concatenate((\
                self.u.reshape((1,nvis)), self.v.reshape((1,nvis)), \
                self.real.reshape((1,nvis)), self.imag.reshape((1,nvis)), \
                self.weights.reshape((1,nvis))), axis=0))

        if type(self.header) != type(None):
            hdu.header = self.header

        hdulist.append(hdu)

        return hdulist

    def read(self, filename=None, usefile=None):
        if (usefile == None):
            f = h5py.File(filename, "r")
        else:
            f = usefile

        u, v, freq, real, imag, weights = None, None, None, None, None, None

        if ('u' in f) and ('v' in f):
            u = f['u'][...].astype(numpy.double)
            v = f['v'][...].astype(numpy.double)
        if ('freq' in f):
            freq = f['freq'][...].astype(numpy.double)
        if ('real' in f) and ('imag' in f):
            real = f['real'][...].astype(numpy.double)
            imag = f['imag'][...].astype(numpy.double)

            if ('weights' in f):
                weights = f['weights'][...].astype(numpy.double)
            else:
                weights = numpy.ones(real.shape)

        self.__init__(u, v, freq, real, imag, weights)

        if (usefile == None):
            f.close()

    def write(self, filename=None, usefile=None):
        if (usefile == None):
            f = h5py.File(filename, "w")
        else:
            f = usefile

        if (type(self.u) != type(None)) and (type(self.v) != type(None)):
            u_dset = f.create_dataset("u", self.u.shape, dtype='float64')
            u_dset[...] = self.u
            v_dset = f.create_dataset("v", self.v.shape, dtype='float64')
            v_dset[...] = self.v
        if (type(self.freq) != type(None)):
            freq_dset = f.create_dataset("freq", self.freq.shape, \
                    dtype='float64')
            freq_dset[...] = self.freq
        if (type(self.real) != type(None)) and (type(self.imag) != type(None)):
            real_dset = f.create_dataset("real", self.real.shape, \
                    dtype='float64')
            real_dset[...] = self.real
            imag_dset = f.create_dataset("imag", self.imag.shape, \
                    dtype='float64')
            imag_dset[...] = self.imag

            if type(self.weights) != type(None):
                if numpy.product(self.weights == numpy.ones(self.real.shape)) \
                        == 0:
                    weights_dset = f.create_dataset("weights", \
                            self.weights.shape, dtype='float64')
                    weights_dset[...] = self.weights

        if (usefile == None):
            f.close()

def average(data, gridsize=256, binsize=None, radial=False, log=False, \
        logmin=None, logmax=None, mfs=False, mode="continuum"):

    if mfs:
        vis = freqcorrect(data)
        u = vis.u
        v = vis.v
        uvdist = vis.uvdist
        freq = vis.freq
        real = vis.real
        imag = vis.imag
        weights = vis.weights
    else:
        u = data.u.copy()
        v = data.v.copy()
        uvdist = data.uvdist.copy()
        freq = data.freq.copy()
        real = data.real.copy()
        imag = data.imag.copy()
        weights = data.weights
    
    # Set the weights equal to 0 when the point is flagged (i.e. weight < 0)
    weights = numpy.where(weights < 0,0.0,weights)
    # Set the weights equal to 0 when the real and imaginary parts are both 0
    weights[(real == 0) & (imag == 0)] = 0.0
    
    good_data = uvdist != 0.0
    u = u[good_data]
    v = v[good_data]
    uvdist = uvdist[good_data]
    real = real[good_data,:]
    imag = imag[good_data,:]
    weights = weights[good_data,:]

    # Set some parameter numbers for future use.

    nuv = u.size
    nfreq = freq.size
    
    if mode == "continuum":
        nchannels = 1
    elif mode == "spectralline":
        nchannels = freq.size

    # Average over the U-V plane by creating bins to average over.
    
    if radial:
        if log:
            temp = numpy.linspace(numpy.log10(logmin), numpy.log10(logmax), \
                    gridsize+1)
            new_u_final = 10**((temp[1:] + temp[0:-1])/2)
        else:
            new_u_final = numpy.linspace(binsize/2,(gridsize-0.5)*binsize,gridsize)
        new_u_final = new_u_final.reshape((1,gridsize))
        new_v_final = numpy.zeros((1,gridsize))
        new_u = numpy.zeros((1,gridsize,nchannels))
        new_v = numpy.zeros((1,gridsize,nchannels))
        new_real = numpy.zeros((1,gridsize,nchannels))
        new_imag = numpy.zeros((1,gridsize,nchannels))
        new_weights = numpy.zeros((1,gridsize,nchannels))

        if log:
            dtemp = temp[1] - temp[0]
            i = numpy.round((numpy.log10(uvdist)- numpy.log10(logmin))/dtemp - \
                    0.5).astype(numpy.uint32)
            j = numpy.zeros(uvdist.size).astype(numpy.uint32)
        else:
            i = numpy.round(uvdist/binsize).astype(numpy.uint32)
            j = numpy.zeros(uvdist.size).astype(numpy.uint32)
    else:
        if gridsize%2 == 0:
            uu = numpy.linspace(-gridsize*binsize/2, (gridsize/2-1)*binsize, \
                    gridsize)
            vv = numpy.linspace(-gridsize*binsize/2, (gridsize/2-1)*binsize, \
                    gridsize)
        else:
            uu = numpy.linspace(-(gridsize-1)*binsize/2, \
                    (gridsize-1)*binsize/2, gridsize)
            vv = numpy.linspace(-(gridsize-1)*binsize/2, \
                    (gridsize-1)*binsize/2, gridsize)

        #new_u, new_v = numpy.meshgrid(uu, vv)
        new_u = numpy.zeros((gridsize,gridsize,nchannels))
        new_v = numpy.zeros((gridsize,gridsize,nchannels))
        new_real = numpy.zeros((gridsize,gridsize,nchannels))
        new_imag = numpy.zeros((gridsize,gridsize,nchannels))
        new_weights = numpy.zeros((gridsize,gridsize,nchannels))

        if gridsize%2 == 0:
            i = numpy.round(u/binsize+gridsize/2.).astype(numpy.uint32)
            j = numpy.round(v/binsize+gridsize/2.).astype(numpy.uint32)
        else:
            i = numpy.round(u/binsize+(gridsize-1)/2.).astype(numpy.uint32)
            j = numpy.round(v/binsize+(gridsize-1)/2.).astype(numpy.uint32)

    good_i = numpy.logical_and(i >= 0, i < gridsize)
    good_j = numpy.logical_and(j >= 0, j < gridsize)
    good = numpy.logical_and(good_i, good_j)
    if good.sum() < good.size:
        print("WARNING: uv.grid was supplied with a gridsize and binsize that do not cover the full range of the input data in the uv-plane and is cutting baselines that are outside of this grid. Make sure to check your results carefully.")

    u = u[good]
    v = v[good]
    real = real[good,:]
    imag = imag[good,:]
    weights = weights[good,:]

    i = i[good]
    j = j[good]

    if mode == "continuum":
        l = numpy.repeat(0, nfreq)
    else:
        l = numpy.arange(nfreq)

    nuv = u.size

    new_u, new_v, new_real, new_imag, new_weights = grid_engine(u, v, real, \
            imag, weights, i, j, l, new_u, new_v, new_real, new_imag, \
            new_weights, nuv, nfreq)

    good_data = new_weights != 0.0
    new_real[good_data] = new_real[good_data] / new_weights[good_data]
    new_imag[good_data] = new_imag[good_data] / new_weights[good_data]
    if not radial:
        new_u[good_data] = new_u[good_data] / new_weights[good_data]
        new_v[good_data] = new_v[good_data] / new_weights[good_data]
    else:
        new_u = new_u_final
        new_v = new_v_final

    good_data = numpy.any(good_data, axis=2)
    if not radial:
        new_u = (new_u*new_weights).sum(axis=2)[good_data] / \
                new_weights.sum(axis=2)[good_data]
        new_v = (new_v*new_weights).sum(axis=2)[good_data] / \
                new_weights.sum(axis=2)[good_data]
    else:
        new_u = new_u[good_data]
        new_v = new_v[good_data]

    good_data = numpy.dstack([good_data for m in range(nchannels)])

    if mode == "continuum":
        freq = numpy.array([data.freq.sum()/data.freq.size])

    return Visibilities(new_u, new_v, freq, \
            new_real[good_data].reshape((new_u.size,nchannels)),\
            new_imag[good_data].reshape((new_u.size,nchannels)), \
            new_weights[good_data].reshape((new_u.size,nchannels)))

@numba.jit(fastmath=True, nopython=True)
def grid_engine(u, v, real, imag, weights, i, j, l, new_u, new_v, new_real, \
        new_imag, new_weights, nuv, nfreq):

    for k in range(nuv):
        for n in range(nfreq):
            new_u[j[k],i[k],l[n]] += u[k]*weights[k,n]
            new_v[j[k],i[k],l[n]] += v[k]*weights[k,n]
            new_real[j[k],i[k],l[n]] += real[k,n]*weights[k,n]
            new_imag[j[k],i[k],l[n]] += imag[k,n]*weights[k,n]
            new_weights[j[k],i[k],l[n]] += weights[k,n]

    return new_u, new_v, new_real, new_imag, new_weights


def grid(data, gridsize=256, binsize=2000.0, convolution="pillbox", \
        mfs=False, channel=None, imaging=False, weighting="natural", \
        robust=2, npixels=0, mode="continuum"):
    
    if mfs:
        vis = freqcorrect(data)
        u = vis.u
        v = vis.v
        freq = vis.freq
        real = vis.real
        imag = vis.imag
        weights = vis.weights.copy()
    else:
        u = data.u
        v = data.v
        if channel != None:
            freq = numpy.array([data.freq[channel]])
            real = data.real[:,channel].reshape((data.real.shape[0],1))
            imag = data.imag[:,channel].reshape((data.real.shape[0],1))
            weights = data.weights[:,channel]. \
                    reshape((data.real.shape[0],1))
        else:
            freq = data.freq
            real = data.real
            imag = data.imag
            weights = data.weights.copy()
    
    # Set the weights equal to 0 when the point is flagged (i.e. weight < 0)
    weights = numpy.where(weights < 0, 0.0, weights)
    # Set the weights equal to 0 when the real and imaginary parts are both 0
    weights[(real==0) & (imag==0)] = 0.0

    # Set some parameter numbers for future use.
    
    inv_binsize = 1. / binsize
    nuv = u.size
    nfreq = freq.size

    if mode == "continuum":
        nchannels = 1
    elif mode == "spectralline":
        nchannels = freq.size

    # Average over the U-V plane by creating bins to average over.
    
    if gridsize%2 == 0:
        uu = numpy.linspace(-gridsize*binsize/2, (gridsize/2-1)*binsize, \
                gridsize)
        vv = numpy.linspace(-gridsize*binsize/2, (gridsize/2-1)*binsize, \
                gridsize)
    else:
        uu = numpy.linspace(-(gridsize-1)*binsize/2, (gridsize-1)*binsize/2, \
                gridsize)
        vv = numpy.linspace(-(gridsize-1)*binsize/2, (gridsize-1)*binsize/2, \
                gridsize)

    new_u, new_v = numpy.meshgrid(uu, vv)
    new_real = numpy.zeros((gridsize,gridsize,nchannels))
    new_imag = numpy.zeros((gridsize,gridsize,nchannels))
    new_weights = numpy.zeros((gridsize,gridsize,nchannels))

    # Get the indices for binning.

    i = numpy.zeros((nuv, nfreq), dtype=numpy.uint32)
    j = numpy.zeros((nuv, nfreq), dtype=numpy.uint32)

    mean_freq = numpy.mean(freq)
    inv_freq = 1./mean_freq

    if gridsize%2 == 0:
        i = numpy.array(u.reshape((u.size,1))*freq*inv_freq/binsize+\
                gridsize/2., dtype=numpy.uint32)
        j = numpy.array(v.reshape((v.size,1))*freq*inv_freq/binsize+\
                gridsize/2., dtype=numpy.uint32)
    else:
        i = numpy.array(u.reshape((u.size,1))*freq*inv_freq/binsize+ \
                    (gridsize-1)/2., dtype=numpy.uint32)
        j = numpy.array(v.reshape((v.size,1))*freq*inv_freq/binsize+ \
                    (gridsize-1)/2., dtype=numpy.uint32)

    if mode == "continuum":
        f = numpy.repeat(0, nfreq)
    else:
        f = numpy.arange(nfreq)
    
    if convolution == "pillbox":
        convolve_func = ones
        ninclude = 3
    elif convolution == "expsinc":
        convolve_func = exp_sinc
        ninclude = 6

    if ninclude%2 == 0:
        ninclude_min = numpy.uint32(ninclude*0.5-1)
        ninclude_max = numpy.uint32(ninclude*0.5)
    else:
        ninclude_min = numpy.uint32((ninclude-1)*0.5)
        ninclude_max = numpy.uint32((ninclude-1)*0.5)

    # Check whether any data falls outside of the grid, and exclude.

    good_i = numpy.logical_and(i >= 0, i < gridsize)
    good_j = numpy.logical_and(j >= 0, j < gridsize)
    good = numpy.logical_and(good_i, good_j)
    if good.sum() < good.size:
        print("WARNING: uv.grid was supplied with a gridsize and binsize that do not cover the full range of the input data in the uv-plane and is cutting baselines that are outside of this grid. Make sure to check your results carefully.")

    not_good = numpy.logical_not(good)
    weights[not_good] == 0.
    i[not_good], j[not_good] = 0, 0

    # If we are using a non-uniform weighting scheme, adjust the data weights.

    if weighting in ["uniform","superuniform","robust"]:
        binned_weights = numpy.ones((gridsize,gridsize,nchannels))

        npix = npixels
        if weighting == "superuniform":
            npix = 3

        weight_binner(weights, binned_weights, i, j, f, nuv, nfreq, npix, \
                gridsize)

        if weighting in ["uniform","superuniform"]:
            f1 = 0.
            f2 = numpy.ones(nfreq)
        elif weighting == "robust":
            f1 = 1.
            f2 = (5*10**(-robust))**2 / \
                    ((binned_weights**2).sum(axis=(0,1)) / weights.sum(axis=0))

        weight_divider(weights, binned_weights, i, j, f1, f2, nuv, nfreq)

    # Now actually go through and calculate the new visibilities.

    new_real, new_imag, new_weights = gridder(u, v, freq, real, imag, weights, \
            new_u, new_v, new_real, new_imag, new_weights, i, j, f, nuv, nfreq,\
            ninclude_min, ninclude_max, gridsize, convolve_func, inv_freq, \
            inv_binsize)

    # If we are making an image, normalize the weights.

    if imaging:
        for n in range(nchannels):
            new_real[:,:,n] /= new_weights[:,:,n].sum()
            new_imag[:,:,n] /= new_weights[:,:,n].sum()
            new_weights[:,:,n] /= new_weights[:,:,n].sum()
    else:
        good = new_weights > 0
        new_real[good] = new_real[good] / new_weights[good]
        new_imag[good] = new_imag[good] / new_weights[good]

    if mode == "continuum":
        freq = numpy.array([data.freq.sum()/data.freq.size])
    
    return Visibilities(new_u.reshape(gridsize**2), new_v.reshape(gridsize**2),\
            freq, new_real.reshape((gridsize**2,nchannels)), \
            new_imag.reshape((gridsize**2,nchannels)), \
            new_weights.reshape((gridsize**2,nchannels)))

@numba.jit(fastmath=True, nopython=True)
def weight_binner(weights, binned_weights, i, j, f, nuv, nfreq, npix, gridsize):

    for k in range(nuv):
        for n in range(nfreq):
            if npix > j[k,n]:
                lmin = 0
            else:
                lmin = j[k,n] - npix

            lmax = min(j[k,n]+npix+1, gridsize)

            if npix > i[k,n]:
                mmin = 0
            else:
                mmin = i[k,n] - npix

            mmax = min(i[k,n]+npix+1, gridsize)

            for l in range(lmin, lmax):
                for m in range(mmin, mmax):
                    binned_weights[l,m,f[n]] += weights[k,n]

    return binned_weights

@numba.jit(fastmath=True, nopython=True)
def weight_divider(weights, binned_weights, i, j, a, b, nuv, nfreq):

    for k in range(nuv):
        for n in range(nfreq):
            l = j[k,n]
            m = i[k,n]

            weights[k,n] /= (a + b[n] * binned_weights[l,m,n])


@numba.jit(fastmath=True, nopython=True)
def gridder(u, v, freq, real, imag, weights, new_u, new_v, new_real, new_imag, \
        new_weights, i, j, f, nuv, nfreq, ninclude_min, ninclude_max, gridsize,\
        convolve_func, inv_freq, inv_binsize):

    for k in range(nuv):
        for n in range(nfreq):
            if ninclude_min > j[k,n]:
                lmin = 0
            else:
                lmin = j[k,n] - ninclude_min

            lmax = min(j[k,n]+ninclude_max+1, gridsize)

            if ninclude_min > i[k,n]:
                mmin = 0
            else:
                mmin = i[k,n] - ninclude_min

            mmax = min(i[k,n]+ninclude_max+1, gridsize)

            for l in range(lmin, lmax):
                for m in range(mmin, mmax):
                    convolve = convolve_func( (u[k]*freq[n]*inv_freq-\
                            new_u[l,m])*inv_binsize, \
                            (v[k]*freq[n]*inv_freq - new_v[l,m]) * inv_binsize)

                    new_real[l,m,f[n]] += real[k,n]*weights[k,n]*convolve
                    new_imag[l,m,f[n]] += imag[k,n]*weights[k,n]*convolve
                    new_weights[l,m,f[n]] += weights[k,n]*convolve

    return new_real, new_imag, new_weights

@numba.jit(fastmath=True, nopython=True)
def sinc(x):

    xp = x * numpy.pi

    return 1. - xp**2/6. + xp**4/120. - xp**6/5040. + xp**8/362880. - \
            xp**10/39916800. + xp**12/6227020800. - xp**14/1307674368000. + \
            xp**16/355687428096000.

@numba.jit(fastmath=True, nopython=True)
def exp(x):

    return 1 + x + x**2/2. + x**3/6. + x**4/24. + x**5/120.

@numba.jit(fastmath=True, nopython=True)
def exp_sinc(u, v):
    
    inv_alpha1 = 1. / 1.55
    inv_alpha2 = 1. / 2.52
    norm = 2.350016262343186
    m = 6
    
    if (abs(u) >= m * 0.5) or (abs(v) >= m * 0.5):
        return 0.

    arr = sinc(u * inv_alpha1) * \
            sinc(v * inv_alpha1) * \
            exp(-1 * (u * inv_alpha2)**2) * \
            exp(-1 * (v * inv_alpha2)**2) / norm

    return arr

@numba.jit(fastmath=True, nopython=True)
def ones(u, v):
    
    m = 1

    if (abs(u) >= m * 0.5) or (abs(v) >= m * 0.5):
        return 0.

    arr = 1.0

    return arr

def freqcorrect(data, freq=None):

    if freq != None:
        new_freq = numpy.array([freq])
    else:
        new_freq = numpy.array([data.freq.mean()])

    inv_freq = 1./new_freq[0]
    scale = data.freq * inv_freq

    new_u = (data.u.reshape((data.u.size,1))*scale).reshape((data.real.size,))
    new_v = (data.v.reshape((data.v.size,1))*scale).reshape((data.real.size,))

    new_real = data.real.reshape((data.real.size,1))
    new_imag = data.imag.reshape((data.imag.size,1))
    new_weights = data.weights.reshape((data.weights.size,1))

    return Visibilities(new_u, new_v, new_freq, new_real, new_imag, \
            new_weights)

def chisq(data, model):
    chi_squared = chisq_calc(data.real, data.imag, data.weights, model.real, \
            model.imag, data.real.size)

    return chi_squared

def chisq_calc(data_real, \
        data_imag, \
        data_weights, \
        model_real, \
        model_imag, nuv):
    chisq = 0

    for i in range(nuv):
        diff1 = data_real[i,0] - model_real[i,0]
        diff2 = data_imag[i,0] - model_imag[i,0]
        chisq += (diff1*diff1 + diff2*diff2) * data_weights[i,0]

    return chisq
