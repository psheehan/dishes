import numpy
cimport numpy
import h5py
import astropy
cimport cython

cdef class VisibilitiesObject:
    cdef public numpy.ndarray u, v, freq, real, imag, weights, \
            uvdist, amp, phase
    cdef public numpy.ndarray baseline
    cdef public str array_name

    def __init__(self, numpy.ndarray[double, ndim=1] u=None, \
            numpy.ndarray[double, ndim=1] v=None, \
            numpy.ndarray[double, ndim=1] freq=None, \
            numpy.ndarray[double, ndim=2] real=None, \
            numpy.ndarray[double, ndim=2] imag=None, \
            numpy.ndarray[double, ndim=2] weights=None, \
            baseline=None, array_name="CARMA"):

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


class Visibilities(VisibilitiesObject):

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

    def __reduce__(self):
        return (self.rebuild, (self.u, self.v, self.freq, self.real, self.imag,\
                self.weights))

    def rebuild(u, v, freq, real, imag, weights):
        return Visibilities(u, v, freq, real, imag, weights)

    def read(self, filename=None, usefile=None):
        if (usefile == None):
            f = h5py.File(filename, "r")
        else:
            f = usefile

        u, v, freq, real, imag, weights = None, None, None, None, None, None

        if ('u' in f) and ('v' in f):
            u = f['u'].value.astype(numpy.double)
            v = f['v'].value.astype(numpy.double)
        if ('freq' in f):
            freq = f['freq'].value.astype(numpy.double)
        if ('real' in f) and ('imag' in f):
            real = f['real'].value.astype(numpy.double)
            imag = f['imag'].value.astype(numpy.double)

            if ('weights' in f):
                weights = f['weights'].value.astype(numpy.double)
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
            u_dset = f.create_dataset("u", self.u.shape, dtype='float32')
            u_dset[...] = self.u
            v_dset = f.create_dataset("v", self.v.shape, dtype='float32')
            v_dset[...] = self.v
        if (type(self.freq) != type(None)):
            freq_dset = f.create_dataset("freq", self.freq.shape, \
                    dtype='float32')
            freq_dset[...] = self.freq
        if (type(self.real) != type(None)) and (type(self.imag) != type(None)):
            real_dset = f.create_dataset("real", self.real.shape, \
                    dtype='float32')
            real_dset[...] = self.real
            imag_dset = f.create_dataset("imag", self.imag.shape, \
                    dtype='float32')
            imag_dset[...] = self.imag

            if type(self.weights) != type(None):
                if numpy.product(self.weights == numpy.ones(self.real.shape)) \
                        == 0:
                    weights_dset = f.create_dataset("weights", \
                            self.weights.shape, dtype='float32')
                    weights_dset[...] = self.weights

        if (usefile == None):
            f.close()

def average(data, gridsize=256, binsize=None, radial=False, log=False, \
        logmin=None, logmax=None, mfs=False, mode="continuum"):

    cdef numpy.ndarray[double, ndim=3] new_real, new_imag, new_weights
    cdef numpy.ndarray[unsigned int, ndim=1] i, j
    cdef unsigned int k, l

    
    cdef numpy.ndarray[double, ndim=1] u, v, freq, uvdist
    cdef numpy.ndarray[double, ndim=2] real, imag, weights

    if mfs:
        vis = freqcorrect(data)
        u = vis.u.copy()
        v = vis.v.copy()
        uvdist = vis.uvdist.copy()
        freq = vis.freq.copy()
        real = vis.real.copy()
        imag = vis.imag.copy()
        weights = vis.weights.copy()
    else:
        u = data.u.copy()
        v = data.v.copy()
        uvdist = data.uvdist.copy()
        freq = data.freq.copy()
        real = data.real.copy()
        imag = data.imag.copy()
        weights = data.weights.copy()
    
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

    cdef int nuv = u.size
    cdef int nfreq = freq.size
    
    cdef int nchannels
    if mode == "continuum":
        nchannels = 1
    elif mode == "spectralline":
        nchannels = freq.size

    # Average over the U-V plane by creating bins to average over.
    
    if radial:
        if log:
            temp = numpy.linspace(numpy.log10(logmin), numpy.log10(logmax), \
                    gridsize+1)
            new_u = 10**((temp[1:] + temp[0:-1])/2)
        else:
            new_u = numpy.linspace(binsize/2,(gridsize-0.5)*binsize,gridsize)
        new_u = new_u.reshape((1,gridsize))
        new_v = numpy.zeros((1,gridsize))
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

        new_u, new_v = numpy.meshgrid(uu, vv)
        new_real = numpy.zeros((gridsize,gridsize,nchannels))
        new_imag = numpy.zeros((gridsize,gridsize,nchannels))
        new_weights = numpy.zeros((gridsize,gridsize,nchannels))

        if gridsize%2 == 0:
            i = numpy.round(u/binsize+gridsize/2.).astype(numpy.uint32)
            j = numpy.round(v/binsize+gridsize/2.).astype(numpy.uint32)
        else:
            i = numpy.round(u/binsize+(gridsize-1)/2.).astype(numpy.uint32)
            j = numpy.round(v/binsize+(gridsize-1)/2.).astype(numpy.uint32)
    
    for k in range(nuv):
        if mode == "continuum":
            for l in range(nfreq):
                new_real[j[k],i[k],0] += real[k,l]*weights[k,l]
                new_imag[j[k],i[k],0] += imag[k,l]*weights[k,l]
                new_weights[j[k],i[k],0] += weights[k,l]
        elif mode == "spectralline":
            for l in range(nfreq):
                new_real[j[k],i[k],l] += real[k,l]*weights[k,l]
                new_imag[j[k],i[k],l] += imag[k,l]*weights[k,l]
                new_weights[j[k],i[k],l] += weights[k,l]
    
    good_data = new_weights != 0.0
    good_data = numpy.all(good_data, axis=2)
    good_data = numpy.dstack([good_data for m in range(nchannels)])
    new_u = new_u[good_data[:,:,0]]
    new_v = new_v[good_data[:,:,0]]
    
    if mode == "continuum":
        freq = numpy.array([data.freq.sum()/data.freq.size])
    
    return Visibilities(new_u,new_v,freq,\
            (new_real[good_data] / new_weights[good_data]).reshape(\
            (new_u.size,nchannels)),\
            (new_imag[good_data] / new_weights[good_data]).reshape(\
            (new_u.size,nchannels)), \
            new_weights[good_data].reshape((new_u.size,nchannels)))

@cython.boundscheck(False)
def grid(data, gridsize=256, binsize=2000.0, convolution="pillbox", \
        mfs=False, channel=None, imaging=False, weighting="natural", \
        robust=2, mode="continuum"):
    
    cdef numpy.ndarray[double, ndim=1] u, v, freq
    cdef numpy.ndarray[double, ndim=2] real, imag, weights, new_u, new_v
    cdef numpy.ndarray[double, ndim=3] new_real, new_imag, new_weights
    cdef numpy.ndarray[unsigned int, ndim=2] i, j
    cdef unsigned int k, l, m, n, ninclude_min, ninclude_max, lmin, lmax, \
            mmin, mmax, ll, mm, ninclude
    cdef double mean_freq, inv_freq

    if mfs:
        vis = freqcorrect(data)
        u = vis.u.copy()
        v = vis.v.copy()
        freq = vis.freq.copy()
        real = vis.real.copy()
        imag = vis.imag.copy()
        weights = vis.weights.copy()
    else:
        u = data.u.copy()
        v = data.v.copy()
        if channel != None:
            freq = numpy.array([data.freq[channel]])
            real = data.real[:,channel].copy().reshape((data.real.shape[0],1))
            imag = data.imag[:,channel].copy().reshape((data.real.shape[0],1))
            weights = data.weights[:,channel].copy(). \
                    reshape((data.real.shape[0],1))
        else:
            freq = data.freq.copy()
            real = data.real.copy()
            imag = data.imag.copy()
            weights = data.weights.copy()
    
    # Set the weights equal to 0 when the point is flagged (i.e. weight < 0)
    weights = numpy.where(weights < 0, 0.0, weights)
    # Set the weights equal to 0 when the real and imaginary parts are both 0
    weights[(real==0) & (imag==0)] = 0.0

    # Set some parameter numbers for future use.
    
    cdef double convolve
    cdef double inv_binsize = 1. / binsize
    cdef int nuv = u.size
    cdef int nfreq = freq.size

    cdef int nchannels
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

    if gridsize%2 == 0:
        for k in range(nuv):
            for n in range(nfreq):
                i[k,n] = numpy.round(u[k]*freq[n]*inv_freq/binsize+\
                        gridsize/2.).astype(numpy.uint32)
                j[k,n] = numpy.round(v[k]*freq[n]*inv_freq/binsize+\
                        gridsize/2.).astype(numpy.uint32)
    else:
        for k in range(nuv):
            for n in range(nfreq):
                i[k,n] = numpy.round(u[k]*freq[n]*inv_freq/binsize+ \
                        (gridsize-1)/2.).astype(numpy.uint32)
                j[k,n] = numpy.round(v[k]*freq[n]*inv_freq/binsize+ \
                        (gridsize-1)/2.).astype(numpy.uint32)
    
    if convolution == "pillbox":
        convolve_func = ones
        ninclude = 3
    elif convolution == "expsinc":
        convolve_func = exp_sinc
        ninclude = 7

    ninclude_min = numpy.uint32((ninclude-1)*0.5)
    ninclude_max = numpy.uint32((ninclude-1)*0.5)

    # Now actually go through and calculate the new visibilities.

    for k in range(nuv):
        for n in range(nfreq):
            if ninclude_min > j[k,n]:
                lmin = 0
            else:
                lmin = j[k,n] - ninclude_min

            lmax = int_min(j[k,n]+ninclude_max+1, gridsize)

            if ninclude_min > i[k,n]:
                mmin = 0
            else:
                mmin = i[k,n] - ninclude_min

            mmax = int_min(i[k,n]+ninclude_max+1, gridsize)

            for l in range(lmin, lmax):
                for m in range(mmin, mmax):
                    convolve = convolve_func( (u[k]-new_u[l,m])*inv_binsize, \
                            (v[k] - new_v[l,m]) * inv_binsize)

                    if mode == "continuum":
                        new_real[l,m,0] += real[k,n]*weights[k,n]*convolve
                        new_imag[l,m,0] += imag[k,n]*weights[k,n]*convolve
                        new_weights[l,m,0] += weights[k,n]*convolve
                    elif mode == "spectralline":
                        new_real[l,m,n] += real[k,n]*weights[k,n]*convolve
                        new_imag[l,m,n] += imag[k,n]*weights[k,n]*convolve
                        new_weights[l,m,n] += weights[k,n]*convolve

    # If we have a special weighting scheme, fix the sums..

    if weighting == "uniform":
        for l in range(gridsize):
            for m in range(gridsize):
                for n in range(nchannels):
                    if new_weights[l,m,n] == 0:
                        continue

                    new_real[l,m,n] /= new_weights[l,m,n]
                    new_imag[l,m,n] /= new_weights[l,m,n]

    elif weighting == "robust":
        f = numpy.sqrt((5 * 10**(-robust))**2 / \
                ((new_weights**2).sum()/weights.sum()))

        for l in range(gridsize):
            for m in range(gridsize):
                for n in range(nchannels):
                    if new_weights[l,m,n] == 0:
                        continue

                    new_real[l,m,n] /= (1+new_weights[l,m,n]*f**2)
                    new_imag[l,m,n] /= (1+new_weights[l,m,n]*f**2)

    if not imaging:
        good = new_weights > 0
        new_real[good] = new_real[good] / new_weights[good]
        new_imag[good] = new_imag[good] / new_weights[good]

    if mode == "continuum":
        freq = numpy.array([data.freq.sum()/data.freq.size])
    
    return Visibilities(new_u.reshape(gridsize**2), new_v.reshape(gridsize**2),\
            freq, new_real.reshape((gridsize**2,nchannels)), \
            new_imag.reshape((gridsize**2,nchannels)), \
            new_weights.reshape((gridsize**2,nchannels)))

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b
cdef inline double int_abs(double a): return -a if a < 0 else a

def test(x):
    return expsinc(x)

cdef unsigned int find_in_arr(double val, numpy.ndarray[double, ndim=1] arr, \
        unsigned int n):

    cdef unsigned int lmin = 0
    cdef unsigned int lmax = n-1
    cdef int not_found = 1
    cdef unsigned int l, ltest

    while not_found:
        ltest = <unsigned int>(lmax-lmin)/2+lmin

        if ((val >= arr[ltest]) and (val <= arr[ltest+1])):
            l = ltest
            not_found = 0
        else:
            if (val < arr[ltest]):
                lmax = ltest
            else:
                lmin = ltest

    return l

cdef numpy.ndarray x_arr = numpy.linspace(-10,10,301)
cdef numpy.ndarray xx_arr = (x_arr[0:-1] + x_arr[1:])/2

cdef numpy.ndarray expsinc_arr = numpy.exp(-(xx_arr/2.52)**2) * \
        numpy.sinc(xx_arr/1.55)

cdef double expsinc(double a):

    cdef unsigned int i = find_in_arr(a, x_arr, 301)

    return expsinc_arr[i]

cdef double exp_sinc(double u, double v):
    
    cdef double inv_alpha1 = 1. / 1.55
    #cdef double inv_alpha1 = 0.6451612903225806
    cdef double inv_alpha2 = 1. / 2.52
    #cdef double inv_alpha2 = 0.3968253968253968
    cdef int m = 6
    
    if (int_abs(u) >= m * 0.5) or (int_abs(v) >= m * 0.5):
        return 0.

    """
    cdef double arr = numpy.sinc(u * inv_alpha1) * \
            numpy.sinc(v * inv_alpha1)* \
            numpy.exp(-1 * (u * inv_alpha2)**2) * \
            numpy.exp(-1 * (v * inv_alpha2)**2)
    """
    cdef double arr = expsinc(u) * expsinc(v)

    return arr

cdef double ones(double u, double v):
    
    cdef int m = 1

    if (int_abs(u) >= m * 0.5) or (int_abs(v) >= m * 0.5):
        return 0.

    cdef double arr = 1.0

    return arr

@cython.boundscheck(False)
def freqcorrect(data, freq=None):
    cdef numpy.ndarray[double, ndim=1] new_u, new_v, new_freq
    cdef numpy.ndarray[double, ndim=2] new_real, new_imag, new_weights

    if freq != None:
        new_freq = numpy.array([freq])
    else:
        new_freq = numpy.array([data.freq.mean()])

    new_u = numpy.empty((data.real.size,))
    new_v = numpy.empty((data.real.size,))
    new_real = numpy.empty((data.real.size, 1))
    new_imag = numpy.empty((data.real.size, 1))
    new_weights = numpy.empty((data.real.size, 1))

    cdef int nuv = data.u.size
    cdef int nfreq = data.freq.size
    cdef unsigned int i, j, index
    cdef double inv_freq = 1./new_freq[0]
    cdef numpy.ndarray[double, ndim=1] scale = data.freq * inv_freq

    for i in range(nuv):
        for j in range(nfreq):
            index = <unsigned int>(j + i*nfreq)
            new_u[index] = data.u[<unsigned int>i] * scale[<unsigned int>j]
            new_v[index] = data.v[<unsigned int>i] * scale[<unsigned int>j]

    new_real = data.real.reshape((data.real.size,1))
    new_imag = data.imag.reshape((data.imag.size,1))
    new_weights = data.weights.reshape((data.weights.size,1))

    return Visibilities(new_u, new_v, new_freq, new_real, new_imag, \
            new_weights)

def chisq(data, model):
    chi_squared = chisq_calc(data.real, data.imag, data.weights, model.real, \
            model.imag, data.real.size)

    return chi_squared

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cdef float chisq_calc(numpy.ndarray[double, ndim=2] data_real, \
        numpy.ndarray[double, ndim=2] data_imag, \
        numpy.ndarray[double, ndim=2] data_weights, \
        numpy.ndarray[double, ndim=2] model_real, \
        numpy.ndarray[double, ndim=2] model_imag, int nuv):
    cdef double chisq = 0
    cdef double diff1, diff2
    cdef unsigned int i

    for i in range(nuv):
        diff1 = data_real[<unsigned int>i,0] - model_real[<unsigned int>i,0]
        diff2 = data_imag[<unsigned int>i,0] - model_imag[<unsigned int>i,0]
        chisq += (diff1*diff1 + diff2*diff2) * data_weights[<unsigned int>i,0]

    return chisq
