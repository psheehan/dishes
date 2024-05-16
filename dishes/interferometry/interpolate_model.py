import numpy
import scipy.interpolate
import astropy.units as units
from .libinterferometry import Visibilities
from galario import double
try:
    import trift
except:
    pass

def interpolate_model(u, v, freq, model, nthreads=1, dRA=0., dDec=0., \
        code="galario", nxy=1024, dxy=0.01):

    if code == "galario":
        real = []
        imag = []

        double.threads(nthreads)

        dxy = (model.x[1] - model.x[0])*units.arcsec.to(units.radian)

        for i in range(len(model.freq)):
            vis = double.sampleImage(model.image[::-1,:,i,0].copy(order='C'), \
                    dxy, u, v, dRA=dRA*units.arcsec.to(units.radian), dDec=dDec*units.arcsec.to(units.radian))

            real.append(vis.real.reshape((u.size,1)))
            imag.append(-vis.imag.reshape((u.size,1)))

        real = numpy.concatenate(real, axis=1)
        imag = numpy.concatenate(imag, axis=1)

    elif code == "galario-unstructured":
        real = []
        imag = []

        double.threads(nthreads)

        vol = None
        for i in range(len(model.freq)):
            vis = double.sampleUnstructuredImage(model.x*units.arcsec.to(units.radian), \
                    -model.y*units.arcsec.to(units.radian), model.image[:,i].copy(order='C'), nxy, \
                    dxy*units.arcsec.to(units.radian), u, v, dRA=dRA*units.arcsec.to(units.radian), dDec=dDec*units.arcsec.to(units.radian))

            real.append(vis.real.reshape((u.size,1)))
            imag.append(-vis.imag.reshape((u.size,1)))

        real = numpy.concatenate(real, axis=1)
        imag = numpy.concatenate(imag, axis=1)

    elif code == "trift":
        vis = trift.cpu.trift(model.x*units.arcsec.to(units.radian), model.y*units.arcsec.to(units.radian), \
                model.image, u, v, dRA*units.arcsec.to(units.radian), dDec*units.arcsec.to(units.radian), \
                nthreads=nthreads, mode="extended")

        real, imag = vis.real, vis.imag

    return Visibilities(u, v, freq, real, imag, numpy.ones(real.shape))
