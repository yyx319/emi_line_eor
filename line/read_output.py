import numpy as np
from scipy.io import FortranFile as ff


def read_flux(filepath, ndir):
    file = ff(filepath)
    flux_a = []
    for i in range(ndir):
        aper, flux = file.read_reals(dtype='f8') 
        flux_a.append( flux )
    flux_a = np.array(flux_a)
    return flux_a

def read_spectrum(filepath, ndir):
    file = ff(filepath)
    spectrum_a = []
    for i in range(ndir):
        res = file.read_ints()
        file.read_reals(dtype='f8')
        spectrum = file.read_reals(dtype='f8') 
        spectrum_a.append( spectrum )
    spectrum_a = np.array(spectrum_a)
    return spectrum_a

def read_image(filepath, ndir):
    file = ff(filepath)
    image_a = []
    for i in range(ndir):
        res = file.read_ints()
        file.read_reals(dtype='f8')
        file.read_reals(dtype='f8')
        image = file.read_reals(dtype='f8') 
        image_a.append( image )
    image_a = np.array(image_a)
    return image_a

def read_output(file_fd, ndir, mech='col'):
    flux_a = read_flux('%s/%s_flux.00000'%(file_fd, mech), ndir)
    spectrum_a = read_spectrum('%s/%s_spectrum.00000'%(file_fd, mech), ndir)
    image_a = read_image('%s/%s_image.00000'%(file_fd, mech), ndir)
    return flux_a, spectrum_a, image_a

