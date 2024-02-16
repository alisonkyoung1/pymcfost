#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy import units
import os
from astropy import constants
import numpy as np

### GET MCFOST OUTPUT ###

def get_mcfost_output(datadir):
    """

    Function take output of mcfost, all in one file, and converts it into
    required parameters for astrochem input file.

    """
    from astropy.io import fits

    datadir = os.path.normpath(os.path.expanduser(datadir))

    ## Gas density ##
    hdu = fits.open(datadir + "/gas_density.fits.gz")
    gas_mass_density = (hdu[0].data * units.g *
                        units.cm**-3)   # g cm^-3
    gas_number_density = gas_mass_density / (2.0 * constants.u.to('g'))  # cm-3
    hdu.close()

    ## Dust temperature ##
    hdu = fits.open(datadir + "/Temperature.fits.gz")
    temperature = hdu[0].data * units.K
    hdu.close()

    ## UV field ##
    hdu = fits.open(datadir + "/UV_field.fits.gz")
    chi = hdu[0].data  # Habing units; convert into Draine units
    hdu.close()

    ## Grain mass ##
    hdu = fits.open(datadir + "/dust_mass_density.fits.gz")
    dust_mass_density = (hdu[0].data * units.g *
                         units.cm**-3)    # g cm^-3
    dust_gas_mass_ratio = dust_mass_density / gas_mass_density
    hdu.close()

    ## Grain size ##
    hdu = fits.open(datadir + "/grain_sizes.fits.gz")
    grain_sizes = hdu[0].data * units.um  # um
    hdu.close()
    hdu = fits.open(datadir + "/dust_particle_density.fits.gz")
    dust_number_density = (hdu[0].data *
                           units.m**-3)   # m^-3 per grain size bin
    hdu.close()
    average_grain_size = np.sqrt(sum(grain_sizes**2 *
                                     dust_number_density) /
                                 sum(dust_number_density)).to('um')

    # X-ray ionization (Bai & Goodman 2009) ##
    zeta1 = 6e-12 * units.s**-1  # s-1 (Tx = 3 keV)
    N1 = 1.5e21 * units.cm**-2   # cm
    zeta2 = 1e-15 * units.s**-1  # s-1
    N2 = 7e23 * units.cm**-2     # cm
    Lx29 = 5.                    # 10^29 erg s-1
    alpha = 0.4
    beta = 0.65
    hdu = fits.open(datadir + "/Column_density.fits.gz")
    column_mass_density_h = (hdu[0].data[0, :, :] *
                             units.g * units.cm**-2)   # g cm^-2, from the star
    column_density_h = column_mass_density_h / \
        (2.0 * constants.u.to('g'))  # cm-2
    column_mass_density_v = (hdu[0].data[1, :, :] *
                             units.g * units.cm**-2)   # g cm^-2, from the disk surface
    column_density_v = column_mass_density_v / \
        (2.0 * constants.u.to('g'))  # cm-2
    hdu.close()
    hdu = fits.open(datadir + "/grid.fits.gz")
    radius_au = hdu[0].data  # au
    hdu.close()
    zeta_xray = Lx29 / radius_au**2.2 * \
        (zeta1 * (np.exp(-(column_density_v/N1)**alpha)) +
         zeta2 * (np.exp(-(column_density_h/N2)**beta)))  # s-1 per H
    zeta_xray /= 2                                               # s-1 per H2

    # Cosmic-ray ionization (Bai & Goodman 2009) ##
    zeta0 = 1.3e-17 * units.s**-1  # "standard" value
    zeta_cr = zeta0 * \
        (np.exp(-(column_mass_density_v / (96 * units.g * units.cm**-2))))

    class Model:
        pass
    model = Model()
    model.gas_num_den = gas_number_density
    model.temperature = temperature
    model.chi = chi
    model.dust_gas_m_ratio = dust_gas_mass_ratio
    model.effective_grain_size = average_grain_size
    model.zeta_x = zeta_xray
    model.cosmic = zeta_cr
    model.density = gas_mass_density
    model.grain_abundance = dust_number_density

    return (model)