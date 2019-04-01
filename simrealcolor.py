
import sys
import os
import math
import numpy
import logging
import time
import galsim
import matplotlib.pyplot as plt

def main(argv):

    cat_file_name = 'real_galaxy_catalog_example.fits'
    dir = 'data'
    # Make output directory if not already present.
    if not os.path.isdir('output'):
        os.mkdir('output')
    cube_file_name = os.path.join('output','cube_real.fits')
    psf_file_name = os.path.join('output','psf_real.fits')

    random_seed = 1512413
    sky_level = 1.e6        # ADU / arcsec^2
    pixel_scale = 0.16      # arcsec
    gal_flux = 1.e5         # arbitrary choice, makes nice (not too) noisy images
    psf_inner_fwhm = 0.1    # arcsec
    psf_outer_fwhm = 0.6    # arcsec
    psf_inner_fraction = 0.8  # fraction of total PSF flux in the inner Gaussian
    psf_outer_fraction = 0.2  # fraction of total PSF flux in the inner Gaussian
    redshift = 0.05
    ngal = 30  
    
    # SED ----------------------------------------------------- 
    path,filename = os.path.split(__file__)
    datapath      = os.path.abspath(os.path.join(path,"data/"))
    outpath       = os.path.abspath(os.path.join(path,"output/"))
    random_seed   = 1234567
    rng           = galsim.BaseDeviate(random_seed)

    # read in SEDs
    SED_names     = ['CWW_E_ext','CWW_Sbc_ext','CWW_Scd_ext','CWW_Im_ext']
    SEDs          = {}
    for SED_name in SED_names:
        SED_filename   = os.path.join(datapath,'{0}.sed'.format(SED_name))
        SED            = galsim.SED(SED_filename,wave_type='Ang')
        SEDs[SED_name] = SED.withFluxDensity(target_flux_density=1.0,wavelength=500)

    filter_names = 'gri'
    filters      = {}

    for filter_name in filter_names:
        filter_filename      = os.path.join(datapath,'LSST_{0}.dat'.format(filter_name))
        filters[filter_name] = galsim.Bandpass(filter_filename)
        filters[filter_name] = filters[filter_name].thin(rel_err=1e-4)

    #--------------------------------------------------------------------
    real_galaxy_catalog = galsim.RealGalaxyCatalog(cat_file_name, dir=dir)

    psf = galsim.Gaussian(fwhm = psf_inner_fwhm, flux = psf_inner_fraction)
    # Draw the PSF with no noise.
    psf_image = psf.drawImage(scale = pixel_scale)
    # write to file
    psf_image.write(psf_file_name)
    # Build the images
    for k in range(ngal):

        # Initialize the random number generator we will be using.
        rng = galsim.UniformDeviate(random_seed+k)

        gal = galsim.RealGalaxy(real_galaxy_catalog, index = k)

        # Set the flux
        #gal = gal.withFlux(gal_flux)
        SED       = SEDs['CWW_E_ext'].atRedshift(redshift)
        #SED       = SEDs['CWW_Sbc_ext'].atRedshift(redshift)
	galcol= galsim.Chromatic(gal,SED)
        final = galsim.Convolve([galcol, psf])

        dx = rng() - 0.5
        dy = rng() - 0.5

        # Draw the profile
	galcol  = galsim.Chromatic(final,SED)
	for filter_name,filter_ in filters.iteritems():
	    img   = galsim.ImageF(78,78,scale=pixel_scale)
            final.drawImage(filter_,image=img)
	    out_filename = os.path.join(outpath,'realcolor_{0}_{1}.fits'.format(filter_name,str(k)))
            galsim.fits.write(img,out_filename)



if __name__ == "__main__":
    main(sys.argv)
