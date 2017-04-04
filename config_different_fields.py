"""
Program for retrieving the data specific to every field (zeropoint, pixel scale and all the like)
By Fernando Buitrago (fbuitrago@gmail.com)
"""

import numpy as np
import pdb


def config_different_fields(galaxy,camera,header):

    if camera == 'ACS':
        pixel_scale = 0.03 #[arcsec/pix]
        #detecting the mosaic's patch of the sky
        if galaxy.find("goodss") != -1:
            zeropoint = 25.94
            exptime = 840211.
            stars = np.array(["tiny_tim_iband","goodss_hlf_star_i_band_fer"],dtype=str) #don't add extension .fits
        elif galaxy.find("goodsn") != -1:
            zeropoint = 25.947 #CANDELS website
            exptime = 4400. #CANDELS website
            stars = np.array(["tiny_tim_iband","goodsn_star_i_band_fer"],dtype=str) #don't add extension .fits
        elif galaxy.find("cos") != -1:
            zeropoint = 25.947
            exptime = 6900.
            stars = np.array(["tiny_tim_iband","cosmos_star_i_band_fer"],dtype=str) #don't add extension .fits
        elif galaxy.find("uds") != -1:
            zeropoint = 25.94333
            exptime = 5700.
            stars = np.array(["tiny_tim_iband","uds_star_i_band_fer"],dtype=str) #don't add extension .fits
        elif galaxy.find("aegis") != -1:
            zeropoint = 25.94333 #CANDELS website
            exptime = 4392. #CANDELS website
            stars = np.array(["tiny_tim_iband","egs_star_i_band_fer"],dtype=str) #don't add extension .fits
    elif camera == 'WFC3':
        pixel_scale = 0.06 #[arcsec/pix]
        #detecting the mosaic's patch of the sky
        if galaxy.find("goodss") != -1:
            zeropoint = 25.94
            exptime = header["EXPTIME"] #exptime = 596277.
            stars = np.array(["gs_deep_f160w_v0.5_psf"],dtype=str) #don't add extension .fits
        elif galaxy.find("goodsn") != -1:
            zeropoint = 25.9463 #CANDELS website
            exptime = header["EXPTIME"] #exptime = 1200. #CANDELS website
            stars = np.array(["gs_deep_f160w_v0.5_psf"],dtype=str) #don't add extension .fits
        elif galaxy.find("cos") != -1:
            zeropoint = 25.9463
            exptime = header["EXPTIME"] #exptime = 3200.
            stars = np.array(["gs_deep_f160w_v0.5_psf"],dtype=str) #don't add extension .fits
        elif galaxy.find("uds") != -1:
            zeropoint = 25.96
            exptime = header["EXPTIME"] #exptime = 3300.
            stars = np.array(["gs_deep_f160w_v0.5_psf"],dtype=str) #don't add extension .fits
        elif galaxy.find("aegis") != -1:
            zeropoint = 25.96 #CANDELS website
            exptime = header["EXPTIME"] #exptime = 1600. #CANDELS website
            stars = np.array(["gs_deep_f160w_v0.5_psf"],dtype=str) #don't add extension .fits
     
    return(zeropoint, pixel_scale, exptime, stars)
