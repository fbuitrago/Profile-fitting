"""
Program for cutting stamps for individual galaxies from a fits mosaic
Input
- Science image
- Weight image already converted into rms image
Output:
- The galaxy stamp themselves
- A text file containing the positions in pixel of those stamps in the total mosaic
By Fernando Buitrago (fbuitrago@gmail.com)
"""

#import montage_wrapper as montage
import numpy as np
import os
import pdb
from astropy.io import fits
from astropy import wcs
from astropy.io import ascii
import montage_wrapper as montage


def ids_as_str(ids):
    #it returns a numpy array of strings, when receiving the ids as float numbers
    ids = np.array(ids, dtype=np.str)
    ids = np.char.rstrip(ids, '.0')
    return ids


which_sample = "low_z" #"high_z"
path_img = ["/mnt/disk1/fb/","/mnt/disk1/fb/", #GOODS-S
            "/mnt/disk1/fb/","/mnt/disk1/fb/", #GOODS-N
            "/mnt/disk2/fb/candels/cosmos/","/mnt/disk1/fb/", #COSMOS
            "/mnt/disk2/fb/candels/before_jun_2013/uds/uds_dr10/","/mnt/disk1/fb/", #UDS
            "/mnt/disk1/fb/","/mnt/disk1/fb/"] #EGS
path_rms = ["/mnt/disk1/fb/","/mnt/disk1/fb/", #GOODS-S
            "/mnt/disk1/fb/","/mnt/disk1/fb/", #GOODS-N
            "/mnt/disk2/fb/candels/cosmos/","/mnt/disk1/fb/", #COSMOS
            "/mnt/disk2/fb/candels/before_jun_2013/uds/uds_dr10/","/mnt/disk1/fb/", #UDS
            "/mnt/disk1/fb/","/mnt/disk1/fb/"] #EGS
path_cat = "/mnt/disk1/fb/massive_disks/"
name_img = ["hlsp_hlf_hst_wfc3-60mas_goodss_f160w_v1.5_sci.fits","hlsp_hlf_hst_acs-30mas_goodss_f814w_v1.5_sci.fits", #GOODS-S
            "hlsp_candels_hst_wfc3_gn-tot-60mas_f160w_v1.0_drz.fits","goodsn_all_acs_wfc_f814w_030mas_v2.0_drz.fits", #GOODS-N
            "hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits","hlsp_candels_hst_acs_cos-tot_f814w_v1.0_drz.fits", #COSMOS
            "hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits","hlsp_candels_hst_acs_uds-tot_f814w_v1.0_drz.fits", #UDS
            "hlsp_candels_hst_wfc3_egs-tot-60mas_f160w_v1.0_drz.fits","egs_all_acs_wfc_f814w_030mas_v1.1_drz.fits"] #EGS
name_rms = ["hlsp_hlf_hst_wfc3-60mas_goodss_f160w_v1.5_wht.fits","hlsp_hlf_hst_acs-30mas_goodss_f814w_v1.5_wht.fits", #GOODS-S
            "hlsp_candels_hst_wfc3_gn-tot-60mas_f160w_v1.0_rms.fits","goodsn_all_acs_wfc_f814w_030mas_v2.0_rms.fits", #GOODS-N
            "hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_rms.fits","hlsp_candels_hst_acs_cos-tot_f814w_v1.0_rms.fits", #COSMOS
            "hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_wht.fits","hlsp_candels_hst_acs_uds-tot_f814w_v1.0_wht.fits", #UDS
            "hlsp_candels_hst_wfc3_egs-tot-60mas_f160w_v1.0_rms.fits","egs_all_acs_wfc_f814w_030mas_v1.1_rms.fits"] #EGS


name_cat = "massive_"+which_sample+"_bin.cat"

#reading the catalog
cat = ascii.read(path_cat+name_cat)
gal    = np.array(cat[0][:],dtype = str)
ra     = np.array(cat[1][:],dtype = float)
dec    = np.array(cat[2][:],dtype = float)
region = np.array(cat[3][:],dtype = str)
#--------------------------------
#gal,ra,dec = np.loadtxt(path_cat+name_cat, dtype = float, unpack = True) #OJO: the names are also floats!
#gal = ids_as_str(gal)
#--------------------------------
#table=numpy.genfromtxt("massive_low_z_bin_uds.cat", dtype='U', comments='#')

#for each mosaic to cut
for kk in range(len(name_img)):
    
    print(name_img[kk])
    
    #ACS cutouts are 2x bigger than WFC3 ones
    if name_img[kk].find("30mas") != -1 or name_img[kk].find("acs") != -1:
        side = 399 #[pix]
        filt = "i_band"
    else:
        side = 199 #[pix]
        filt = "h_band"
    
    #detecting the mosaic's patch of the sky
    if name_img[kk].find("goodss") != -1:
        field = "goodss"
    elif name_img[kk].find("gn") != -1 or name_img[kk].find("goodsn") != -1:
        field = "goodsn"
    elif name_img[kk].find("cos") != -1:
        field = "cosmos"
    elif name_img[kk].find("uds") != -1:
        field = "uds"
    elif name_img[kk].find("egs") != -1:
        field = "aegis"
        
    #obtaining the galaxy coordinates in pixels
    img = fits.open(path_img[kk]+name_img[kk])
    ww = wcs.WCS(img[0].header)
    coo = []
    indices = [] #to know which galaxies from the whole gal vector
    #selecting galaxies that belong to our field only
    for ii in range(len(gal)):
        if region[ii] == field:
            coo.append([ra[ii],dec[ii]])
            indices.append(ii) 
    coo = np.array(coo, dtype = float)
    pix_coo = ww.wcs_world2pix(coo,1)
    indices = np.array(indices, dtype = int)
    
    print(side,filt,field,gal[indices])
    if indices.size == 0: continue

    #creating folders for the stamps and the rms stamps
    new_galaxy_stamps_folder = "./galaxy_images"+"_"+filt+"_"+field+"/"
    if not os.path.exists(new_galaxy_stamps_folder):
        os.mkdir(new_galaxy_stamps_folder)
    new_sigma_stamps_folder = "./sigma_images"+"_"+filt+"_"+field+"/"
    if not os.path.exists(new_sigma_stamps_folder):
        os.mkdir(new_sigma_stamps_folder)
        
    #calculating the coordinates in case I need to paste the galaxies back in the image
    x0 = pix_coo[:,0] - (side/2.)
    x1 = pix_coo[:,0] + (side/2.)
    y0 = pix_coo[:,1] - (side/2.)
    y1 = pix_coo[:,1] + (side/2.)

    #Round elements of the array to the nearest integer
    x0 = np.rint(x0); x0 = x0.astype(int)
    x1 = np.rint(x1); x1 = x1.astype(int)
    y0 = np.rint(y0); y0 = y0.astype(int)
    y1 = np.rint(y1); y1 = y1.astype(int)

    #creating the galaxy stamps
    for ii in range(len(indices)):
        montage.mSubimage_pix(path_img[kk]+name_img[kk],new_galaxy_stamps_folder+gal[indices[ii]]+".fits"    , x0[ii], y0[ii], side)
        montage.mSubimage_pix(path_rms[kk]+name_rms[kk],new_sigma_stamps_folder +gal[indices[ii]]+"_rms.fits", x0[ii], y0[ii], side)

        #identifying wrong pixels in the rms image (either NaN or Inf; also the ones equal to 0, i.e. beyond borders)
        new_img = fits.open(new_sigma_stamps_folder+gal[indices[ii]]+"_rms.fits")
        new_rms = np.array(new_img[0].data)
        #---
        bad_rms_pixels = np.isfinite(new_rms)
        mask_bad_rms_pixels = np.logical_not(bad_rms_pixels)
        #---
        #beyond_borders_pixels = new_rms == 0.
        #mask_beyond_borders = np.logical_not(beyond_borders_pixels)
        #---
        mask_bad_rms_pixels = mask_bad_rms_pixels #or mask_beyond_borders
        
        #in case we are managing the weight image instead of the rms image
        if name_rms[kk].find("wht") != -1:
            new_rms = 1./np.sqrt(new_rms)
        else:
            new_rms = new_rms

        #---repeating this part while I cannot flag the pixels with 0 value before
        bad_rms_pixels = np.isfinite(new_rms)
        mask_bad_rms_pixels = np.logical_not(bad_rms_pixels)
        #---

        #writing the rms image
        new_rms[mask_bad_rms_pixels] = 99999999. #flagging bad pixels
        fits.writeto(       new_sigma_stamps_folder+gal[indices[ii]]+"_rms.fits",new_rms,new_img[0].header,clobber=True)
            
    #saving the coordinates where I cut the galaxies from
    ascii.write([gal[indices],x0,x1,y0,y1],names=['id','x0','x1','y0','y1'],output=filt+"_"+field+"_"+which_sample+".coo",Writer = ascii.CommentedHeader)
