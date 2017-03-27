"""
Program for fitting Sersic functions using GALFIT
By Fernando Buitrago (fbuitrago@gmail.com)
"""

from astropy.io import fits
import numpy as np
import os
import pdb
import my_python_library
from astropy.table import Table
#===========================
from run_sextractor import run_sextractor
from read_sextractor_output import read_sextractor_output
from detecting_central_obj import detecting_central_obj
from removing_far_away_objs import removing_far_away_objs
from create_mask import create_mask
from create_central_mask import create_central_mask
from create_constraints_file import create_constraints_file
from removing_faint_objs import removing_faint_objs
from write_galfit_script import write_galfit_script
from is_sex_file_empty import is_sex_file_empty
from config_different_fields import config_different_fields


#definitions=======================
SExtractor_config_file = "wfc3_H.sex"
zeropoint = 25.96
single_or_double = 1 #2 #single or double component fit
bulge_disk_decomp = False
flag_mask_central = True #if False we skip the creation of central mask (a mask containing only the central obj, good for overcrowded fields)
enlarge_radius_close_objs = 6.
enlarge_masks = 3.
threshold_for_faint = 3. #mag greater than the target galaxy
pixel_scale = 0.06 #[arcsec/pix]
exptime = 1. #it will be override by the EXPTIME keyword in the header (if it exists)
different_fields = True
camera = "WFC3"
#==================================

#PSF stars=========================
stars = np.array(["star_h_band_candels"],dtype=str) #don't add extension .fits
#==================================

#knowing the galaxies to run GALFIT
cwd = os.getcwd()
path, last_folder = os.path.split(cwd)
cpu_number = last_folder.replace("cpu","")
catalog_name = "./aux_files/gals_cpu"+cpu_number+".cat"
#reading the catalog
tt = Table.read(catalog_name, names=('ids','nothing'), format='ascii.commented_header')
ids = tt['ids']

#for each galaxy
for ii in range(ids.size):
    galaxy = ids[ii]
    
    #reading the galaxy
    img = fits.open("../../galaxy_images/"+galaxy+".fits")
    header = img[0].header
    x_size  = header['NAXIS1']
    y_size  = header['NAXIS2']
    
    #if using different fields
    if different_fields == True:
        zeropoint, pixel_scale, exptime, stars = config_different_fields(galaxy, camera, header)

    #generating subfolders for each galaxy
    if not os.path.exists(galaxy):
        os.mkdir(galaxy)

    run_sextractor(galaxy,SExtractor_config_file,zeropoint)

    flag_empty = is_sex_file_empty("./"+galaxy+"/"+galaxy+"_ex.txt")
    if flag_empty == False:
        pass
    else:
        continue 

    sex_output = read_sextractor_output(galaxy)

    target = detecting_central_obj(sex_output,x_size/2.,y_size/2.)
    #detecting_central_obj_from_cat()

    filter_objs = removing_far_away_objs(sex_output,target,enlarge_radius_close_objs)

    filter_objs = removing_faint_objs(sex_output,target,filter_objs,threshold_for_faint)

    #remove_five_sigma()

    create_mask        (galaxy,sex_output,target,enlarge_masks,filter_objs)
    masks = []
    if flag_mask_central == True:
        create_central_mask(galaxy,sex_output,target,enlarge_masks)
        masks = ["normal_mask","central_mask"]
    else:
        masks = ["normal_mask"]
    for mask in masks: create_constraints_file(galaxy,single_or_double,target,filter_objs,mask) #1 or 2 constraints file depending on the mask I use

    for star_cont in range(len(stars)):
        for mask_cont in range(len(masks)):
            galfit_script = write_galfit_script(galaxy,stars[star_cont],masks[mask_cont],x_size,y_size,zeropoint,exptime,pixel_scale,sex_output,target,single_or_double,bulge_disk_decomp,filter_objs)
            os.chdir(cwd+"/"+galaxy+"/") #changing the working directory for galfit to store there all its output
            os.system("galfit "+galfit_script)
            os.chdir(cwd)
            
