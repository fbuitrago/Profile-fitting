"""
Program for dividing the objects from a master catalog among the different subtrees according to the number of cpus we have
By Fernando Buitrago (fbuitrago@gmail.com)
"""

import numpy as np
import os
from shutil import copyfile
from astropy.table import Table
import pdb


def ids_as_str(ids):
    #it returns a numpy array of strings, when receiving the ids as float numbers
    ids = np.array(ids, dtype=np.str)
    ids[:] = [ele[:-2] for ele in ids]
    return ids


no_cpus = 3
catalog_file = "massive_high_z_bin.cat"
base_path = "./high_z_galaxy_images/"


ids = os.listdir(base_path+"./galaxy_images/")
ids[:] = [ele[:-5] for ele in ids] #I remove the .fits extension
ids = np.array(ids, dtype = str)
ids = np.sort(ids)

#tt = Table.read("massive_low_z_bin.cat", format = "ascii.commented_header", names=('ids', 'ra', 'dec', 'field'))
#ids = ids_as_str(tt['ids'])
#ra  = tt['ra']
#dec = tt['dec']
#region = tt['field']; region = np.array(region, dtype=str)

no_galaxies = ids.size
no_gals_per_directory = int(no_galaxies/no_cpus)

start_value = 0
for ii in range(1,no_cpus+1): #+1 because range upper limit otherwise is not inclusive (it does -1)

    cpu_folder = "./cpu"+str(ii)
    cpu_cat    = cpu_folder+"/aux_files/gals_cpu"+str(ii)+".cat"

    #generating the subfolders
    if not os.path.exists(cpu_folder):
        os.mkdir(cpu_folder)
    #if not os.path.exists(cpu_folder + "/galaxy_images"):
        #os.mkdir(cpu_folder + "/galaxy_images")
    #if not os.path.exists(cpu_folder + "/sigma_images"):
        #os.mkdir(cpu_folder + "/sigma_images")
    if not os.path.exists(cpu_folder + "/aux_files"):
        os.mkdir(cpu_folder + "/aux_files")

    #sending files to folders
    final_value = no_gals_per_directory*ii
    for jj in range(start_value,final_value): #+1 because range upper limit otherwise is not inclusive (it does -1)

        if ii is not no_cpus:
            #copyfile("./galaxy_images/"+ids[jj]+".fits"    ,cpu_folder+"/galaxy_images/"+ids[jj]+".fits")
            #copyfile("./sigma_images/" +ids[jj]+"_rms.fits",cpu_folder+"/sigma_images/" +ids[jj]+"_rms.fits")
            if jj == final_value-1: #for the last element
                data = Table([ids[start_value:final_value],np.zeros(final_value-start_value)], names=['gal','nothing'])
                data.write(cpu_cat, format = "ascii.commented_header")
                start_value = final_value
        else:
            final_value = no_galaxies #because I will need to add all the restant galaxies to the final folder
            #copyfile("./galaxy_images/"+ids[jj]+".fits"    ,cpu_folder+"/galaxy_images/"+ids[jj]+".fits")
            #copyfile("./sigma_images/" +ids[jj]+"_rms.fits",cpu_folder+"/sigma_images/" +ids[jj]+"_rms.fits")
            if jj == final_value-1: #for the last element 
                data = Table([ids[start_value:final_value],np.zeros(final_value-start_value)], names=['gal','nothing'])
                data.write(cpu_cat, format = "ascii.commented_header")
                #start_value = final_value -> This is not necessary, as it is the final folder
        
