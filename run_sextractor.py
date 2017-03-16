import os

def run_sextractor(galaxy,SEx_file,zeropoint):
	path_to_galaxy = "./"+galaxy+"/"
	image =          "./galaxy_images/"+galaxy+".fits"
	sigma_image =    "./sigma_images/" +galaxy+"_rms.fits"
	SEx_catalog =    path_to_galaxy+galaxy+"_ex.txt"
	SEx_checkimage = path_to_galaxy+galaxy+"_ex.fits"
	
	order = '/usr/local/scisoft/bin/sex '+image+' -c '+'../'+SEx_file+' -CATALOG_NAME '+SEx_catalog+ \
	        ' -CHECKIMAGE_TYPE SEGMENTATION'+' -CHECKIMAGE_NAME '+ SEx_checkimage+' -MAG_ZEROPOINT '+str(zeropoint) #+ \
	        #' -WEIGHT_TYPE MAP_RMS'+' -WEIGHT_IMAGE '+sigma_image

	print(order)
	os.system(order)