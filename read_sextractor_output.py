import numpy as np
import my_python_library
import pdb

def read_sextractor_output(galaxy):
	path_to_galaxy = "./"+galaxy+"/"

	objs,ra,dec,mag,magerr,x_image,y_image,flux_radius,ellipticity,theta,kron,a_image = \
		np.loadtxt(path_to_galaxy+galaxy+"_ex.txt", dtype = float, unpack = True)
		
	sex_output = np.array([objs,ra,dec,mag,magerr,x_image,y_image,flux_radius,ellipticity,theta,kron,a_image],float)

	return(sex_output)
