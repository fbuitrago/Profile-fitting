import numpy as np
from calculating_radius_to_fit import calculating_radius_to_fit
import pdb

def removing_far_away_objs(sex_output,target,enlarge_radius_close_objs):
	x_sextractor = sex_output[5, :]
	y_sextractor = sex_output[6, :]
	a_image      = sex_output[11,:]

	xc = x_sextractor[target]
	yc = y_sextractor[target]

	#detecting central obj
	dists = np.sqrt( ((xc-x_sextractor)**2.) + ((yc-y_sextractor)**2.) )

	#assuming only objects within a certain radius from the primary target
	radfit = calculating_radius_to_fit(a_image[target],enlarge_radius_close_objs)
	filter_dists = dists <= radfit

	#only_good = np.delete(sex_output, XXX, axis=1)
	return(filter_dists)
