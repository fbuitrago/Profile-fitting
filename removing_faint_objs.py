import numpy as np

def removing_faint_objs(sex_output,target,filter_objs,threshold):

	if sex_output.ndim > 1:
		mags = sex_output[3, :]
	else:
		mags = np.array([sex_output[3]])

	new_filter = mags < (mags[target]+threshold)

	output = filter_objs & new_filter

	return(output)


