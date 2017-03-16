import pdb

def removing_faint_objs(sex_output,target,filter_objs,threshold):

	mags = sex_output[3, :]

	new_filter = mags < (mags[target]+threshold)

	output = filter_objs & new_filter

	return(output)


