def calculating_radius_to_fit(sex_semimajor_ax,enlarge):

	#at least semimajor axis = 1 pix
	if sex_semimajor_ax <= 0.:
		sex_semimajor_ax = 1.

	return(enlarge*sex_semimajor_ax)