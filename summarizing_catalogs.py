import numpy as np
from astropy.table import Table
import os

no_cpus = 3
single_or_double = 1  #single (1) or double (2) component fit
camera = ""
name_output_cat = "gama_viking_z_single_sersic_fits.cat"
path_for_objects = "/single_sersic_fits/cpu"

#defining vectors for the final results
vector_gal_name          = np.array([])
vector_zz                = np.array([])
vector_index             = np.array([])
vector_flag_obj_detected = np.array([])
if single_or_double == 1:
    vector_mag       = np.array([])
    vector_re        = np.array([])
    vector_re_kpc    = np.array([])
    vector_nn        = np.array([])
    vector_ar        = np.array([])
    vector_pa        = np.array([])
    vector_magerr    = np.array([])
    vector_reerr     = np.array([])
    vector_reerr_kpc = np.array([])
    vector_nnerr     = np.array([])
    vector_arerr     = np.array([])
    vector_paerr     = np.array([])
    vector_mag_exist    = np.array([])
    vector_re_exist     = np.array([])
    vector_nn_exist     = np.array([])
    vector_ar_exist     = np.array([])
    vector_pa_exist     = np.array([])
else:
    vector_mag1        = np.array([])
    vector_re1         = np.array([])
    vector_re1_kpc     = np.array([])
    vector_nn1         = np.array([])
    vector_ar1         = np.array([])
    vector_pa1         = np.array([])
    vector_magerr1     = np.array([])
    vector_reerr1      = np.array([])
    vector_reerr1_kpc  = np.array([])
    vector_nnerr1      = np.array([])
    vector_arerr1      = np.array([])
    vector_paerr1      = np.array([])
    vector_mag2        = np.array([])
    vector_re2         = np.array([])
    vector_re2_kpc     = np.array([])
    vector_nn2         = np.array([])
    vector_ar2         = np.array([])
    vector_pa2         = np.array([])
    vector_magerr2     = np.array([])
    vector_reerr2      = np.array([])
    vector_reerr2_kpc  = np.array([])
    vector_nnerr2      = np.array([])
    vector_arerr2      = np.array([])
    vector_paerr2      = np.array([])
    vector_mag1_exist   = np.array([])
    vector_re1_exist    = np.array([])
    vector_nn1_exist    = np.array([])
    vector_ar1_exist    = np.array([])
    vector_pa1_exist    = np.array([])
    vector_mag2_exist   = np.array([])
    vector_re2_exist    = np.array([])
    vector_nn2_exist    = np.array([])
    vector_ar2_exist    = np.array([])
    vector_pa2_exist    = np.array([])

cwd = os.getcwd()

for ii in range(no_cpus):
    folder_for_objects = cwd + path_for_objects + str(ii+1)

    if single_or_double == 1:

        tt = Table.read(folder_for_objects+"/results_cpu"+str(ii+1)+".cat",\
             names = ["gal_id","zz","mag","mag_error","mag_exist","re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err","n","n_error","n_exist",\
                      "ar","ar_error","ar_exist","pa","pa_error","pa_exist","index","flag_obj_detected"],\
             format="ascii.commented_header")

        vector_gal_name          = np.concatenate((vector_gal_name,         np.array(tt["gal_id"])))
        vector_zz                = np.concatenate((vector_zz,               np.array(tt["zz"])))
        vector_index             = np.concatenate((vector_index,            np.array(tt["index"])))
        vector_flag_obj_detected = np.concatenate((vector_flag_obj_detected,np.array(tt["flag_obj_detected"])))
        vector_mag               = np.concatenate((vector_mag,              np.array(tt["mag"])))
        vector_re                = np.concatenate((vector_re,               np.array(tt["re_pix"])))
        vector_re_kpc            = np.concatenate((vector_re_kpc,           np.array(tt["re_kpc"])))
        vector_nn                = np.concatenate((vector_nn,               np.array(tt["n"])))
        vector_ar                = np.concatenate((vector_ar,               np.array(tt["ar"])))
        vector_pa                = np.concatenate((vector_pa,               np.array(tt["pa"])))
        vector_magerr            = np.concatenate((vector_magerr,           np.array(tt["mag_error"])))
        vector_reerr             = np.concatenate((vector_reerr,            np.array(tt["re_pix_error"])))
        vector_reerr_kpc         = np.concatenate((vector_reerr_kpc,        np.array(tt["re_kpc_err"])))
        vector_nnerr             = np.concatenate((vector_nnerr,            np.array(tt["n_error"])))
        vector_arerr             = np.concatenate((vector_arerr,            np.array(tt["ar_error"])))
        vector_paerr             = np.concatenate((vector_paerr,            np.array(tt["pa_error"])))
        vector_mag_exist         = np.concatenate((vector_mag_exist,        np.array(tt["mag_exist"])))
        vector_re_exist          = np.concatenate((vector_re_exist,         np.array(tt["re_exist"])))
        vector_nn_exist          = np.concatenate((vector_nn_exist,         np.array(tt["n_exist" ])))
        vector_ar_exist          = np.concatenate((vector_ar_exist,         np.array(tt["ar_exist"])))
        vector_pa_exist          = np.concatenate((vector_pa_exist,         np.array(tt["pa_exist"])))

    else:
        tt = Table.read(folder_for_objects+"/results_cpu"+str(ii+1)+".cat",\
             names = ["gal_id","zz","mag1","mag_error1","mag1_exist","re_pix1",\
                      "re_pix_error1","re_kpc1","re1_exist","re_kpc_error1","n1","n_error1","n1_exist",\
                      "ar1","ar_error1","ar1_exist","pa1","pa_error1","pa1_exist","mag2",\
                      "mag_error2","mag2_exist","re_pix2","re_pix_error2","re2_exist","re_kpc2","re_kpc_error2",\
                      "n2","n_error2","n2_exist","ar2","ar_error2","ar2_exist","pa2",\
                      "pa_error2","pa2_exist","index","flag_obj_detected"],\
             format="ascii.commented_header")

        vector_gal_name           = np.concatenate((vector_gal_name,          np.array(tt["gal_id"])))
        vector_zz                 = np.concatenate((vector_zz,                np.array(tt["zz"])))
        vector_index             = np.concatenate((vector_index,              np.array(tt["index"])))
        vector_flag_obj_detected  = np.concatenate((vector_flag_obj_detected, np.array(tt["flag_obj_detected"])))
        vector_mag1               = np.concatenate((vector_mag1,              np.array(tt["mag1"])))
        vector_re1                = np.concatenate((vector_re1,               np.array(tt["re_pix1"])))
        vector_re1_kpc            = np.concatenate((vector_re1_kpc,           np.array(tt["re_kpc1"])))
        vector_nn1                = np.concatenate((vector_nn1,               np.array(tt["n1"])))
        vector_ar1                = np.concatenate((vector_ar1,               np.array(tt["ar1"])))
        vector_pa1                = np.concatenate((vector_pa1,               np.array(tt["pa1"])))
        vector_magerr1            = np.concatenate((vector_magerr1,           np.array(tt["mag_error1"])))
        vector_reerr1             = np.concatenate((vector_reerr1,            np.array(tt["re_pix_error1"])))
        vector_reerr1_kpc         = np.concatenate((vector_reerr1_kpc,        np.array(tt["re_kpc_error1"])))
        vector_nnerr1             = np.concatenate((vector_nnerr1,            np.array(tt["n_error1"])))
        vector_arerr1             = np.concatenate((vector_arerr1,            np.array(tt["ar_error1"])))
        vector_paerr1             = np.concatenate((vector_paerr1,            np.array(tt["pa_error1"])))
        vector_mag2               = np.concatenate((vector_mag2,              np.array(tt["mag2"])))
        vector_re2                = np.concatenate((vector_re2,               np.array(tt["re_pix2"])))
        vector_re2_kpc            = np.concatenate((vector_re2_kpc,           np.array(tt["re_kpc2"])))
        vector_nn2                = np.concatenate((vector_nn2,               np.array(tt["n2"])))
        vector_ar2                = np.concatenate((vector_ar2,               np.array(tt["ar2"])))
        vector_pa2                = np.concatenate((vector_pa2,               np.array(tt["pa2"])))
        vector_magerr2            = np.concatenate((vector_magerr2,           np.array(tt["mag_error2"])))
        vector_reerr2             = np.concatenate((vector_reerr2,            np.array(tt["re_pix_error2"])))
        vector_reerr2_kpc         = np.concatenate((vector_reerr2_kpc,        np.array(tt["re_kpc_error2"])))
        vector_nnerr2             = np.concatenate((vector_nnerr2,            np.array(tt["n_error2"])))
        vector_arerr2             = np.concatenate((vector_arerr2,            np.array(tt["ar_error2"])))
        vector_paerr2             = np.concatenate((vector_paerr2,            np.array(tt["pa_error2"])))
        vector_mag1_exist         = np.concatenate((vector_mag1_exist,        np.array(tt["mag1_exist"])))
        vector_re1_exist          = np.concatenate((vector_re1_exist,         np.array(tt["re1_exist"])))
        vector_nn1_exist          = np.concatenate((vector_nn1_exist,         np.array(tt["n1_exist" ])))
        vector_ar1_exist          = np.concatenate((vector_ar1_exist,         np.array(tt["ar1_exist"])))
        vector_pa1_exist          = np.concatenate((vector_pa1_exist,         np.array(tt["pa1_exist"])))
        vector_mag2_exist         = np.concatenate((vector_mag2_exist,        np.array(tt["mag2_exist"])))
        vector_re2_exist          = np.concatenate((vector_re2_exist,         np.array(tt["re2_exist"])))
        vector_nn2_exist          = np.concatenate((vector_nn2_exist,         np.array(tt["n2_exist" ])))
        vector_ar2_exist          = np.concatenate((vector_ar2_exist,         np.array(tt["ar2_exist"])))
        vector_pa2_exist          = np.concatenate((vector_pa2_exist,         np.array(tt["pa2_exist"])))
        

#writing the final catalogs
if single_or_double == 1:
    data = Table( [vector_gal_name,vector_zz,\
                   vector_mag,vector_magerr,vector_mag_exist,\
                   vector_re, vector_reerr, vector_re_exist,vector_re_kpc,vector_reerr_kpc,\
                   vector_nn, vector_nnerr, vector_nn_exist,\
                   vector_ar, vector_arerr, vector_ar_exist,\
                   vector_pa, vector_paerr, vector_pa_exist,\
                   vector_index,vector_flag_obj_detected],\
                   names=("gal_id","zz","mag","mag_error","mag_exist","re_pix","re_pix_error","re_exist","re_kpc","re_kpc_err","n","n_error","n_exist",\
                         "ar","ar_error","ar_exist","pa","pa_error","pa_exist","index_of_best_fit","flag_obj_detected") )
else:
    data = Table( [vector_gal_name,vector_zz,\
                   vector_mag1,vector_magerr1,vector_mag1_exist,\
                   vector_re1, vector_reerr1, vector_re1_exist,vector_re1_kpc,vector_reerr1_kpc,\
                   vector_nn1, vector_nnerr1, vector_nn1_exist,\
                   vector_ar1, vector_arerr1, vector_ar1_exist,\
                   vector_pa1, vector_paerr1, vector_pa1_exist,\
                   vector_mag2,vector_magerr2,vector_mag2_exist,\
                   vector_re2, vector_reerr2, vector_re2_exist,vector_re2_kpc,vector_reerr2_kpc,\
                   vector_nn2, vector_nnerr2, vector_nn2_exist,
                   vector_ar2, vector_arerr2, vector_ar2_exist,
                   vector_pa2, vector_paerr2, vector_pa2_exist,\
                   vector_index,vector_flag_obj_detected], \
                   names=("gal_id","zz","mag1","mag_error1","mag1_exist","re_pix1",\
                          "re_pix_error1","re1_exist","re_kpc1","re_kpc_error1","n1","n_error1","n1_exist",\
                          "ar1","ar_error1","ar1_exist","pa1","pa_error1","pa1_exist","mag2",\
                          "mag_error2","mag2_exist","re_pix2","re_pix_error2","re2_exist","re_kpc2","re_kpc_error2",\
                          "n2","n_error2","n2_exist","ar2","ar_error2","ar2_exist","pa2",\
                          "pa_error2","pa2_exist","index_of_best_fit","flag_obj_detected") )
data.write(name_output_cat, format='ascii.commented_header')
