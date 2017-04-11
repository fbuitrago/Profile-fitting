from astropy.io import fits
import my_python_library
import numpy as np
import os
import pdb
import string
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
#===========================
from config_different_fields_no_headers import config_different_fields_no_headers


def detect_bad(param):
    flag_special=0
    flag_fixed=0

    if string.find(param,"*") != -1:
        flag_special=1
        param=string.split(param,"*")
        param=param[1]
    if string.find(param,"[") != -1:
        flag_fixed=1
        param=string.split(param,"[")
        param=param[1]
        param=string.split(param,"]")
        param=param[0]

    param=string.split(param)[0]    
    param=float(param)

    return(param,flag_special,flag_fixed)


#ASSUMPTION
#All the even analyses are done with the normal masks

#CONSTANTS
different_fields = True
camera = "WFC3"
flag_mask_central = True ##if False we skip the creation of central mask (a mask containing only the central obj, good for overcrowded fields)
single_or_double = 1  #single (1) or double (2) component fit
ini_conds = ""#grid initial conditions "_1","_2",...

#moving distances from kpc to arcsec
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

#knowing the galaxies to run GALFIT
cwd = os.getcwd()
path, last_folder = os.path.split(cwd)
cpu_number = last_folder.replace("cpu","")
catalog_name = "./aux_files/gals_cpu"+cpu_number+".cat"
#reading the catalog
tt = Table.read(catalog_name, names=('ids','zz'), format='ascii.commented_header')
galaxies = tt['ids']
zz       = tt['zz']

#defining the masks
masks = []
if flag_mask_central == True:
    masks = ["normal_mask","central_mask"]
else:
    masks = ["normal_mask"]
masks = np.array(masks,dtype=str)

#defining vectors for the final results
vector_gal_name          = np.array([])
vector_zz                = np.array([])
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

#reading the headers
gal_list = []
for ii in range(len(galaxies)):
    #I start with an empty dictionary
    if single_or_double == 1:
        gal = {'id': -99, 'zz': [], 'chi2': [], 'chi2nu': [], 'xc': [], 'yc': [], 'mag': [], 'mag_flag_special': [], 'mag_flag_fixed': [],
               're': [], 're_flag_special': [], 're_flag_fixed': [], 'nn': [], 'nn_flag_special': [], 'nn_flag_fixed': [], 
               'ar': [], 'ar_flag_special': [], 'ar_flag_fixed': [], 'pa': [], 'pa_flag_special': [], 'pa_flag_fixed': []}
    else:
        gal = {'id': -99, 'zz': [], 'chi2': [], 'chi2nu': [], 'xc1': [], 'yc1': [], 'mag1': [], 'mag_flag_special1': [], 'mag_flag_fixed1': [],
               're1': [], 're_flag_special1': [], 're_flag_fixed1': [], 'nn1': [], 'nn_flag_special1': [], 'nn_flag_fixed1': [], 
               'ar1': [], 'ar_flag_special1': [], 'ar_flag_fixed1': [], 'pa1': [], 'pa_flag_special1': [], 'pa_flag_fixed1': [],
                                      'xc2': [], 'yc2': [], 'mag2': [], 'mag_flag_special2': [], 'mag_flag_fixed2': [],
               're2': [], 're_flag_special2': [], 're_flag_fixed2': [], 'nn2': [], 'nn_flag_special2': [], 'nn_flag_fixed2': [], 
               'ar2': [], 'ar_flag_special2': [], 'ar_flag_fixed2': [], 'pa2': [], 'pa_flag_special2': [], 'pa_flag_fixed2': []}
               
    gal['id'] = galaxies[ii]
    gal['zz'] = zz[ii]
    #if using different fields
    if different_fields == True:
        zeropoint, pix_scale, stars = config_different_fields_no_headers(gal['id'], camera)
    
    for jj in range(len(stars)):
        for kk in range(len(masks)):
            #for cont in range(len(ini_conds)):
            filename = "./"+galaxies[ii]+"/"+galaxies[ii]+"_"+stars[jj]+"_"+masks[kk]+".fits"
            print(filename)
            if os.path.exists(filename) == True:
                data = fits.open(filename)
                hdr = data[2].header

                gal['chi2'].append(float(hdr['CHISQ']))
                gal['chi2nu'].append(float(hdr["CHI2NU"]))
                if single_or_double == 1:
                    gal["xc"].append(hdr["2_XC"].split()[0])
                    gal["yc"].append(hdr["2_YC"].split()[0])
                    param,fs,fx = detect_bad(hdr["2_MAG"]); gal["mag"].append(float(param)); gal["mag_flag_special"].append(float(fs)); gal["mag_flag_fixed"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["2_RE"]);  gal["re"].append(float(param));  gal["re_flag_special"].append(float(fs));  gal["re_flag_fixed"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["2_N"]);   gal["nn"].append(float(param));  gal["nn_flag_special"].append(float(fs));  gal["nn_flag_fixed"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["2_AR"]);  gal["ar"].append(float(param));  gal["ar_flag_special"].append(float(fs));  gal["ar_flag_fixed"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["2_PA"]);  gal["pa"].append(float(param));  gal["pa_flag_special"].append(float(fs));  gal["pa_flag_fixed"].append(float(fx))
                else:
                    gal["xc1"].append(float(hdr["2_XC"].split()[0]))
                    gal["yc1"].append(float(hdr["2_YC"].split()[0]))
                    param,fs,fx = detect_bad(hdr["2_MAG"]); gal["mag1"].append(float(param)); gal["mag_flag_special1"].append(float(fs)); gal["mag_flag_fixed1"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["2_RE"]);  gal["re1"].append(float(param));  gal["re_flag_special1"].append(float(fs));  gal["re_flag_fixed1"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["2_N"]);   gal["nn1"].append(float(param));  gal["nn_flag_special1"].append(float(fs));  gal["nn_flag_fixed1"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["2_AR"]);  gal["ar1"].append(float(param));  gal["ar_flag_special1"].append(float(fs));  gal["ar_flag_fixed1"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["2_PA"]);  gal["pa1"].append(float(param));  gal["pa_flag_special1"].append(float(fs));  gal["pa_flag_fixed1"].append(float(fx))
                    gal["xc2"].append(float(hdr["3_XC"].split()[0]))
                    gal["yc2"].append(float(hdr["3_YC"].split()[0]))
                    param,fs,fx = detect_bad(hdr["3_MAG"]); gal["mag2"].append(float(param)); gal["mag_flag_special2"].append(float(fs)); gal["mag_flag_fixed2"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["3_RE"]);  gal["re2"].append(float(param));  gal["re_flag_special2"].append(float(fs));  gal["re_flag_fixed2"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["3_N"]);   gal["nn2"].append(float(param));  gal["nn_flag_special2"].append(float(fs));  gal["nn_flag_fixed2"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["3_AR"]);  gal["ar2"].append(float(param));  gal["ar_flag_special2"].append(float(fs));  gal["ar_flag_fixed2"].append(float(fx))
                    param,fs,fx = detect_bad(hdr["3_PA"]);  gal["pa2"].append(float(param));  gal["pa_flag_special2"].append(float(fs));  gal["pa_flag_fixed2"].append(float(fx))
            else:
                if single_or_double == 1:
                    gal["xc"] .append(float('NaN'))
                    gal["yc"] .append(float('NaN'))
                    gal["mag"].append(float('NaN')); gal["mag_flag_special"].append(float('NaN')); gal["mag_flag_fixed"].append(float('NaN'))
                    gal["re"] .append(float('NaN'));  gal["re_flag_special"].append(float('NaN'));  gal["re_flag_fixed"].append(float('NaN'))
                    gal["nn"] .append(float('NaN'));  gal["nn_flag_special"].append(float('NaN'));  gal["nn_flag_fixed"].append(float('NaN'))
                    gal["ar"] .append(float('NaN'));  gal["ar_flag_special"].append(float('NaN'));  gal["ar_flag_fixed"].append(float('NaN'))
                    gal["pa"] .append(float('NaN'));  gal["pa_flag_special"].append(float('NaN'));  gal["pa_flag_fixed"].append(float('NaN'))
                else:
                    gal["xc1"] .append(float('NaN'))
                    gal["yc1"] .append(float('NaN'))
                    gal["mag1"].append(float('NaN')); gal["mag_flag_special1"].append(float('NaN')); gal["mag_flag_fixed1"].append(float('NaN'))
                    gal["re1"] .append(float('NaN'));  gal["re_flag_special1"].append(float('NaN'));  gal["re_flag_fixed1"].append(float('NaN'))
                    gal["nn1"] .append(float('NaN'));  gal["nn_flag_special1"].append(float('NaN'));  gal["nn_flag_fixed1"].append(float('NaN'))
                    gal["ar1"] .append(float('NaN'));  gal["ar_flag_special1"].append(float('NaN'));  gal["ar_flag_fixed1"].append(float('NaN'))
                    gal["pa1"] .append(float('NaN'));  gal["pa_flag_special1"].append(float('NaN'));  gal["pa_flag_fixed1"].append(float('NaN'))
                    gal["xc2"] .append(float('NaN'))
                    gal["yc2"] .append(float('NaN'))
                    gal["mag1"].append(float('NaN')); gal["mag_flag_special2"].append(float('NaN')); gal["mag_flag_fixed2"].append(float('NaN'))
                    gal["re1"] .append(float('NaN'));  gal["re_flag_special2"].append(float('NaN'));  gal["re_flag_fixed2"].append(float('NaN'))
                    gal["nn1"] .append(float('NaN'));  gal["nn_flag_special2"].append(float('NaN'));  gal["nn_flag_fixed2"].append(float('NaN'))
                    gal["ar1"] .append(float('NaN'));  gal["ar_flag_special2"].append(float('NaN'));  gal["ar_flag_fixed2"].append(float('NaN'))
                    gal["pa1"] .append(float('NaN'));  gal["pa_flag_special2"].append(float('NaN'));  gal["pa_flag_fixed2"].append(float('NaN'))
    gal_list.append(gal)

#creating vectors to contain the information in the headers
galaxy = {}
for gal_dict in gal_list:
    galaxy.update(gal_dict)

    #converting everything into numpy vectors (_array suffix)
    zz_array     = np.array( galaxy["zz"] )
    chi2_array   = np.array( galaxy["chi2"] )
    chi2nu_array = np.array( galaxy["chi2nu"] )
    if single_or_double == 1:
        mag_flag_special_array = np.array(galaxy["mag_flag_special"])
        re_flag_special_array  = np.array(galaxy["re_flag_special"])
        nn_flag_special_array  = np.array(galaxy["nn_flag_special"])
        ar_flag_special_array  = np.array(galaxy["ar_flag_special"])
        pa_flag_special_array  = np.array(galaxy["pa_flag_special"])             
    else:
        mag_flag_special1_array = np.array(galaxy["mag_flag_special1"])
        re_flag_special1_array  = np.array(galaxy["re_flag_special1"])
        nn_flag_special1_array  = np.array(galaxy["nn_flag_special1"])
        ar_flag_special1_array  = np.array(galaxy["ar_flag_special1"])
        pa_flag_special1_array  = np.array(galaxy["pa_flag_special1"]) 
        mag_flag_special2_array = np.array(galaxy["mag_flag_special2"])
        re_flag_special2_array  = np.array(galaxy["re_flag_special2"])
        nn_flag_special2_array  = np.array(galaxy["nn_flag_special2"])
        ar_flag_special2_array  = np.array(galaxy["ar_flag_special2"])
        pa_flag_special2_array  = np.array(galaxy["pa_flag_special2"]) 

    #to see which analyses converge from the ones with the normal mask
    if single_or_double == 1:
        is_good_analysis = mag_flag_special_array[::2] + mag_flag_special_array[::2] + mag_flag_special_array[::2] + mag_flag_special_array[::2] + mag_flag_special_array[::2]
    else:
        is_good_analysis = mag_flag_special1_array[::2] + mag_flag_special1_array[::2] + mag_flag_special1_array[::2] + mag_flag_special1_array[::2] + mag_flag_special1_array[::2] + \
                           mag_flag_special2_array[::2] + mag_flag_special2_array[::2] + mag_flag_special2_array[::2] + mag_flag_special2_array[::2] + mag_flag_special2_array[::2]

    if np.any( np.logical_not(is_good_analysis) ): #if any of the normal fits went well (NB. Good analyses have the value 0, that's why I use the logical not)
        analyses_to_take = np.ones(len(is_good_analysis)*2) #I create an array full of False elements
        cont_in_is_good_analysis = 0

        for index_ele_to_take,ele_to_take in enumerate(analyses_to_take):
            if not (index_ele_to_take % 2): #for those analyses that were done with the normal mask
                analyses_to_take[index_ele_to_take] = is_good_analysis[cont_in_is_good_analysis] #I add the value of is_good_analysis for this very element
                cont_in_is_good_analysis = cont_in_is_good_analysis+1
        #the final analyses_to_take has the form [is_good_analysis[0], 1, is_good_analysis[1], 1, ...]
        #by doing the logical_not all the good analyses (0 value) will became 1 and the others 0, and I can convert everything to booleans

        analyses_to_take = np.logical_not(analyses_to_take)

        analyses_to_take_boolean = analyses_to_take.astype(bool) #I convert it to boolean
        chi2_to_take = chi2_array[analyses_to_take_boolean] #I take into account only chi2_to_take
        if len(chi2_to_take) > 1:
            filter_chi2 = chi2_array == chi2_to_take.min() #I get the minimun in chi2_to_take
        else:
            filter_chi2 = analyses_to_take_boolean #if only a single valid element, I take it as the best analysis
        filter_result = filter_chi2
    else:
        if len(chi2nu_array) != 0:
            filter_chi2nu = chi2nu_array == chi2nu_array.min() #if no good analysis for the normal masks, I take all masks
            filter_result = filter_chi2nu


    index_best_result = np.where(filter_result)
    index_best_result = index_best_result[0] #this way, I remove the tuple I get with the output of np.where
    if index_best_result.size > 1: index_best_result = index_best_result[0] #this way I am sure I get the first element if there are several
    index_best_result = int(index_best_result)  #to be sure it is an integer
    index = index_best_result #I will repeat this index a lot, so I make the name shorter


    vector_gal_name = np.append(vector_gal_name,galaxy["id"])
    vector_zz       = np.append(vector_zz      ,galaxy["zz"])
    kpc_per_arcsec = 1./cosmo.arcsec_per_kpc_proper(galaxy["zz"])
    #is there any object in the central pixels of the image?
    img = fits.open("../../galaxy_images/"+galaxy["id"]+".fits")
    hdr = img[0].header
    xc = hdr["NAXIS1"]/2.
    yc = hdr["NAXIS2"]/2.
    #-------
    if single_or_double == 1:
        dist_target = np.sqrt( (xc-float(galaxy["xc"][index]))**2. + (yc-float(galaxy["yc"][index]))**2. )
    else:
        dist_target = np.sqrt( (xc-float(galaxy["xc1"][index]))**2. + (yc-float(galaxy["yc1"][index]))**2. )
    #-------
    if dist_target >  5.:
        flag_obj_detected = False
    else:
        flag_obj_detected = True
    vector_flag_obj_detected = np.append(vector_flag_obj_detected,flag_obj_detected)

    if flag_obj_detected == True:
        if single_or_double == 1:
            vector_mag = np.append(vector_mag, galaxy["mag"][index] )
            vector_re  = np.append(vector_re,  galaxy["re"][index] )
            vector_nn  = np.append(vector_nn,  galaxy["nn"][index] )
            vector_ar  = np.append(vector_ar,  galaxy["ar"][index] )
            vector_pa  = np.append(vector_pa,  galaxy["pa"][index] )
            vector_magerr = np.append(vector_magerr, np.std(galaxy["mag"]) )
            vector_reerr  = np.append(vector_reerr,  np.std(galaxy["re"]) )
            vector_nnerr  = np.append(vector_nnerr,  np.std(galaxy["nn"]) )
            vector_arerr  = np.append(vector_arerr,  np.std(galaxy["ar"]) )
            vector_paerr  = np.append(vector_paerr,  np.std(galaxy["pa"]) )
            vector_re_kpc     = np.append(vector_re_kpc,      galaxy["re"][index]*pix_scale*kpc_per_arcsec)
            vector_reerr_kpc  = np.append(vector_reerr_kpc,  np.std(galaxy["re"])*pix_scale*kpc_per_arcsec)
        else:
            vector_mag1 = np.append(vector_mag1, galaxy["mag1"][index] )
            vector_re1  = np.append(vector_re1,  galaxy["re1"][index] )
            vector_nn1  = np.append(vector_nn1,  galaxy["nn1"][index])
            vector_ar1  = np.append(vector_ar1,  galaxy["ar1"][index] )
            vector_pa1  = np.append(vector_pa1,  galaxy["pa1"][index] )
            vector_mag1err = np.append(vector_mag1err, np.std(galaxy["mag1"]) )
            vector_reerr1  = np.append(vector_reerr1,  np.std(galaxy["re1"]) )
            vector_nnerr1  = np.append(vector_nnerr1,  np.std(galaxy["nn1"]) )
            vector_arerr1  = np.append(vector_arerr1,  np.std(galaxy["ar1"]) )
            vector_paerr1  = np.append(vector_paerr1,  np.std(galaxy["pa1"]) )
            vector_mag2 = np.append(vector_mag2, galaxy["mag2"][index] )
            vector_re2 = np.append(vector_re2,  galaxy["re2"][index] )
            vector_nn2 = np.append(vector_nn2,  galaxy["nn2"][index] )
            vector_ar2 = np.append(vector_ar2,  galaxy["ar2"][index] )
            vector_pa2 = np.append(vector_pa2,  galaxy["pa2"][index] )
            vector_magerr2 = np.append(vector_magerr2, np.std(galaxy["mag2"]) )
            vector_reerr2 = np.append(vector_reerr2,  np.std(galaxy["re2"]) )
            vector_nnerr2 = np.append(vector_nnerr2,  np.std(galaxy["nn2"]) )
            vector_arerr2 = np.append(vector_arerr2,  np.std(galaxy["ar2"]) )
            vector_paerr2 = np.append(vector_paerr2,  np.std(galaxy["pa2"]) )
            vector_re1_kpc     = np.append(vector_re_kpc1,      galaxy["re1"][index]*pix_scale*kpc_per_arcsec)
            vector_re2_kpc     = np.append(vector_re_kpc2,      galaxy["re2"][index]*pix_scale*kpc_per_arcsec)
            vector_reerr1_kpc  = np.append(vector_reerr1_kpc,  np.std(galaxy["re1"])*pix_scale*kpc_per_arcsec)
            vector_reerr2_kpc  = np.append(vector_reerr2_kpc,  np.std(galaxy["re2"])*pix_scale*kpc_per_arcsec)
    else:
        if single_or_double == 1:
            vector_mag = np.append(vector_mag, float('NaN') )
            vector_re  = np.append(vector_re,  float('NaN') )
            vector_nn  = np.append(vector_nn,  float('NaN') )
            vector_ar  = np.append(vector_ar,  float('NaN') )
            vector_pa  = np.append(vector_pa,  float('NaN') )
            vector_magerr = np.append(vector_magerr, float('NaN') )
            vector_reerr  = np.append(vector_reerr,  float('NaN') )
            vector_nnerr  = np.append(vector_nnerr,  float('NaN') )
            vector_arerr  = np.append(vector_arerr,  float('NaN') )
            vector_paerr  = np.append(vector_paerr,  float('NaN') )
            vector_re_kpc = np.append(vector_re_kpc, float('NaN') )
            vector_reerr_kpc = np.append(vector_reerr_kpc, float('NaN') )
        else:
            vector_mag1 = np.append(vector_mag1, float('NaN') )
            vector_re1  = np.append(vector_re1,  float('NaN') )
            vector_nn1  = np.append(vector_nn1,  float('NaN') )
            vector_ar1  = np.append(vector_ar1,  float('NaN') )
            vector_pa1  = np.append(vector_pa1,  float('NaN') )
            vector_mag1err = np.append(vector_mag1err, float('NaN') )
            vector_reerr1  = np.append(vector_reerr1,  float('NaN') )
            vector_nnerr1  = np.append(vector_nnerr1,  float('NaN') )
            vector_arerr1  = np.append(vector_arerr1,  float('NaN') )
            vector_paerr1  = np.append(vector_paerr1,  float('NaN') )
            vector_mag2 = np.append(vector_mag2,  float('NaN') )
            vector_re2  = np.append(vector_re2,  float('NaN') )
            vector_nn2  = np.append(vector_nn2,  float('NaN') )
            vector_ar2  = np.append(vector_ar2,  float('NaN') )
            vector_pa2  = np.append(vector_pa2,  float('NaN') )
            vector_magerr2 = np.append(vector_magerr2, float('NaN') )
            vector_reerr2  = np.append(vector_reerr2,  float('NaN') )
            vector_nnerr2  = np.append(vector_nnerr2,  float('NaN') )
            vector_arerr2  = np.append(vector_arerr2,  float('NaN') )
            vector_paerr2  = np.append(vector_paerr2,  float('NaN') )
            vector_re1_kpc     = np.append(vector_re1_kpc,     float('NaN') )
            vector_reerr1_kpc  = np.append(vector_reerr1_kpc,  float('NaN') )
            vector_re2_kpc     = np.append(vector_re2_kpc,     float('NaN') )
            vector_reerr2_kpc  = np.append(vector_reerr2_kpc,  float('NaN') )


#writing the final catalogs
if single_or_double == 1:
    data = Table( [vector_gal_name,vector_zz,vector_mag,vector_magerr,vector_re,vector_reerr,vector_re_kpc,vector_reerr_kpc,vector_nn,vector_nnerr,\
                   vector_ar,vector_arerr,vector_pa,vector_paerr,vector_flag_obj_detected], \
                   names=("gal_id","zz","mag","mag_error","re_pix","re_pix_error","re_kpc","re_kpc_err","n","n_error",\
                         "ar","ar_error","pa","pa_error","flag_obj_detected") )
else:
    data = Table( [vector_gal_name,vector_zz,vector_mag1,vector_magerr1,vector_re1,vector_reerr1,vector_re1_kpc,vector_reerr1_kpc,vector_nn1,vector_nnerr,\
                   vector_ar1,vector_arerr1,vector_pa1,vector_paerr1,vector_mag2,vector_magerr2,vector_re2,vector_reerr2,\
                   vector_re2_kpc,vector_reerr2_kpc,vector_nn2,vector_nnerr2,vector_ar2,vector_arerr2,vector_pa2,vector_paerr2,vector_flag_obj_detected], \
                   names=("gal_id","zz","mag1","mag_error1","re_pix1","re_pix_error1","re_pix1_kpc","re_kpc_error1","n1","n_error1",\
                         "ar1","ar_error1","pa1","pa_error1", "mag2","mag_error2","re_pix2","re_pix_error2",\
                         "n2","n_error2","ar2","ar_error2","pa2","pa_error2","flag_obj_detected") )
data.write('results_cpu'+cpu_number+'.cat', format='ascii.commented_header')
