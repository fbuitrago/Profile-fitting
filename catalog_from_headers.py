from astropy.io import fits
import my_python_library
import numpy as np
import os
import pdb
import string
from astropy.table import Table


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


#knowing how many elements my datacube will have
stars = np.array(["star_h_band_candels"],dtype=str) #CHANGE IF NECESSARY #don't add extension .fits
#=============================
flag_mask_central = True #CHANGE IF NECESSARY #if False we skip the creation of central mask (a mask containing only the central obj, good for overcrowded fields)
masks = []
if flag_mask_central == True:
    masks = ["normal_mask","central_mask"]
else:
    masks = ["normal_mask"]
masks = np.array(masks,dtype=str)
#==============================
single_or_double = 1  #CHANGE IF NECESSARY #single (1) or double (2) component fit
#==============================
ini_conds = "" #CHANGE IF NECESSARY #grid initial conditions "_1","_2",...

#knowing the galaxies to run GALFIT
cwd = os.getcwd()
path, last_folder = os.path.split(cwd)
cpu_number = last_folder.replace("cpu","")
ids,ra,dec = np.loadtxt("./aux_files/gals_cpu"+cpu_number+".cat", dtype = float, unpack = True)
galaxies = my_python_library.ids_as_str(ids)

#reading the headers
gal_list = []
for ii in range(len(galaxies)):
    #I start with an empty dictionary
    if single_or_double == 1:
        gal = {'id': -99, 'chi2': [], 'chi2nu': [], 'xc': [], 'yc': [], 'mag': [], 'mag_flag_special': [], 'mag_flag_fixed': [],
               're': [], 're_flag_special': [], 're_flag_fixed': [], 'nn': [], 'nn_flag_special': [], 'nn_flag_fixed': [], 
               'ar': [], 'ar_flag_special': [], 'ar_flag_fixed': [], 'pa': [], 'pa_flag_special': [], 'pa_flag_fixed': []}
    else:
        gal = {'id': -99, 'chi2': [], 'chi2nu': [], 'xc1': [], 'yc1': [], 'mag1': [], 'mag_flag_special1': [], 'mag_flag_fixed1': [],
               're1': [], 're_flag_special1': [], 're_flag_fixed1': [], 'nn1': [], 'nn_flag_special1': [], 'nn_flag_fixed1': [], 
               'ar1': [], 'ar_flag_special1': [], 'ar_flag_fixed1': [], 'pa1': [], 'pa_flag_special1': [], 'pa_flag_fixed1': [],
                                      'xc2': [], 'yc2': [], 'mag2': [], 'mag_flag_special2': [], 'mag_flag_fixed2': [],
               're2': [], 're_flag_special2': [], 're_flag_fixed2': [], 'nn2': [], 'nn_flag_special2': [], 'nn_flag_fixed2': [], 
               'ar2': [], 'ar_flag_special2': [], 'ar_flag_fixed2': [], 'pa2': [], 'pa_flag_special2': [], 'pa_flag_fixed2': []}
    gal['id'] = galaxies[ii]
    for jj in range(len(stars)):
        for kk in range(len(masks)):
            #for cont in range(len(ini_conds)):
            filename = "./"+galaxies[ii]+"/"+galaxies[ii]+"_"+stars[jj]+"_"+masks[kk]+".fits"
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

    chi2_vector = np.array( galaxy["chi2"] )
    filter_chi2 = chi2_vector == chi2_vector.min()
    index_best_chi2 = np.where(filter_chi2)

    if single_or_double == 1:
        is_good_analysis = gal["mag_flag_special"] + gal["re_flag_special"] + gal["nn_flag_special"] + gal["ar_flag_special"] + gal["pa_flag_special"]
    else:
        is_good_analysis = gal["mag_flag_special1"] + gal["re_flag_special1"] + gal["nn_flag_special1"] + gal["ar_flag_special1"] + gal["pa_flag_special1"] + \
                           gal["mag_flag_special2"] + gal["re_flag_special2"] + gal["nn_flag_special2"] + gal["ar_flag_special2"] + gal["pa_flag_special2"]

    #if possible, choose among the fits with the normal masks. If no good fit there, use chi2nu to choose the best fit among all fits
    if ( (int(index_best_chi2[0]) % 2)==0 and is_good_analysis==0 ):
        filter_result = filter_chi2
    else:
        chi2nu_vector = np.array( galaxy["chi2nu"] )
        filter_chi2nu = chi2nu_vector == chi2nu_vector.min()
        filter_result = filter_chi2nu
    index_best_result = np.where(filter_result)
    index = int(index_best_result[0])

    gal_name = galaxy["id"]
    img = fits.open("../galaxy_images/"+gal_name+".fits")
    hdr = img[0].header
    xc = hdr["NAXIS1"]
    yc = hdr["NAXIS2"]
    if single_or_double == 1:
        dist_target = np.sqrt( (xc-galaxy["xc"][index])**2. + (yc-galaxy["yc"][index])**2. )
    else:
        dist_target = np.sqrt( (xc-galaxy["xc1"][index])**2. + (yc-galaxy["yc1"][index])**2. )
    if dist_target >  5.:
        flag_obj_detected = False
    else:
        flag_obj_detected = True

    if flag_obj_detected == True:
        if single_or_double == 1:
            mag.append( galaxy["mag"][index] )
            re .append( galaxy["re"][index] )
            nn .append( galaxy["nn"][index] )
            ar .append( galaxy["ar"][index] )
            pa .append( galaxy["pa"][index] )
            magerr.append( np.std(galaxy["mag"]) )
            reerr .append( np.std(galaxy["re"]) )
            nnerr .append( np.std(galaxy["nn"]) )
            arerr .append( np.std(galaxy["ar"]) )
            paerr .append( np.std(galaxy["pa"]) )
        else:
            mag1.append( galaxy["mag1"][index] )
            re1 .append( galaxy["re1"][index] )
            nn1 .append( galaxy["nn1"][index])
            ar1 .append( galaxy["ar1"][index] )
            pa1 .append( galaxy["pa1"][index] )
            magerr1.append( np.std(galaxy["mag1"]) )
            reerr1 .append( np.std(galaxy["re1"]) )
            nnerr1 .append( np.std(galaxy["nn1"]) )
            arerr1 .append( np.std(galaxy["ar1"]) )
            paerr1 .append( np.std(galaxy["pa1"]) )
            mag2.append( galaxy["mag2"][index] )
            re2 .append( galaxy["re2"][index] )
            nn2 .append( galaxy["nn2"][index] )
            ar2 .append( galaxy["ar2"][index] )
            pa2 .append( galaxy["pa2"][index] )
            magerr2.append( np.std(galaxy["mag2"]) )
            reerr2 .append( np.std(galaxy["re2"]) )
            nnerr2 .append( np.std(galaxy["nn2"]) )
            arerr2 .append( np.std(galaxy["ar2"]) )
            paerr2 .append( np.std(galaxy["pa2"]) )
    else:
        if single_or_double == 1:
            mag.append( float('NaN') )
            re .append( float('NaN') )
            nn .append( float('NaN') )
            ar .append( float('NaN') )
            pa .append( float('NaN') )
            magerr.append( float('NaN') )
            reerr .append( float('NaN') )
            nnerr .append( float('NaN') )
            arerr .append( float('NaN') )
            paerr .append( float('NaN') )
        else:
            mag1.append( float('NaN') )
            re1 .append( float('NaN') )
            nn1 .append( float('NaN') )
            ar1 .append( float('NaN') )
            pa1 .append( float('NaN') )
            magerr1.append( float('NaN') )
            reerr1 .append( float('NaN') )
            nnerr1 .append( float('NaN') )
            arerr1 .append( float('NaN') )
            paerr1 .append( float('NaN') )
            mag2.append( float('NaN') )
            re2 .append( float('NaN') )
            nn2 .append( float('NaN') )
            ar2 .append( float('NaN') )
            pa2 .append( float('NaN') )
            magerr2.append( float('NaN') )
            reerr2 .append( float('NaN') )
            nnerr2 .append( float('NaN') )
            arerr2 .append( float('NaN') )
            paerr2 .append( float('NaN') )

#writing the final catalogs
if single_or_double == 1:
    data = Table( [gal_name,mag,magerr,re,reerr,nn,nnerr,ar,arerr,pa,paerr],
                  names=("gal_id","mag","mag_error","re_pix","re_pix_error","n","n_error","ar","ar_error","pa","pa_error") , Writer = ascii.CommentedHeader)
else:
    data = Table( [gal_name1,mag1,magerr1,re1,reerr1,nn1,nnerr,ar1,arerr1,pa1,paerr1,mag2,magerr2,re2,reerr2,nn2,nnerr2,ar2,arerr2,pa2,paerr2],
                  names=("gal_id1","mag1","mag_error1","re_pix1","re_pix_error1","n1","n_error1","ar1","ar_error1","pa1","pa_error1" \
                                   "mag2","mag_error2","re_pix2","re_pix_error2","n2","n_error2","ar2","ar_error2","pa2","pa_error2"), \
                  Writer = ascii.CommentedHeader)
data.write('results_cpu'+cpu_number+'.cat', format='ascii')

pdb.set_trace()
