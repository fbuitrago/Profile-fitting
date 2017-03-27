import numpy as np
import pdb

def write_galfit_script(galaxy,star,mask,x_size,y_size,zeropoint,exptime,pixel_scale,sex_output,target,single_or_double,bulge_disk_decomp,filter_objs):
    mag   = sex_output[3,:]
    xc    = sex_output[5,:]
    yc    = sex_output[6,:]
    ellip = sex_output[8,:]
    theta = sex_output[9,:]
    kron  = sex_output[10,:]
    a_image = sex_output[11,:]

    galfit_zeropoint = zeropoint-2.5*np.log10(exptime)
    re_target = (a_image[target]*kron[target])/2.
    axis_ratio = 1.-ellip
    theta = theta-90. #moving SExtractor reference frame to GALFIT reference frame
    for ii in range(len(theta)): #due to the GALFIT values for the position angle (-90<pa<90)
        if theta[ii] < -90.: theta[ii] = theta[ii]+180.
        if theta[ii] >  90.: theta[ii] = theta[ii]-180.

    galfit_script = str(galaxy)+"_"+str(star)+"_"+str(mask)+".script" #a unique GALFIT script for every combination of input parameters

    file = open("./"+galaxy+"/"+galfit_script, "w")

    file.write("# IMAGE PARAMETERS\n")
    file.write(" A) ../../../galaxy_images/"+str(galaxy)+".fits  # Input Data image (FITS file)\n")
    file.write(" B) "+str(galaxy)+"_"+str(star)+"_"+str(mask)+".fits  # Name for the output image\n")
    file.write(" C) ../../../sigma_images/"+str(galaxy)+"_rms.fits  # Noise image name (made from data if blank or 'none')\n")
    file.write(" D) ../../../star_images/"+str(star)+".fits # Input PSF image and (optional) diffusion kernel\n")
    file.write(" E) 1  # PSF oversampling factor relative to data\n")
    file.write(" F) "+str(galaxy)+"_"+str(mask)+".fits  # Pixel mask (ASCII file or FITS file with non-0 values)\n")
    file.write(" G) constraints"+"_"+str(mask)+".txt  # Parameter constraint file (ASCII)\n")
    file.write(" H) 1         "+str(x_size)+ "  1          "+str(y_size)+ "  # Image region to fit (xmin xmax ymin ymax)\n")
    file.write(" I) "+str(x_size)+"      "+str(y_size)+"  # Size of convolution box (x y)\n")
    file.write(" J) "+str(galfit_zeropoint)+"  # Magnitude photometric zeropoint\n")
    file.write(" K) "+str(pixel_scale)+"     "+str(pixel_scale)+"  # Plate scale (dx dy)\n")
    file.write(" O) regular  # Display type (regular, curses, both)\n")
    file.write(" P) 0  # Create ouput only? (1=yes; 0=optimize)\n")
    file.write(" S) 0  # Modify/create objects interactively?\n")
    file.write("\n")
    file.write("\n")
    file.write("# Component number: 1\n")
    file.write(" 0) sky                    #  Component type\n")
    file.write(" 1) 0.00      1       #  Sky background at center of fitting region [ADUs]\n")
    file.write(" 2) 0.000     1       #  dsky/dx (sky gradient in x)     [ADUs/pix]\n")
    file.write(" 3) 0.000     1       #  dsky/dy (sky gradient in y)     [ADUs/pix]\n")
    file.write(" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
    file.write("\n")
    file.write("# Component number: 2\n")
    file.write(" 0) sersic                 #  Component type\n")
    file.write(" 1) "+str(xc[target])+" "+str(yc[target])+" 1 1  #  Position x, y\n")
    file.write(" 3) "+str(mag[target])+"     1          #  Integrated magnitude\n") 
    file.write(" 4) "+str(re_target)+"      1          #  R_e (effective radius)   [pix]\n")
    file.write(" 5) 1.           1          #  Sersic index n (de Vaucouleurs n=4)\n") 
    file.write(" 9) "+str(axis_ratio[target])+"      1          #  Axis ratio (b/a)  \n")
    file.write("10) "+str(theta[target])+"    1          #  Position angle (PA) [deg: Up=0, Left=90]\n")
    file.write(" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
    file.write("\n")

    #second component for target galaxy
    if single_or_double == 2:
        file.write("# Component number: 3\n")
        file.write(" 0) sersic                 #  Component type\n")
        file.write(" 1) "+str(xc[target])+" "+str(yc[target])+" 1 1  #  Position x, y\n")
        file.write(" 3) "+str(mag[target]+1.)+"     1          #  Integrated magnitude\n") 
        file.write(" 4) "+str(re_target - (re_target/5.))+"      1          #  R_e (effective radius)   [pix]\n")
        if bulge_disk_decomp == True:
            file.write(" 5) 4.           0          #  Sersic index n (de Vaucouleurs n=4)\n") 
        else:
            file.write(" 5) 1.5           1          #  Sersic index n (de Vaucouleurs n=4)\n") 
        file.write(" 9) "+str(axis_ratio[target])+"      1          #  Axis ratio (b/a)  \n")
        file.write("10) "+str(theta[target])+"    1          #  Position angle (PA) [deg: Up=0, Left=90]\n")
        file.write(" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
        file.write("\n")
 
    #galaxy neighbours
    order = single_or_double+1 #the first +1 is for the sky
    for ii in range(len(filter_objs)):
        if (ii != target) & (filter_objs[ii]==True) & (mask.find("normal")!=-1):
            order = order+1 #the first +1 is for the sky
            file.write("# Component number: "+str(order)+"\n")
            file.write(" 0) sersic                 #  Component type\n")
            file.write(" 1) "+str(xc[ii])+" "+str(yc[ii])+" 1 1  #  Position x, y\n")
            file.write(" 3) "+str(mag[ii])+"     1          #  Integrated magnitude\n") 
            file.write(" 4) "+str((a_image[ii]*kron[ii])/4.)+"      1          #  R_e (effective radius)   [pix]\n")
            file.write(" 5) 1.           1          #  Sersic index n (de Vaucouleurs n=4)\n") 
            file.write(" 9) "+str(axis_ratio[ii])+"      1          #  Axis ratio (b/a)  \n")
            file.write("10) "+str(theta[ii])+"    1          #  Position angle (PA) [deg: Up=0, Left=90]\n")
            file.write(" Z) 0                      #  Skip this model in output image?  (yes=1, no=0)\n")
            file.write("\n")

    file.close()
    return(galfit_script)
