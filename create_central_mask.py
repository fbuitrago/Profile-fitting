from astropy.io import fits
import numpy as np
from dist_ellipse import dist_ellipse

def mask_obj(primary_mask,x_size,y_size,xc,yc,axis_ratio,position_angle,a_image,kron_radius,enlargemask):
    new_img = dist_ellipse([x_size,y_size], xc, yc, axis_ratio, position_angle)

    first_size = a_image*kron_radius
    if first_size > 5.:
        final_size = first_size*enlargemask 
    else:
        final_size = 5.*enlargemask

    np.putmask(new_img,new_img < final_size,1) #pixels I will mask
    np.putmask(new_img,new_img > final_size,0) #pixels I will not mask

    np.logical_or(primary_mask,new_img,primary_mask) #v1,v2,result

    return(primary_mask)

def unmask_obj(primary_mask,x_size,y_size,xc,yc,axis_ratio,position_angle,a_image,kron_radius,enlargemask):
    new_img = dist_ellipse([x_size,y_size], xc, yc, axis_ratio, position_angle)

    first_size = a_image*kron_radius
    if first_size > 5.:
        final_size = first_size*enlargemask 
    else:
        final_size = 5.*enlargemask

    np.putmask(new_img,new_img < final_size,0) #pixels I will mask
    np.putmask(new_img,new_img > final_size,1) #pixels I will not mask

    np.logical_and(primary_mask,new_img,primary_mask) #v1,v2,result

    return(primary_mask)

def create_central_mask(galaxy,sex_output,target,enlargemask):
    objs  = sex_output[0,:]
    xc    = sex_output[5,:]
    yc    = sex_output[6,:]
    ellip = sex_output[8,:]
    theta = sex_output[9,:]
    kron  = sex_output[10,:]
    a_image = sex_output[11,:]

    img = fits.open("../../galaxy_images/"+galaxy+".fits")
    header = img[0].header
    x_size = header['NAXIS1']
    y_size = header['NAXIS2']

    primary_mask = np.ones( (x_size,y_size), dtype=bool)

    #unmasking the central region
    primary_mask = unmask_obj(primary_mask,x_size,y_size,xc[target],yc[target],1./(1.-ellip[target]), theta[target],a_image[target],kron[target],enlargemask)

    #masking again all objects

    primary_mask = np.array(primary_mask,dtype=int)
    fits.writeto("./"+galaxy+"/"+galaxy+'_central_mask.fits',primary_mask,header,clobber=True)
