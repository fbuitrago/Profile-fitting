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

def create_mask(galaxy,sex_output,target,enlargemask,filter_objs):
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

    primary_mask = np.zeros( (x_size,y_size), dtype=bool)

    #masking all objects
    for ii in range(len(objs)):
        primary_mask = mask_obj(primary_mask,x_size, y_size, xc[ii], yc[ii], 1./(1.-ellip[ii]), theta[ii],a_image[ii],kron[ii],enlargemask)

    #unmasking the central region
    primary_mask = unmask_obj(primary_mask,x_size,y_size,xc[target],yc[target],1./(1.-ellip[target]), theta[target],a_image[target],kron[target],enlargemask)

    #unmasking selected neighbours
    for ii in range(len(filter_objs)):
        if (ii != target) & (filter_objs[ii]==True):
            primary_mask = unmask_obj(primary_mask,x_size,y_size,xc[ii],yc[ii],1./(1.-ellip[ii]), theta[ii],a_image[ii],kron[ii],enlargemask)

    primary_mask = np.array(primary_mask,dtype=int)
    fits.writeto("./"+galaxy+"/"+galaxy+'_normal_mask.fits',primary_mask,header,clobber=True)


