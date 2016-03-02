"""

Return the MUSE FOV as a patch that can be overplotted on a image


scale - > need to convert arcsec in pixels for a given image

fov - return the field of view
sgs - return the slow guiding sytem area
ao  - return the AO tip-tilt area
xc,yc - central pixel 
lwidth - line width

"""
import matplotlib.patches as patches


def muse_fov(scale,xc=1,yc=1,fov=True,sgs=True,ao=True,lwidth=2,lstyle='-'):

    patch_list=[]

    
    if(fov):
        #add muse fov 
        p=patches.Rectangle((xc-30./scale,yc-30./scale),58./scale,58./scale,\
                                edgecolor='red',linewidth=lwidth,facecolor='none',linestyle=lstyle)
        patch_list.append(p)
       

    if(ao):
        #add tiptilt area
        p=patches.Circle((xc,yc),104.*0.5/scale,edgecolor='blue',linewidth=lwidth,facecolor='none',linestyle=lstyle)
        patch_list.append(p)
        p=patches.Circle((xc,yc),203.*0.5/scale,edgecolor='blue',linewidth=lwidth,facecolor='none',linestyle=lstyle)
        patch_list.append(p)

    if(sgs):
    
        #add guide area
        xguide=xc-30.*0.5/scale
        yguide=yc+37./scale-8.*0.5/scale
        xgw=30./scale
        ygw=8./scale
        p=patches.Rectangle((xguide,yguide),xgw,ygw,edgecolor='green',linewidth=lwidth,facecolor='none',linestyle=lstyle)
        patch_list.append(p)
        
        xguide=xc-30.*0.5/scale
        yguide=yc-37./scale-8.*0.5/scale
        xgw=30./scale
        ygw=8./scale
        p=patches.Rectangle((xguide,yguide),xgw,ygw,edgecolor='green',linewidth=lwidth,facecolor='none',linestyle=lstyle)
        patch_list.append(p)
        
        yguide=yc-30.*0.5/scale
        xguide=xc-37./scale-8.*0.5/scale
        ygw=30./scale
        xgw=8./scale
        p=patches.Rectangle((xguide,yguide),xgw,ygw,edgecolor='green',linewidth=lwidth,facecolor='none',linestyle=lstyle)
        patch_list.append(p)
        
        yguide=yc-30.*0.5/scale
        xguide=xc+37./scale-8.*0.5/scale
        ygw=30./scale
        xgw=8./scale
        p=patches.Rectangle((xguide,yguide),xgw,ygw,edgecolor='green',linewidth=lwidth,facecolor='none',linestyle=lstyle)
        patch_list.append(p)


    
    return patch_list
