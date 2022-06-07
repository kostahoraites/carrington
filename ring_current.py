from carrington import *
#from carrington import fg_grid, mkdir_path, numcurl3d, numjacobian3d, R_EARTH
from pyPlots import plot_vdf
import fieldmodels
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from myutils import *

global CELLSIZE_XYZ
#global CELLSIZE

np.seterr(divide='ignore', invalid='ignore')
# calculate current density j


if __name__ == "__main__":

    #run = 'EGL'             # default values, these may be overwritten depending on how the function is called
    #run = 'EGM'             # default values, these may be overwritten depending on how the function is called
    #run = 'EGN'             # default values, these may be overwritten depending on how the function is called
    #run = 'EGO'             # default values, these may be overwritten depending on how the function is called
    #run = 'EGL'
    #run = 'EGP'
    run = 'EGILIKE'             # default values, these may be overwritten depending on how the function is called
    dim = '3D'

                            # make one plot for a single frame in every delta_nframes
    if run == 'EGL':
       delta_nframes = 20      # EGL, 
    elif run == 'EGP':
       delta_nframes = 1      # EGP
    else: 
       delta_nframes = 1

    # Frame extent for this job given as command-line arguments
    if len(sys.argv)==3:    # Starting and end frames given
        fileIndex_list = range(int(sys.argv[1]) * delta_nframes, int(sys.argv[2]) * delta_nframes, delta_nframes)
    elif len(sys.argv)==2:  # Only starting frame given, generate one frame
        print(int(sys.argv[1]))
        fileIndex_list = range(int(sys.argv[1]), int(sys.argv[1])+1, 1)
    else:    # no inputs given
        # load data
        if run == 'EGL':
            fileIndex_list = [1760]     # EGL
        elif run == 'EGM':
            fileIndex_list = [1247]    # EGM
        elif run == 'EGN':
            fileIndex_list = [488]     # EGN
        elif run == 'EGO':
            fileIndex_list = [151]     # EGO
        elif run == 'EGP':
            fileIndex_list = [299]     # bulk1, EGP (available fileIndex: 269-299)
            #fileIndex_list = [53]     # bulk5, EGP (available fileIndex: 1-53)
        elif run == 'EGILIKE':
            fileIndex_list = [719]     # bulk, EGP (available fileIndex: 718-719)

       # dipole field
    tilt_angle_phi = 0.
    tilt_angle_theta = 0.

    if run == 'EGL':
       dipoleXFull = 9.5565e7 # 15 RE
       dipoleXZero = 2.5e8
    elif run == 'EGP':
       dipoleXFull = 9.5565e7 # 15 RE 
       dipoleXZero = 2.0e8
    else:
       dipoleXFull = None
       dipoleXZero = None

    #fieldmodels.dipole.set_dipole(centerx, centery, centerz, tilt_phi, tilt_theta, mult=1.0, radius_f=None, radius_z=None):
    #dip = fieldmodels.dipole(0,0,0,tilt_angle_phi,tilt_angle_theta)
    dip = fieldmodels.dipole(0,0,0,tilt_angle_phi,tilt_angle_theta, radius_f = dipoleXFull, radius_z = dipoleXZero)
    #mdip = fieldmodels.dipole(80*RE,0,0,tilt_angle_phi,180.-tilt_angle_theta)


    for fileIndex in fileIndex_list:

        filename = get_vlsvfile_fullpath(run, fileIndex)
        f=pt.vlsvfile.VlsvReader(filename)
        fg_b = f.read_variable('fg_b')
        #CELLSIZE = (f.read_parameter('xmax') - f.read_parameter('xmin')) / (fg_b.shape[0] - 1)    ### COMMENT THIS OUT!!!
        CELLSIZE_XYZ = [ (f.read_parameter('xmax') - f.read_parameter('xmin')) / fg_b.shape[0],
                         (f.read_parameter('ymax') - f.read_parameter('ymin')) / fg_b.shape[1],
                         (f.read_parameter('zmax') - f.read_parameter('zmin')) / fg_b.shape[2] ]
        x, y, z = fg_grid(f, fg_b = fg_b)
        zind_temp = np.where(np.abs(z) == np.nanmin(np.abs(z)))[0][0]
        zeq0_ind = 2
        fg_b_analyze = fg_b[:,:,zind_temp-zeq0_ind:zind_temp+zeq0_ind+1,:]   # reduce the data size so that numcurl3d doesn't require too much memory
        J = (1 / mu_0) * numcurl3d(fg_b_analyze, CELLSIZE_XYZ)
        #J_magnitude = (J[:,:,:,0]**2 + J[:,:,:,1]**2 + J[:,:,:,2]**2)**0.5
        Jxy_magnitude = (J[:,:,:,0]**2 + J[:,:,:,1]**2 )**0.5
        B = (fg_b_analyze[:,:,:,0]**2 + fg_b_analyze[:,:,:,1]**2 + fg_b_analyze[:,:,:,2]**2)**0.5
        x2_temp, y2_temp = np.meshgrid(x, y, indexing='ij', sparse=True)
        x2d = x2_temp + (y2_temp * 0)
        y2d = y2_temp + (x2_temp * 0)
        r_2d = (x2d**2 + y2d**2)**0.5
        vg_r_min = 4.7
        Jxy_magnitude[r_2d/R_EARTH < vg_r_min] = 0      # zero out anything outside simulation domain


        shape = fg_b_analyze.shape
        B_dip = np.zeros(shape)
   
        for i in range(shape[0]):   # x
           for j in range(shape[1]):  # y 
              for k in range(shape[2]):  # z
                 for l in range(shape[3]):   #3 v-components
                    B_dip[i,j,k,l] = dip.get(x[i],y[j],z[ zind_temp-zeq0_ind+k ],0,l,0)
                    #B_dip[i,j,k,l] = dip.get_old(x[i],y[j],z[ zind_temp-zeq0_ind+k ],0,l,0)

        B_dip_mag = (B_dip[:,:,:,0]**2 + B_dip[:,:,:,1]**2 + B_dip[:,:,:,2]**2)**0.5
        J_dip = (1 / mu_0) * numcurl3d(B_dip, CELLSIZE_XYZ)
        Jxy_dip_magnitude = (J_dip[:,:,:,0]**2 + J_dip[:,:,:,1]**2 )**0.5
        Jxy_dip_magnitude[r_2d/R_EARTH < vg_r_min] = 0      # zero out anything outside simulation domain

        J_pert = J - J_dip
        Jxy_pert_magnitude = (J_pert[:,:,:,0]**2 + J_pert[:,:,:,1]**2 )**0.5
        Jxy_pert_magnitude[r_2d/R_EARTH < vg_r_min] = 0      # zero out anything outside simulation domain


        if run == 'EGL' or run == 'EGP':
            dq = 7        # index spacing between quivers
        else:
            dq = 1        # index spacing between quivers
        fig, ((ax3, ax4), (ax5, ax6) ) = plt.subplots(2,2)
        xlim = [-15, 15]
        ylim = [-15, 15]
        indx = np.where((x/R_EARTH >= xlim[0]) & (x/R_EARTH <= xlim[1]))[0]
        ix1 = np.nanmin(indx); ix2 = np.nanmax(indx)
        minx = np.nanmin(x[indx]); maxx = np.nanmax(x[indx])
        indy = np.where((y/R_EARTH >= ylim[0]) & (y/R_EARTH <= ylim[1]))[0]
        iy1 = np.nanmin(indy); iy2 = np.nanmax(indy)
        miny = np.nanmin(y[indy]); maxy = np.nanmax(y[indy])


           # plot the Jxy current
        type(Jxy_magnitude)
        type(J)
        Jxy_magnitude.shape
        J.shape
        Jxy_magnitude_plot = np.log10(Jxy_magnitude[ix1:ix2,iy1:iy2,zeq0_ind])
        Jx_q_plot = J[ix1:ix2,iy1:iy2,zeq0_ind, 0]/Jxy_magnitude[ix1:ix2,iy1:iy2, zeq0_ind]
        Jy_q_plot = J[ix1:ix2,iy1:iy2,zeq0_ind, 1]/Jxy_magnitude[ix1:ix2,iy1:iy2, zeq0_ind]
        for i in range(dq-1):
            Jx_q_plot[i+1::dq,:] = np.nan; Jx_q_plot[:,i+1::dq] = np.nan
            Jy_q_plot[i+1::dq,:] = np.nan; Jy_q_plot[:,i+1::dq] = np.nan
        im3 = ax3.pcolormesh(x2d[ix1:ix2,iy1:iy2]/R_EARTH, y2d[ix1:ix2,iy1:iy2]/R_EARTH, Jxy_magnitude_plot, shading = 'auto', cmap = 'plasma', vmin = -15, vmax = -8)
        divider = make_axes_locatable(ax3)
        cax3 = divider.append_axes('right', size='5%', pad=0.05)
        cbar3 = fig.colorbar(im3, cax=cax3, orientation='vertical')
        #cbar3.set_label(r'$log10(Jxy [A/m^2])$')
         # now quiver
        q = ax3.quiver(x2d[ix1:ix2:dq,iy1:iy2:dq]/R_EARTH, y2d[ix1:ix2:dq,iy1:iy2:dq]/R_EARTH, Jx_q_plot[0::dq,0::dq], Jy_q_plot[0::dq,0::dq])   # scale = dq*10?
        ax3.set_title(r'$log10(Jxy [A/m^2])$, z=0, run={}, time={}'.format(run, fileIndex) )
        ax3.set_xlim(xlim)
        ax3.set_ylim(ylim)
        #ax3.set_xlabel(r'x [$R_E$]')
        ax3.set_ylabel(r'y [$R_E$]')
        ax3.set_aspect('equal', adjustable='box')

           # plot the Jxy current for the dipole field
        Jxy_dip_magnitude_plot = np.log10(Jxy_dip_magnitude[ix1:ix2,iy1:iy2,zeq0_ind])
        Jx_dip_q_plot = J_dip[ix1:ix2,iy1:iy2,zeq0_ind, 0]/Jxy_dip_magnitude[ix1:ix2,iy1:iy2, zeq0_ind]
        Jy_dip_q_plot = J_dip[ix1:ix2,iy1:iy2,zeq0_ind, 1]/Jxy_dip_magnitude[ix1:ix2,iy1:iy2, zeq0_ind]
        for i in range(dq-1):
            Jx_dip_q_plot[i+1::dq,:] = np.nan; Jx_dip_q_plot[:,i+1::dq] = np.nan
            Jy_dip_q_plot[i+1::dq,:] = np.nan; Jy_dip_q_plot[:,i+1::dq] = np.nan
        im4 = ax4.pcolormesh(x2d[ix1:ix2,iy1:iy2]/R_EARTH, y2d[ix1:ix2,iy1:iy2]/R_EARTH, Jxy_dip_magnitude_plot, shading = 'auto', cmap = 'plasma', vmin = -15, vmax = -8 )
        divider = make_axes_locatable(ax4)
        cax4 = divider.append_axes('right', size='5%', pad=0.05)
        cbar4 = fig.colorbar(im4, cax=cax4, orientation='vertical')
        cbar4.set_label(r'$log10(Jxy [A/m^2])$')
         # now quiver
        q = ax4.quiver(x2d[ix1:ix2:dq,iy1:iy2:dq]/R_EARTH, y2d[ix1:ix2:dq,iy1:iy2:dq]/R_EARTH, Jx_dip_q_plot[0::dq,0::dq], Jy_dip_q_plot[0::dq,0::dq])   # scale = dq*10?
        ax4.set_title("'', dipole" )
        ax4.set_xlim(xlim)
        ax4.set_ylim(ylim)
        #ax4.set_xlabel(r'x [$R_E$]')
        #ax4.set_ylabel(r'y [$R_E$]')
        ax4.set_aspect('equal', adjustable='box')

           # plot the Jxy current for the perturbed field
        Jxy_pert_magnitude_plot = np.log10( abs(Jxy_pert_magnitude[ix1:ix2,iy1:iy2,zeq0_ind] ) )
        Jx_pert_q_plot = J_pert[ix1:ix2,iy1:iy2,zeq0_ind, 0]/Jxy_pert_magnitude[ix1:ix2,iy1:iy2, zeq0_ind]
        Jy_pert_q_plot = J_pert[ix1:ix2,iy1:iy2,zeq0_ind, 1]/Jxy_pert_magnitude[ix1:ix2,iy1:iy2, zeq0_ind]
        for i in range(dq-1):
            Jx_pert_q_plot[i+1::dq,:] = np.nan; Jx_pert_q_plot[:,i+1::dq] = np.nan
            Jy_pert_q_plot[i+1::dq,:] = np.nan; Jy_pert_q_plot[:,i+1::dq] = np.nan
        im5 = ax5.pcolormesh(x2d[ix1:ix2,iy1:iy2]/R_EARTH, y2d[ix1:ix2,iy1:iy2]/R_EARTH, Jxy_pert_magnitude_plot, shading = 'auto', cmap = 'plasma', vmin = -15, vmax = -8)
        divider = make_axes_locatable(ax5)
        cax5 = divider.append_axes('right', size='5%', pad=0.05)
        cbar5 = fig.colorbar(im5, cax=cax5, orientation='vertical')
        #cbar5.set_label(r'$log10(Jxy [A/m^2])$')
         # now quiver
        q = ax5.quiver(x2d[ix1:ix2:dq,iy1:iy2:dq]/R_EARTH, y2d[ix1:ix2:dq,iy1:iy2:dq]/R_EARTH, Jx_pert_q_plot[0::dq,0::dq], Jy_pert_q_plot[0::dq,0::dq])   # scale = dq*10?
        ax5.set_title("'', perturbed" )
        ax5.set_xlim(xlim)
        ax5.set_ylim(ylim)
        ax5.set_xlabel(r'x [$R_E$]')
        ax5.set_ylabel(r'y [$R_E$]')
        ax5.set_aspect('equal', adjustable='box')

           # dipole relative to total xy field
        Jxy_dip_ratio_magnitude_plot = np.log10( Jxy_magnitude[ix1:ix2,iy1:iy2,zeq0_ind] / Jxy_dip_magnitude[ix1:ix2,iy1:iy2,zeq0_ind]  )
        im6 = ax6.pcolormesh(x2d[ix1:ix2,iy1:iy2]/R_EARTH, y2d[ix1:ix2,iy1:iy2]/R_EARTH, Jxy_dip_ratio_magnitude_plot, shading = 'auto', cmap = 'bwr', vmin = -5, vmax = 5)
        divider = make_axes_locatable(ax6)
        cax6 = divider.append_axes('right', size='5%', pad=0.05)
        cbar6 = fig.colorbar(im6, cax=cax6, orientation='vertical')
        cbar6.set_label(r'$log10(J_{xy, total} / J_{xy, dip})$')
         # now quiver
        ax6.set_title("ratio (dipole/total)" )
        ax6.set_xlim(xlim)
        ax6.set_ylim(ylim)
        ax6.set_xlabel(r'x [$R_E$]')
        #ax6.set_ylabel(r'y [$R_E$]')
        ax6.set_aspect('equal', adjustable='box')

        save_dir = '{}{}/ring_current/{}/'.format(ROOT_DIR, run.upper(), str(fileIndex).zfill(7))  
        filename_plot = '{}ring_current_zeq0_4plot_{}.png'.format(save_dir, str(fileIndex).zfill(5))
        mkdir_path(filename_plot)
        print(filename_plot)
        plt.savefig(filename_plot, dpi = 300)
        plt.close()


            # make standalone plot
        fig, ax3 = plt.subplots()
        im3 = ax3.pcolormesh(x2d[ix1:ix2,iy1:iy2]/R_EARTH, y2d[ix1:ix2,iy1:iy2]/R_EARTH, Jxy_magnitude_plot, shading = 'auto', cmap = 'plasma', vmin = -15, vmax = -8)
        divider = make_axes_locatable(ax3)
        cax3 = divider.append_axes('right', size='5%', pad=0.05)
        cbar3 = fig.colorbar(im3, cax=cax3, orientation='vertical')
        #cbar3.set_label(r'$log10(Jxy [A/m^2])$')
         # now quiver
        q = ax3.quiver(x2d[ix1:ix2:dq,iy1:iy2:dq]/R_EARTH, y2d[ix1:ix2:dq,iy1:iy2:dq]/R_EARTH, Jx_q_plot[0::dq,0::dq], Jy_q_plot[0::dq,0::dq])   # scale = dq*10?
        ax3.set_title(r'$log10(Jxy [A/m^2])$, z=0, run={}, time={}'.format(run, fileIndex) )
        ax3.set_xlim(xlim)
        ax3.set_ylim(ylim)
        ax3.set_xlabel(r'x [$R_E$]')
        ax3.set_ylabel(r'y [$R_E$]')
        ax3.set_aspect('equal', adjustable='box')

        save_dir = '{}{}/ring_current/{}/'.format(ROOT_DIR, run.upper(), str(fileIndex).zfill(7))  
        filename_plot = '{}ring_current_zeq0_{}.png'.format(save_dir, str(fileIndex).zfill(5))
        mkdir_path(filename_plot)
        print(filename_plot)
        plt.savefig(filename_plot, dpi = 300)
        plt.close()



           # Now make same plot in xz plane (y = 0)
           # note I'm reusing and redefining a bunch of variables, bad practice...

        yind_temp = np.where(np.abs(y) == np.nanmin(np.abs(y)))[0][0]
        yeq0_ind = 2
        fg_b_analyze = fg_b[:,yind_temp-yeq0_ind:yind_temp+yeq0_ind+1,:,:]   # reduce the data size so that numcurl3d doesn't require too much memory
        J = (1 / mu_0) * numcurl3d(fg_b_analyze, CELLSIZE_XYZ)
        #J_magnitude = (J[:,:,:,0]**2 + J[:,:,:,1]**2 + J[:,:,:,2]**2)**0.5
        #Jxz_magnitude = (J[:,:,:,0]**2 + J[:,:,:,1]**2 )**0.5
        B = (fg_b_analyze[:,:,:,0]**2 + fg_b_analyze[:,:,:,1]**2 + fg_b_analyze[:,:,:,2]**2)**0.5
        x2_temp, z2_temp = np.meshgrid(x, z, indexing='ij', sparse=True)
        x2d = x2_temp + (z2_temp * 0)
        z2d = z2_temp + (x2_temp * 0)
        r_2d = (x2d**2 + z2d**2)**0.5
        vg_r_min = 4.7
        #Jxz_magnitude[r_2d/R_EARTH < vg_r_min] = 0      # zero out anything outside simulation domain

        fig, ax3 = plt.subplots()
        xlim = [-15, 15]
        zlim = [-15, 15]
        indx = np.where((x/R_EARTH >= xlim[0]) & (x/R_EARTH <= xlim[1]))[0]
        ix1 = np.nanmin(indx); ix2 = np.nanmax(indx)
        minx = np.nanmin(x[indx]); maxx = np.nanmax(x[indx])
        indz = np.where((z/R_EARTH >= zlim[0]) & (z/R_EARTH <= zlim[1]))[0]
        iz1 = np.nanmin(indz); iz2 = np.nanmax(indz)
        minz = np.nanmin(z[indz]); maxz = np.nanmax(z[indz])

           # plot the Jy current in the xz plane
        #Jxz_magnitude_plot = np.log10(Jxz_magnitude[ix1:ix2,yeq0_ind,iz1:iz2])
        #Jx_q_plot = J[ix1:ix2,yeq0_ind,iz1:iz2, 0]/Jxy_magnitude[ix1:ix2,yeq0_ind,iz1:iz2]
        #Jz_q_plot = J[ix1:ix2,yeq0_ind,iz1:iz2, 2]/Jxy_magnitude[ix1:ix2,yeq0_ind,iz1:iz2]
        Jy_plot = J[ix1:ix2,yeq0_ind,iz1:iz2, 1]
        print(r_2d.shape)
        print(Jy_plot.shape)
        print(ix1, ix2, yeq0_ind, iz1, iz2)    # 464 655 2 272 463
        Jy_plot[r_2d[ix1:ix2,iz1:iz2]/R_EARTH < vg_r_min] = 0      # zero out anything outside simulation domain
        #for i in range(dq-1):
        #    Jx_q_plot[i+1::dq,:] = np.nan; Jx_q_plot[:,i+1::dq] = np.nan
        #    Jz_q_plot[i+1::dq,:] = np.nan; Jz_q_plot[:,i+1::dq] = np.nan
        im3 = ax3.pcolormesh(x2d[ix1:ix2,iz1:iz2]/R_EARTH, z2d[ix1:ix2,iz1:iz2]/R_EARTH, Jy_plot, shading = 'auto', cmap = 'bwr', 
                             norm=colors.SymLogNorm(linthresh=2e-15, linscale=0.03, vmin=-1e-7, vmax=1e-7, base =10))
        divider = make_axes_locatable(ax3)
        cax3 = divider.append_axes('right', size='5%', pad=0.05)
        cbar3 = fig.colorbar(im3, cax=cax3, orientation='vertical')
        #cbar3.set_label(r'$log10(Jxz [A/m^2])$')
         # now quiver
        #q = ax3.quiver(x2d[ix1:ix2:dq,iz1:iz2:dq]/R_EARTH, z2d[ix1:ix2:dq,iz1:iz2:dq]/R_EARTH, Jx_q_plot[0::dq,0::dq], Jz_q_plot[0::dq,0::dq])   # scale = dq*10?
        ax3.set_title(r'$Jy [A/m^2]$, y=0, run={}, time={}'.format(run, fileIndex) )
        ax3.set_xlim(xlim)
        ax3.set_ylim(zlim)
        ax3.set_xlabel(r'x [$R_E$]')
        ax3.set_ylabel(r'z [$R_E$]')
        ax3.set_aspect('equal', adjustable='box')

        save_dir = '{}{}/ring_current/{}/'.format(ROOT_DIR, run.upper(), str(fileIndex).zfill(7))
        filename_plot = '{}ring_current_yeq0_{}.png'.format(save_dir, str(fileIndex).zfill(5))
        mkdir_path(filename_plot)
        print(filename_plot)
        plt.savefig(filename_plot, dpi = 300)
        plt.close()






        # Axis limits for VDF plots [in m/s]
        #VDFlim = 3e6
        VDFlim = 4e6
        cutpoint_list = [-2,-1,0,1,2]   #re
        vdf_x_list = [-10,-9,-8,-7]

        for i, x in enumerate(vdf_x_list):
        # Second plot: VDF (slice in the XZ plane; look up the bpara, bpara1, bperp parameters instead of xz in the plot_vdf function)
           cidrequest = f.get_cellid([-x*R_EARTH, 0*R_EARTH, 0])
           cid = plot_vdf.getNearestCellWithVspace(f, cidrequest)
           pt.plot.plot_vdf(filename=filename,cellids=[cid],box=[-VDFlim,VDFlim,-VDFlim,VDFlim],xy=1,fmin=1e-20,fmax=1e-9,axisunit=6,colormap='nipy_spectral',cbulk=1, slicethick=0, outputfile=save_dir+'test_vdf_xy_xeq{}_yeq0_{}_{}.png'.format(x, run,fileIndex,cid))
           pt.plot.plot_vdf(filename=filename,cellids=[cid],box=[-VDFlim,VDFlim,-VDFlim,VDFlim],xz=1,fmin=1e-20,fmax=1e-9,axisunit=6,colormap='nipy_spectral',cbulk=1, slicethick=0, outputfile=save_dir+'test_vdf_xz_xeq{}_yeq0_{}_{}.png'.format(x, run,fileIndex,cid))
           pt.plot.plot_vdf(filename=filename,cellids=[cid],box=[-VDFlim,VDFlim,-VDFlim,VDFlim],yz=1,fmin=1e-20,fmax=1e-9,axisunit=6,colormap='nipy_spectral',cbulk=1, slicethick=0, outputfile=save_dir+'test_vdf_yz_xeq{}_yeq0_{}_{}.png'.format(x, run,fileIndex,cid))

        for i, cutpoint in enumerate(cutpoint_list):
        # Second plot: VDF (slice in the XZ plane; look up the bpara, bpara1, bperp parameters instead of xz in the plot_vdf function)
           x = -9
           cidrequest = f.get_cellid([-x*R_EARTH, cutpoint*R_EARTH, 0])
           cid = plot_vdf.getNearestCellWithVspace(f, cidrequest)
           pt.plot.plot_vdf(filename=filename,cellids=[cid],box=[-VDFlim,VDFlim,-VDFlim,VDFlim],xy=1,fmin=1e-20,fmax=1e-9,axisunit=6,colormap='nipy_spectral',cbulk=1, slicethick=0, outputfile=save_dir+'test_vdf_xy_xeq{}_yeq{}_{}_{}.png'.format(x, cutpoint, run,fileIndex,cid))
           pt.plot.plot_vdf(filename=filename,cellids=[cid],box=[-VDFlim,VDFlim,-VDFlim,VDFlim],xz=1,fmin=1e-20,fmax=1e-9,axisunit=6,colormap='nipy_spectral',cbulk=1, slicethick=0, outputfile=save_dir+'test_vdf_xz_xeq{}_yeq{}_{}_{}.png'.format(x, cutpoint, run,fileIndex,cid))
           pt.plot.plot_vdf(filename=filename,cellids=[cid],box=[-VDFlim,VDFlim,-VDFlim,VDFlim],yz=1,fmin=1e-20,fmax=1e-9,axisunit=6,colormap='nipy_spectral',cbulk=1, slicethick=0, outputfile=save_dir+'test_vdf_yz_xeq{}_yeq{}_{}_{}.png'.format(x, cutpoint, run,fileIndex,cid))

        # NEXT: plot vgy in the y =0 plane, using plot3dslice (or 2d slice?) with field lines overplotted (want to see field line connectivity where ring current is supposed to be)
           dim = '3D'

        # Defining source and output file locations
           #fluxLocation = '/wrk-vakka/group/spacephysics/vlasiator/{}/{}/flux/'.format(dim, run)
           outputLocation = '/wrk-vakka/users/horakons/test_analysator/'
        #operator = 'y'   # (DOESN'T WORK?)y-component. It seems like the data_reducer proton/vg_vy doesn't exist in /proj/horakons/analysator/pyVlsv/reduction.py   ASK ABOUT THIS
        # plot: veleocity vg_v (y?) lines in y=0 plane
           pt.plot.plot_colormap3dslice(filename=filename,var='proton/vg_v', boxre=[-15,15,-15,15],vmin=-2e6,vmax=2e6, lin=1, symlog=0, run=run,colormap='bwr',step=fileIndex,outputdir=save_dir,outputfile='vg_v_yeq{}_{}_{}.png'.format(cutpoint, run,fileIndex),Earth=1, cutpointre = cutpoint)
           pt.plot.plot_colormap3dslice(filename=filename,var='proton/vg_v', operator = 'y', normal = 'y', boxre=[-15,15,-15,15], vmin=-2e5,vmax=2e5, lin=1,symlog=0, run=run,colormap='bwr',step=fileIndex,outputdir=save_dir,outputfile='vg_vy_yeq{}_{}_{}.png'.format(cutpoint, run,fileIndex),Earth=1, streamlines = 'vg_b_vol', cutpointre = cutpoint)
           pt.plot.plot_colormap3dslice(filename=filename,var='proton/vg_v', operator = 'y',normal='z', boxre=[-15,15,-15,15], vmin=-2e5,vmax=2e5, lin=1,symlog=0, run=run,colormap='bwr',step=fileIndex,outputdir=save_dir,outputfile='vg_vy_zeq{}_{}_{}.png'.format(cutpoint, run,fileIndex),Earth=1, streamlines = 'proton/vg_v', cutpointre = cutpoint)
           pt.plot.plot_colormap3dslice(filename=filename,var='proton/vg_v', operator = 'x',normal='z', boxre=[-15,15,-15,15], vmin=-2e5,vmax=2e5, lin=1,symlog=0, run=run,colormap='bwr',step=fileIndex,outputdir=save_dir,outputfile='vg_vx_zeq{}_{}_{}.png'.format(cutpoint, run,fileIndex),Earth=1, streamlines = 'proton/vg_v', cutpointre = cutpoint)
           pt.plot.plot_colormap3dslice(filename=filename,var='proton/vg_v',normal='z', boxre=[-15,15,-15,15], vmin=0,vmax=2e5, lin=1,symlog=0, run=run,colormap='bwr',step=fileIndex,outputdir=save_dir,outputfile='vg_v_zeq{}_{}_{}.png'.format(cutpoint, run,fileIndex),Earth=1, streamlines = 'proton/vg_v', cutpointre=cutpoint)
           #density
           pt.plot.plot_colormap3dslice(filename=filename,var='proton/vg_rho',normal='z', boxre=[-15,15,-15,15], vmin = 1e4, vmax = 1e8, run=run,colormap='plasma',step=fileIndex,outputdir=save_dir,outputfile='vg_rho_zeq{}_{}_{}.png'.format(cutpoint, run,fileIndex),Earth=1, streamlines = 'proton/vg_v', cutpointre=cutpoint)
           pt.plot.plot_colormap3dslice(filename=filename,var='proton/vg_rho', normal = 'y', boxre=[-15,15,-15,15],  vmin = 1e4, vmax = 1e8,run=run,colormap='plasma',step=fileIndex,outputdir=save_dir,outputfile='vg_rho_yeq{}_{}_{}.png'.format(cutpoint, run,fileIndex),Earth=1, streamlines='vg_b_vol', cutpointre = cutpoint)
           #pressure
           pt.plot.plot_colormap3dslice(filename=filename,var='vg_p_perpendicular', normal='z', boxre=[-15,15,-15,15], vmin = 1e-12, vmax = 1e-8, run=run,colormap='plasma',step=fileIndex,outputdir=save_dir,outputfile='vg_p_perp_zeq{}_{}_{}.png'.format(cutpoint, run,fileIndex),Earth=1, streamlines = 'proton/vg_v', cutpointre=cutpoint)
           pt.plot.plot_colormap3dslice(filename=filename,var='vg_p_perpendicular', normal='y', boxre=[-15,15,-15,15], vmin = 1e-12, vmax = 1e-8, run=run,colormap='plasma',step=fileIndex,outputdir=save_dir,outputfile='vg_p_perp_yeq{}_{}_{}.png'.format(cutpoint, run,fileIndex),Earth=1, streamlines='vg_b_vol', cutpointre = cutpoint)






