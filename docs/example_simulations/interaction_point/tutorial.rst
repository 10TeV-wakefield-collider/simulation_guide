.. raw:: html

   <h2 style="color: #1a202c; margin: 1.5em 0 0.5em 0;">Visualizing the beams with Guinea-Pig and WarpX</h2>


.. code:: ipython3

    import os, sys
    import numpy as np 
    import matplotlib.pyplot as plt
    from openpmd_viewer import OpenPMDTimeSeries
    from scipy.constants import c, micro, nano, pi

.. raw:: html

   <h3 style="color: #2d3748; margin: 1.2em 0 0.5em 0;">Collider parameters</h3>


.. code:: ipython3

    sigmaz = 100*micro
    sigmax = 210*nano
    sigmay = 3.1*nano
    npart = 6.24e9
    n0 = npart / (sigmax * sigmay * sigmaz * (2.*pi)**(3./2.))

.. raw:: html

   <h3 style="color: #2d3748; margin: 1.2em 0 0.5em 0;">Simulation parameters</h3>


.. code:: ipython3

    # time
    n_iterations = 256
    iterations = range(n_iterations)
    # box
    nx = 512
    ny = 512
    nz = 512
    Lx = 20*sigmax
    Ly = 20*sigmay
    Lz = 16*sigmaz 
    gridx = np.linspace(-0.5*Lx, 0.5*Lx, nx+1)
    gridy = np.linspace(-0.5*Ly, 0.5*Ly, ny+1)
    gridz = np.linspace(-0.5*Lz, 0.5*Lz, nz+1)
    dx = gridx[1]-gridx[0]
    dy = gridy[1]-gridy[0]
    dz = gridz[1]-gridz[0]
    # particles
    nmacropart = 1e5
    w0 = npart / nmacropart # weight

.. raw:: html

   <h3 style="color: #2d3748; margin: 1.2em 0 0.5em 0;">Function to plot a single Guine-Pig step</h3>


.. code:: ipython3

    def one_step_gp(n, gp_dir):
        '''
        inputs
            n: current timestep
            gp_dir: simulation folder 
        outputs
            H_zx, H_zy: density of the beams integrated along y and x, resp
        '''
        global w0, dx, dy, dz, gridx, gridy, gridz
        
        # get beams' data, columns are:
        # Particle Energy [GeV] | x [um] | y [um] | z [um] | x' [urad] | y' [urad] 
        data1 = np.loadtxt(os.path.join(gp_dir, 'b1.%d' % n))
        data2 = np.loadtxt(os.path.join(gp_dir, 'b2.%d' % n))
    
        # get the number of macroparticles in the beams
        N1 = np.shape(data1)[0]
        N2 = np.shape(data2)[0]
    
        # stack the data together and convert to SI
        x_data = np.hstack((data1[:,1], data2[:,1]))*micro
        y_data = np.hstack((data1[:,2], data2[:,2]))*micro
        z_data = np.hstack((data1[:,3], data2[:,3]))*micro
    
        weights = np.ones(N1+N2) * w0
        weights[-N2:] = - w0 # assign negative weights to the second beam
    
        H_zx, bx, bz = np.histogram2d(x_data, z_data, bins=(gridx, gridz), weights=weights)
        H_zy, by, bz = np.histogram2d(y_data, z_data, bins=(gridy, gridz), weights=weights)
    
        print(np.max(np.abs(H_zx)), np.max(np.abs(H_zy)))
    
        return H_zx/(dz*dx), H_zy/(dz*dy)


.. raw:: html

   <h3 style="color: #2d3748; margin: 1.2em 0 0.5em 0;">Function to plot a single WarpX step</h3>


.. code:: ipython3

    def one_step_wx(n, wx_series):
        '''
        inputs
            n: current timestep
            wx_dir: simulation folder 
        outputs
            H_zx: density of the beams projected on the plane (z,x), integrated along y
            H_zy: density of the beams projected on the plane (z,y), integrated along x
        '''    
        global dx, dy, dz, gridx, gridy, gridz
    
        x1,y1,z1,w1 = wx_series.get_particle(['x','y','z','w'], species='beam1', iteration=n)
        x2,y2,z2,w2 = wx_series.get_particle(['x','y','z','w'], species='beam2', iteration=n)
        w1 = -w1
        X = np.hstack((x1,x2))
        Y = np.hstack((y1,y2))
        Z = np.hstack((z1,z2))
        W = np.hstack((w1,w2))
    
        H_zx, bx, bz = np.histogram2d(X, Z, bins=(gridx, gridz), weights=W)
        H_zy, by, bz = np.histogram2d(Y, Z, bins=(gridy, gridz), weights=W)
    
        print(np.max(np.abs(H_zx)), np.max(np.abs(H_zy)))
    
        return H_zx/(dz*dx), H_zy/(dz*dy)


.. raw:: html

   <h3 style="color: #2d3748; margin: 1.2em 0 0.5em 0;">Generate a video of the colliding beams</h3>


.. code:: ipython3

    !mkdir -p "plots"
    
    gp_dir = "gp"
    wx_dir = "wx"
    
    path=os.path.join(wx_dir, 'diags/trajs')
    series = OpenPMDTimeSeries(path)
    
    v0 = n0 * np.sqrt(sigmax * sigmay) * 0.05
    extent_zx = [gridz[0]/micro, gridz[-1]/micro, gridx[0]/nano, gridx[-1]/nano]
    extent_zy = [gridz[0]/micro, gridz[-1]/micro, gridy[0]/nano, gridy[-1]/nano]
    
    plt.rcParams.update({'font.size': 16})
    
    # loop through the timesteps
    for n in iterations:
    
        # prepare canvas
        fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(16,12), dpi=300, sharex='col', sharey='row')
        
        H_zx, H_zy = one_step_gp(n, gp_dir)
        im=ax[0][0].imshow(H_zx, extent=extent_zx, cmap='seismic', origin='lower', interpolation='nearest', vmin=-v0, vmax=v0)
        im=ax[1][0].imshow(H_zy, extent=extent_zy, cmap='seismic', origin='lower', interpolation='nearest', vmin=-v0, vmax=v0)
        ax[0][0].set_title("Guinea-Pig")
    
        
        H_zx, H_zy = one_step_wx(n, series)
        im=ax[0][1].imshow(H_zx, extent=extent_zx, cmap='seismic', origin='lower', interpolation='nearest', vmin=-v0, vmax=v0)
        im=ax[1][1].imshow(H_zy, extent=extent_zy, cmap='seismic', origin='lower', interpolation='nearest', vmin=-v0, vmax=v0)
        ax[0][1].set_title("WarpX")
    
        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.88, 0.085, 0.015, 0.87])
        fig.colorbar(im, cax=cbar_ax, label='density [arb. units]')
    
        for a in ax.reshape(-1):
            a.set_aspect('auto')
        
        ax[1][0].set_xlabel(r'z [$\mu$m]')
        ax[1][1].set_xlabel(r'z [$\mu$m]')
        
        ax[0][0].set_ylabel(r'x [nm]')
        ax[1][0].set_ylabel(r'y [nm]')
    
        plt.savefig(f"plots/img_{n:04d}",dpi=300, bbox_inches='tight') 
        plt.close("all")

.. code:: ipython3

    ! ffmpeg -framerate 30 -i 'plots/img_%04d.png' -y out.mp4

.. code:: ipython3

    from IPython.display import Video, display
    display(Video("out.mp4", embed=True))
