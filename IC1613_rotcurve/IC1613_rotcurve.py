#Plot IC1613 rotation curve data.
#29.07.21 | Justin Read

if __name__ == "__main__":
    #Import plots library:
    import numpy as np
    import pylab as plt
    from scipy.integrate.quadrature import simps
    from matplotlib import rcParams

    #IC1613 Rhalf & sigstar(Rhalf) [km/s] from Kirby et al. 2014:
    Rhalf_kpc = 1.092
    sigstar = 10.8
    sigstar_err = 1.0
    
    #From Walker et al. 2009:
    Vhalf = sigstar / 0.63
    Vhalferr = sigstar_err / 0.63

    #Other IC1613 parameters (Read et al. 2017):
    dgal_kpc = 740.0
    Mstar = 1.5e7
    Mstar_err = 0.5e7
    Mgas = 8e7
    Mbar = Mstar + Mgas 
    Mbar_err = Mstar_err
    
    #Plot parameters:
    xmin = 0.0
    xmax = 2.5
    ymin = 0.0
    ymax = 40.0
    figsize = 8
    figx = figsize
    figy = figsize
    myfontsize = 25
    mylegendfontsize = 20
    mylinewidth = 3
    mymarkersize = 12

    #Rotation curve data:
    f = open('./Data/1613_RC_CE.txt','r')
    data_in = np.genfromtxt(f)
    f.close()
    data_rot = np.ndarray(shape=(len(data_in[:,0]),4))
    data_rot[:,0] = data_in[:,0]*dgal_kpc/360./60./60.*2.0*np.pi
    data_rot[:,1] = data_in[:,1]
    data_rot[:,2] = data_in[:,2]
    data_rot[:,3] = data_in[:,2]

    #Use LaTeX Fonts:    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    #Set thick axes:
    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(mylinewidth)
        
    #Setup ticks:
    from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
    tick_spacing = 0.01
    ax.minorticks_on()
    ax.tick_params('both', length=15, width=2,\
        which='major',direction='in')
    ax.tick_params('both', length=10, width=1,\
        which='minor',direction='in')
    plt.xticks(fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    
    #Plot the rotation curve data:
    plt.errorbar(data_rot[:,0],data_rot[:,1],\
                 yerr=[np.abs(data_rot[:,2]),np.abs(data_rot[:,3])],\
                 linewidth=2,fmt='o',\
                 color='black',label=r'$i=39^\circ$')

    #Plot the inclination consistent with stellar kinematics:
    inc = 39.0
    inc_true = 20.0
    vc_fac = np.sin(inc/360.0*2.0*np.pi)/\
        np.sin(inc_true/360.0*2.0*np.pi)
    plt.errorbar(data_rot[:,0],data_rot[:,1]*vc_fac,\
                 yerr=[np.abs(data_rot[:,2]),\
                    np.abs(data_rot[:,3])],\
                 linewidth=2,fmt='o',alpha=0.5,\
                 color='black',label=r'$i=20^\circ$')

    #Overlay location of the HI holes:
    plotholes = 'yes'
    inc_sigstar = 'yes'
    if (plotholes == 'yes'):
        holeloclow = np.array([0.,0.57,1.32,0.66])
        holelochigh = np.array([0.85,1.27,1.60,1.13])
        holeloc = (holelochigh - holeloclow)/2. + holeloclow
        holelocerr = (holelochigh - holeloclow)/2.0
        yholepos = 1.5
    for i in range(len(holeloc)):
        yholep = yholepos+i*0.4
        if (i == len(holeloc)-1):
            yholep = yholepos+(i-1)*0.4
        plt.errorbar(holeloc[i],yholep,\
            xerr=holelocerr[i],\
            linewidth=2,color='black')
    plt.annotate('HI holes',xy=(0.75,yholepos),\
             xytext=(0.75,yholepos+(len(holeloc)-1)*0.4),\
                fontsize=mylegendfontsize,\
                color='black')

    if (inc_sigstar == 'yes'):
        plt.errorbar([Rhalf_kpc],[Vhalf],[Vhalferr],fmt='o',\
                     color='blue',markersize=mymarkersize,\
                     linewidth=mylinewidth,\
                     label='Stellar kinematics')
  
    #Overlay Vmax expectation from McGaugh 2012 BTFR
    A = 47.0
    A_err = 6.0
    Vmax_pred = (Mbar/A)**(0.25)
    Vmax_pred_err = Vmax_pred*np.sqrt((0.25*Mbar_err/Mbar)**2.0+\
        (0.25*A_err/A)**2.0)
    pnts = 20
    xarr = np.linspace(xmin,xmax,pnts)
    yarr_min = np.linspace(Vmax_pred-Vmax_pred_err,\
        Vmax_pred-Vmax_pred_err,pnts)
    yarr_max = np.linspace(Vmax_pred+Vmax_pred_err,\
        Vmax_pred+Vmax_pred_err,pnts)        
    plt.fill_between(xarr,yarr_max,yarr_min,\
                 facecolor='black',edgecolor='none',alpha=0.25,\
                 zorder=0,label='BTFR')
    
    #Axis labels:
    plt.xlabel(r'Radius\,(kpc)',\
               fontsize=myfontsize)
    plt.ylabel(r'Circular speed\,(km\,s$^{-1}$)',\
               fontsize=myfontsize)

    #Plot range:
    plt.ylim([ymin+0.0001,ymax])
    plt.xlim([xmin,xmax])
    
    plt.legend(loc='upper left',fontsize=mylegendfontsize,\
        facecolor='white', framealpha=1)
    plt.savefig('IC1613_rotcurve.pdf',bbox_inches='tight')
