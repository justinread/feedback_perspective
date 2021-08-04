#Plot data on burstiness of star formation 
#in dwarf galaxies.
#29.07.21 | Justin Read

def cosmo_cfunc(M200,h):
    #This is c200 from Dutton & Maccio 2014:
    c = 10.0**(0.905 - 0.101*np.log10(M200/(1e12/h)))
    return c

def Rhalf_func(M200,rhocrit):
    r200 = (3./4.*M200/(np.pi*200.0*rhocrit))**(1./3.)
    
    #Rhalf from Kravtsov 2014:
    Rhalf = 0.015 * r200
    return Rhalf
     
def tdyn_func(M200,rhocrit,h):
    G = 6.67e-11
    kpc = 3.086e19
    Msun = 1.989e30
    Gyr = 365.*24.*60.*60.*1e9

    #Calculate r200 and Rhalf:
    r200 = (3./4.*M200/(np.pi*200.0*rhocrit))**(1./3.)*kpc
    Rhalf = Rhalf_func(M200,rhocrit)*kpc

    #And tdyn at Rhalf assuming NFW profile
    #with negligible baryonic contribution:
    c = cosmo_cfunc(M200,h)
    gcon=1./(np.log(1.+c)-c/(1.+c))
    deltachar=200.0*c**3.*gcon/3.
    rs=r200/c
    manal_Rhalf = M200 * gcon * Msun * \
        (np.log(1.0 + Rhalf/rs)-Rhalf/rs/(1.0+Rhalf/rs))
    tdyn = 2.0*np.pi*np.sqrt(Rhalf**3.0/(G*manal_Rhalf))/Gyr
    
    return tdyn

if __name__ == "__main__":
    #Import plots library:
    import numpy as np
    import pylab as plt
    from scipy.integrate.quadrature import simps
    from matplotlib import rcParams
    
    #Plot parameters:
    xmin = 7.8
    xmax = 10.5
    ymin = 0.0
    ymax = 1.0
    xmin2 = 6
    xmax2 = 9
    ymin2 = 0.0
    ymax2 = 320
    xmin3 = 6
    xmax3 = 9
    ymin3 = 0.0
    ymax3 = 110.0

    figsize = 8
    figx = figsize*3
    figy = figsize
    myfontsize = 25
    mylegendfontsize = 20
    mylinewidth = 3
    mymarkersize = 12

    #Load in the Kauffmann data:
    f = open('./Data/Kauffmann2014_cont_SFH.txt','r')
    data_in = np.genfromtxt(f)
    f.close()
    data_cont = np.ndarray(shape=(len(data_in[:,0]),2))
    data_cont[:,0] = data_in[:,0]
    data_cont[:,1] = data_in[:,1]

    f = open('./Data/Kauffmann2014_burst_SFH.txt','r')
    data_in = np.genfromtxt(f)
    f.close()
    data_burst = np.ndarray(shape=(len(data_in[:,0]),2))
    data_burst[:,0] = data_in[:,0]
    data_burst[:,1] = data_in[:,1]
    
    #Use LaTeX Fonts:    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=myfontsize)

    #Set up subplots:
    fig, (ax0, ax1, ax2) = \
        plt.subplots(nrows=1, ncols=3,
           figsize=(figx, figy))   

    ### Burst fraction plot ###

    #Set thick axes:
    for axis in ['top','bottom','left','right']:
        ax0.spines[axis].set_linewidth(mylinewidth)
       
    #Setup ticks:
    from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
    tick_spacing = 0.01
    ax0.minorticks_on()
    ax0.tick_params('both', length=15, width=2,\
        which='major',direction='in')
    ax0.tick_params('both', length=10, width=1,\
        which='minor',direction='in')
    
    #Plot the data:
    ax0.plot(data_cont[:,0],data_cont[:,1],\
                 linewidth=mylinewidth,\
                 color='red',label=r'Kauffmann 2014 continuous SF')
    ax0.plot(data_burst[:,0],data_burst[:,1],\
                 linewidth=mylinewidth,\
                 color='blue',label=r'Kauffmann 2014 bursty SF')
    
    #Axis labels:
    ax0.set_xlabel(r'Log$_{10}$[Stellar mass\,(M$_\odot$)]',\
               fontsize=myfontsize)
    ax0.set_ylabel(r'Fraction',\
               fontsize=myfontsize)

    #Plot range:
    ax0.set_ylim([ymin,ymax])
    ax0.set_xlim([xmin,xmax])
    
    ax0.legend(loc='upper right',\
        fontsize=mylegendfontsize,\
        facecolor='white', framealpha=1)
 
    ### Burst duration plot ### 
    
    #Set thick axes:
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(mylinewidth)
       
    #Setup ticks:
    ax1.minorticks_on()
    ax1.tick_params('both', length=15, width=2,\
        which='major',direction='in')
    ax1.tick_params('both', length=10, width=1,\
        which='minor',direction='in')
        
    #Plot Kauffmann 2014 calculation:
    ax1.plot(np.log10(1.0e8),2.0e8/1.0e6,'*',color='blue',\
        markersize=mymarkersize*2,label='Kauffmann 2014')
        
    #Plot Emami et al. 2019 (builds on Weisz et al. 2012):
    f = open('./Data/Emami2019_burst.txt','r')
    data_E = np.genfromtxt(f)
    f.close()
    LogMstarE = (data_E[:,0]+data_E[:,1])/2.0
    LogMstarE_err = LogMstarE-data_E[:,0]
    PeriodE = (data_E[:,2]+data_E[:,3])/2.0
    PeriodE_err = PeriodE - data_E[:,2]
    AmpE = data_E[:,4]
    ax1.errorbar(LogMstarE,PeriodE,\
        xerr=LogMstarE_err,yerr=PeriodE_err,\
        markersize=mymarkersize,\
        fmt='o',color='green',label='Emami et al. 2019')

    #Overplot an estimate of the mean dynamical time
    #expected in LCDM at ~Rhalf. Uses M*-M200 relation
    #from Read et al. 2017, the Rhalf-r200 relation
    #from Kravtsov 2014 and the M200-c200 relation
    #from Dutton & Maccio 2014.
    h = 0.7          #Dimensionless Hubble 
    rhocrit = 136.0  #Msun kpc^-3
    f = open('./Data/mstar-mhalo-field.txt','r')
    data = np.genfromtxt(f)
    f.close()
    
    Mstar = np.logspace(xmin2,xmax2,250)
    M200lowarr_MS = data[:,0][np.where(data[:,1])]
    Mstarlow = data[:,1][np.where(data[:,1])]
    M200higharr_MS = data[:,0][np.where(data[:,2])]
    Mstarhigh = data[:,2][np.where(data[:,2])]
    M200_Mstarhi = np.interp(Mstar,Mstarlow,M200lowarr_MS)
    M200_Mstarlo = np.interp(Mstar,Mstarhigh,M200higharr_MS)
    M200_Mstar = (M200_Mstarhi+M200_Mstarlo)/2.0
    tdyn_lo = tdyn_func(M200_Mstarlo,rhocrit,h) * 1000.0
    tdyn_hi = tdyn_func(M200_Mstarhi,rhocrit,h) * 1000.0

    ax1.fill_between(np.log10(Mstar),tdyn_lo,\
        tdyn_hi,facecolor='black',edgecolor='none',alpha=0.5,\
        label='Dynamical time at $R_{1/2}$')
    ax1.plot(np.log10(Mstar),tdyn_lo,linewidth=2,\
        color='black')
    ax1.plot(np.log10(Mstar),tdyn_hi,linewidth=2,\
        color='black')
        
    #Mark impulsive versus adiabatic transition:
    ax1.annotate('Impulsive',xy=(6.1,135),\
                rotation=9,\
                fontsize=mylegendfontsize,\
                color='black')
    ax1.annotate('Adiabatic',xy=(6.1,157.5),\
                rotation=9,\
                fontsize=mylegendfontsize,\
                color='black')
 
    #Axis labels:
    ax1.set_xlabel(r'Log$_{10}$[Stellar mass\,(M$_\odot$)]',\
               fontsize=myfontsize)
    ax1.set_ylabel(r'Burst timescale (Myrs)',\
               fontsize=myfontsize)

    #Plot range:
    ax1.set_ylim([ymin2,ymax2])
    ax1.set_xlim([xmin2,xmax2])
    
    #Legend:
    handles, labels = ax1.get_legend_handles_labels()
    order = [0,2,1]
    ax1.legend([handles[idx] for idx in order],\
        [labels[idx] for idx in order],\
        loc='upper left',\
        fontsize=mylegendfontsize,\
        facecolor='white', framealpha=1)
 
    ### Burst amplitude plot ### 
    
    #Set thick axes:
    for axis in ['top','bottom','left','right']:
        ax2.spines[axis].set_linewidth(mylinewidth)
       
    #Setup ticks:
    ax2.minorticks_on()
    ax2.tick_params('both', length=15, width=2,\
        which='major',direction='in')
    ax2.tick_params('both', length=10, width=1,\
        which='minor',direction='in')
      
    #Plot Kauffmann 2014 calculation:
    ax2.plot(np.log10(1.0e8),25.0,'*',color='blue',\
        markersize=mymarkersize*2,label='Kauffmann 2014')
    ax2.plot(np.log10(1.0e10),10.0,'*',color='blue',\
        markersize=mymarkersize*2)
    
    #Plot Emami et al. 2019:
    ax2.errorbar(LogMstarE,AmpE,\
        xerr=LogMstarE_err,\
        markersize=mymarkersize,\
        fmt='o',color='green',label='Emami et al. 2019')
  
    #Axis labels:
    ax2.set_xlabel(r'Log$_{10}$[Stellar mass\,(M$_\odot$)]',\
               fontsize=myfontsize)
    ax2.set_ylabel(r'Burst amplitude',\
               fontsize=myfontsize)

    #Plot range:
    ax2.set_ylim([ymin3,ymax3])
    ax2.set_xlim([xmin3,xmax3])
    
#    ax2.legend(loc='upper right',\
#        fontsize=mylegendfontsize,\
#        facecolor='white', framealpha=1)
 
    #Save the figure: 
    plt.savefig('Bursty_SF.pdf',bbox_inches='tight')
