import numpy as np
import matplotlib.pyplot as plt

# Plot style

plt.rcParams["font.family"] = "Times New Roman"
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)

fsize = 18
tsize = 15
tdir = 'in'
major = 5.0
minor = 3.0
lwidth = 1.8
lhandle = 2.0

plt.style.use('default')
plt.rcParams['font.size'] = fsize
plt.rcParams['legend.fontsize'] = tsize
plt.rcParams['text.usetex']=True
plt.rcParams['xtick.minor.visible'], plt.rcParams['xtick.top'] = True,True
plt.rcParams['ytick.minor.visible'], plt.rcParams['ytick.right'] = True,True
plt.rcParams['xtick.direction'], plt.rcParams['ytick.direction'] = 'in','in' 

left = b1 = 0.15
width = 0.75
h1 = 0.75
height = 0.75

# Read in data

name, mv, mvu, mvd, feh, fe, tag = np.loadtxt('Mass_metallicity.cat', unpack=True, dtype ='U8,f,f,f,f,f,U8')

mw = tag == 'MW'
m31 = tag == 'M31'
iso = tag == 'Iso'
can = tag == 'Cen'

# Plot data

fig=plt.figure(figsize=(8,6))
ax=fig.add_axes([left,b1,width,height])
ax.set_ylim(-3.2, -0.4)
ax.set_xlim(-0.1,-17.2)
ax.set_ylabel('[Fe/H]', labelpad = 10)
ax.set_xlabel('$M_V$', labelpad = 10,)

ax.errorbar(mv, feh, yerr = fe, xerr =[mvu,mvd] , fmt = 'None', ecolor = '0.5',elinewidth = 0.5, zorder=0)

ax.scatter(mv[iso],feh[iso],c='tomato',marker='*',s=150,label='LG isolated')
ax.scatter(mv[mw],feh[mw],c='mediumorchid',s=80, zorder=80,label='Milky Way')
ax.scatter(mv[m31],feh[m31],c='turquoise',marker='p',s=70)
ax.scatter(mv[can],feh[can],c='gold',marker='s',s=50,zorder=10,label='Cen A')

# Sgr tides
ax.arrow(-15.5,-0.53, 1.8,0, color='k', head_width=0.05,head_length=0.15,length_includes_head=True)
ax.scatter(-15.5,-0.53,c='darkorchid',marker='x',s=80)

# Kirby+2013 relation

x = np.arange(1e1,5e9,1e4)
y = -1.68+0.29*np.log10(x/1e6)
y1 = y+0.16
y2 = y-0.16

xmv = -2.5*np.log10(x)+4.83
plt.plot(xmv,y,ls='--',c='k',zorder=0,label='Kirby+13 relation')

fill_kwargs = {"zorder":0}
plt.fill_between(xmv,y1,y2,color='0.8',alpha=0.5,**fill_kwargs)

ax.legend(loc='upper left')

plt.savefig('Mv-FeH.pdf')


# Plus sims

fs = np.genfromtxt('Simulations.dat',dtype='U24,f,f',
                  names='id,mv,feh',
                  usecols=(0,1,2))

id_edge=fs['id'][0:17]
mv_edge=fs['mv'][0:17]
feh_edge=fs['feh'][0:17]

id_fire=fs['id'][17:21]
mv_fire=-2.5*np.log10(fs['mv'][17:21])+4.83 #M_* assuming M/L=1
feh_fire=fs['feh'][17:21]
                      
id_justice=fs['id'][21:]
mv_justice=fs['mv'][21:]
feh_justice=fs['feh'][21:]

ax.scatter(mv_edge,feh_edge,marker='o',edgecolor='k',lw=2,color='None',s=90,zorder=100,label='EDGE')
ax.scatter(mv_fire,feh_fire,marker='^',edgecolor='r',lw=2,color='None',s=90,zorder=100,label='FIRE')
ax.scatter(mv_justice,feh_justice,marker='s',edgecolor='b',lw=2,color='None',s=90,zorder=100,label='Justice League')

plt.savefig('Mv-FeH-sims.pdf')


