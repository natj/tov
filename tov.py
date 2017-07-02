import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import units as cgs
from math import pi
from polytropes import monotrope, polytrope
from crust import SLyCrust
from eoslib import get_eos, glue_crust_and_core, eosLib
from scipy.integrate import odeint
from label_line import label_line

import math

#from matplotlib import cm
import palettable as pal
cmap = pal.colorbrewer.qualitative.Set1_6.mpl_colormap
#cmap = pal.cmocean.sequential.Matter_8.mpl_colormap #best so far
#cmap = pal.wesanderson.Zissou_5.mpl_colormap

#--------------------------------------------------

class tov:

    def __init__(self, peos):
        self.physical_eos = peos

    def tov(self, y, r):
        P, m = y
        eden = self.physical_eos.edens_inv( P )

        dPdr = -cgs.G*(eden + P/cgs.c**2)*(m + 4.0*pi*r**3*P/cgs.c**2)
        dPdr = dPdr/(r*(r - 2.0*cgs.G*m/cgs.c**2))
        dmdr = 4.0*pi*r**2*eden

        return [dPdr, dmdr]

    def tovsolve(self, rhoc):

        N = 800
        r = np.linspace(1e0, 18e5, N)
        P = self.physical_eos.pressure( rhoc )
        eden = self.physical_eos.edens_inv( P )
        m = 4.0*pi*r[0]**3*eden

        psol = odeint(self.tov, [P, m], r, rtol=1.0e-4, atol=1.0e-4)
        return r, psol[:,0], psol[:,1]

    def mass_radius(self):
        N = 100
        mcurve = np.zeros(N)
        rcurve = np.zeros(N)
        rhocs = np.logspace(14.0, 16.0, N)

        mass_max = 0.0
        j = 0
        for rhoc in rhocs:
            rad, press, mass = self.tovsolve(rhoc)

            rad  /= 1.0e5 #cm to km
            mass /= cgs.Msun

            mstar = mass[-1]
            rstar = rad[-1]
            for i, p in enumerate(press):
                if p > 0.0:
                    mstar = mass[i]
                    rstar = rad[i]
            mcurve[j] = mstar
            rcurve[j] = rstar

            j += 1
            if mass_max < mstar:
                mass_max = mstar
            else:
                break

        return mcurve[:j], rcurve[:j], rhocs[:j]


#--------------------------------------------------
def main(argv):

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=7)
    

    fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
    #fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
    gs = plt.GridSpec(1, 1)

    ax = plt.subplot(gs[0, 0])
    ax.minorticks_on()
    ax.set_xlim(9.0, 16.0)
    ax.set_ylim(0.0, 3.0)

    ax.set_xlabel(r'Radius $R$ (km)')
    ax.set_ylabel(r'Mass $M$ (M$_{\odot}$)')


    # test single eos separately
    if False: 
        key = 'SLy'
        #key = 'BGN1H1'
        #key = 'ALF2'
        #key = 'ENG'
        #key = 'PAL6'

        dense_eos = get_eos(key)
        eos = glue_crust_and_core( SLyCrust, dense_eos )
        t = tov(eos)

        mass, rad, rho = t.mass_radius()
        print mass
        print rad
        ax.plot(rad, mass)


    if True:
        i = 0
        for key, value in eosLib.iteritems():

            print "Solving TOVs for:", key, "(",i, ")"

            dense_eos = get_eos(key)
            eos = glue_crust_and_core( SLyCrust, dense_eos )
            t = tov(eos)
            mass, rad, rhoc = t.mass_radius()
            
            linestyle='solid'
            col = 'k'
            if value[4] == 'npem':
                col = 'k'
            if value[4] == 'meson':
                col = 'b'
            if value[4] == 'hyperon':
                col = 'g'
            if value[4] == 'quark':
                col = 'r'

            l, = ax.plot(rad, mass, color=col, linestyle=linestyle, alpha = 0.9)


            # labels for lines
            near_y = 1.45
            near_x = None
            rotation_offset=180.0
            if key == 'APR3':
                near_x = 11.0
                near_y = None
            if key == 'ENG':
                near_x = 11.2
                near_y = None
            if key == 'ALF2':
                near_y = 1.0
                rotation_offset = 0.0
            if key == 'MPA1':
                near_x = 12.0
                near_y = None
            if key == 'MS1b':
                near_y = 0.4
            if key == 'MS1':
                near_y = 2.3

            label_line(l, key, near_y=near_y, near_x=near_x, rotation_offset=rotation_offset)



            i += 1


    #plot approximate central pressure isocurves
    #twice nuclear saturation density
    if False:
        x = [11.0, 15.8]
        y = [0.4, 2.5]
        ax.plot(x, y, "r--")

        txt = ax.text(15.5, 2.35, r'$2 \rho_n$', rotation=32, ha='center', va='center', size=8)
        txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))




if __name__ == "__main__":
    main(sys.argv)
    plt.subplots_adjust(left=0.15, bottom=0.16, right=0.98, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig('mr.pdf')



