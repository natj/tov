import units as cgs
from math import pi

from polytropes import monotrope
from polytropes import polytrope

#import numpy as np

#Dictionary from Read et al 2009 
# all M_max < 2Msun commented out
eosLib = {
#    'PAL6'  :[ 34.380,  2.227,  2.189,  2.159, 'npem' ],
    'SLy'   :[ 34.384,  3.005,  2.988,  2.851, 'npem' ],
#    'APR1'  :[ 33.943,  2.442,  3.256,  2.908, 'npem' ],
#    'APR2'  :[ 34.126,  2.643,  3.014,  2.945, 'npem' ],
    'APR3'  :[ 34.392,  3.166,  3.573,  3.281, 'npem' ],
    'APR4'  :[ 34.269,  2.830,  3.445,  3.348, 'npem' ],
#    'FPS'   :[ 34.283,  2.985,  2.863,  2.600, 'npem' ],
    'WFF1'  :[ 34.031,  2.519,  3.791,  3.660, 'npem' ],
    'WFF2'  :[ 34.233,  2.888,  3.475,  3.517, 'npem' ],
#    'WFF3'  :[ 34.283,  3.329,  2.952,  2.589, 'npem' ],
#    'BBB2'  :[ 34.331,  3.418,  2.835,  2.832, 'npem' ],
#    'BPAL12':[ 34.358,  2.209,  2.201,  2.176, 'npem' ],
    'ENG'   :[ 34.437,  3.514,  3.130,  3.168, 'npem' ],
    'MPA1'  :[ 34.495,  3.446,  3.572,  2.887, 'npem' ],
    'MS1'   :[ 34.858,  3.224,  3.033,  1.325, 'npem' ],
#    'MS2'   :[ 34.605,  2.447,  2.184,  1.855, 'npem' ],
    'MS1b'  :[ 34.855,  3.456,  3.011,  1.425, 'npem' ],
#    'PS'    :[ 34.671,  2.216,  1.640,  2.365, 'meson' ],
#    'GS1a'  :[ 34.504,  2.350,  1.267,  2.421, 'meson' ],
#    'GS2a'  :[ 34.642,  2.519,  1.571,  2.314, 'meson' ],
#    'BGN1H1':[ 34.623,  3.258,  1.472,  2.464, 'hyperon' ],
#    'GNH3'  :[ 34.648,  2.664,  2.194,  2.304, 'hyperon' ],
#    'H1'    :[ 34.564,  2.595,  1.845,  1.897, 'hyperon' ],
#    'H2'    :[ 34.617,  2.775,  1.855,  1.858, 'hyperon' ],
#    'H3'    :[ 34.646,  2.787,  1.951,  1.901, 'hyperon' ],
    'H4'    :[ 34.669,  2.909,  2.246,  2.144, 'hyperon' ],
#    'H5'    :[ 34.609,  2.793,  1.974,  1.915, 'hyperon' ],
#    'H6a'   :[ 34.593,  2.637,  2.121,  2.064, 'hyperon' ],
#    'H7'    :[ 34.559,  2.621,  2.048,  2.006, 'hyperon' ],
#    'PCL2'  :[ 34.507,  2.554,  1.880,  1.977, 'hyperon' ],
#    'ALF1'  :[ 34.055,  2.013,  3.389,  2.033, 'quark' ],
    'ALF2'  :[ 34.616,  4.070,  2.411,  1.890, 'quark' ],
#    'ALF3'  :[ 34.283,  2.883,  2.653,  1.952, 'quark' ],
#    'ALF4'  :[ 34.314,  3.009,  3.438,  1.803, 'quark' ]
}


#EoS class using dense eos from Read et al (2009) 
def get_eos(key):

    #read eos table for parameters
    ll = eosLib[ key ]

    #dense eos starting pressure
    p1 = 10**ll[0]

    #polytrope indices
    g1 = ll[1] 
    g2 = ll[2]
    g3 = ll[3]

    #transition densities
    r1 = 2.8e14
    r2 = 10**14.7
    r3 = 10**15.0

    #scaling constants
    K1 = p1/(r2**g1)
    K2 = K1 * r2**(g1-g2)
    K3 = K2 * r3**(g2-g3)

    tropes = [monotrope(K1, g1),
              monotrope(K2, g2),
              monotrope(K3, g3) ]
    trans = [r1, r2, r3]

    dense_eos = polytrope(tropes, trans)

    return dense_eos

# Smoothly glue core to SLy crust
# for polytropic eos we can just unpack 
# and repack the piecewise presentation
def glue_crust_and_core(crust, core):

    #unpack crust and core
    tropes_crust = crust.tropes
    trans_crust  = crust.transitions

    tropes_core = core.tropes
    trans_core  = core.transitions

    #find transition depth
    rho_tr = (tropes_core[0].K / tropes_crust[-1].K )**( 1.0/( tropes_crust[-1].G - tropes_core[0].G ) )
    #print "Transition from core to crust at", rho_tr, np.log10(rho_tr), crust.edens_inv( crust.pressure( rho_tr ) )/cgs.GeVfm_per_dynecm
    trans_core[0] = rho_tr
    #trans_crust[-1] = rho_tr

    #repack
    tropes = tropes_crust + tropes_core
    trans  = trans_crust  + trans_core

    for trope in tropes:
        trope.a = 0.0

    eos = polytrope( tropes, trans )
    return eos






##################################################
# Some simple phenomenological mono/bitrope EoSs
def simple_eos():

    K1 = 3.99873692e-8 
    G1 = 1.35692395 
    r1 = 0

    K2 = 2.23872092e-10
    G2 = 3 
    #a = 0.010350691 * c *c
    r2 = 1.4172900e14

    m1 = monotrope(K1*(cgs.c**2), G1)
    m2 = monotrope(K2, G2)
    eos = polytrope([m1, m2], [r1, r2])

    return eos

def simple_eos2():
    Gamma0 = 5.0/3.0 
    K0 = (3.0*pi**2)**(2.0/3.0)*cgs.hbar**2/(5.0*cgs.mn**(8.0/3.0))
    Gamma1 = 3.0 # high densities: stiffer equation of state
    #Gamma1 = 2.5 # high densities: softer equation of state
    rho1 = 5e14
    P1 = K0*rho1**Gamma0
    K1 = P1/rho1**Gamma1

    m1 = monotrope(K0, Gamma0)
    m2 = monotrope(K1, Gamma1)
    eos = polytrope([m1, m2], [0.0, rho1])

    return eos

def simpler_eos():
    K = 1.982e-6
    G = 2.75
    m = monotrope(K, G)
    eos = polytrope([m], [0.0])
    return eos


# If run as a script, lets visualize the library
if __name__ == "__main__":
    import matplotlib
    import matplotlib.pyplot as plt
    from crust import SLyCrust
    import numpy as np
    from label_line import label_line


    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=7)

    fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
    #fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
    gs = plt.GridSpec(1, 1)

    ax = plt.subplot(gs[0, 0])
    ax.minorticks_on()
    ax.set_xlim(2.0e13, 2.0e16)
    ax.set_ylim(1.0e32, 1.0e39)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(r'Density $\rho$ (g cm$^{-3}$)')
    ax.set_ylabel(r'Pressure $P$ (dyne cm$^{-2}$)')


    #rho_ND
    rhond = 4.0e11

    #normal nuclear saturation density
    rhon = 2.8e14

    #inner crust
    col = 'lightgrey'
    xmin = rhond
    xmax = 0.5*rhon
    ax.fill_between([xmin, xmax], [1e10, 1e10], [1e42, 1e42], facecolor=col, color=None, alpha=1.0, edgecolor=col)

    #different ns structures
    y_text = 1.0e38
    ax.text(5.0e13, y_text, 'Inner\ncrust', rotation=0, ha='center', va='center', size=8)
    ax.text(6.0e14, y_text, 'Core', rotation=0, ha='center', va='center', size=8)

    
    lstyle = 'dotted'
    ax.plot([rhon, rhon], [1.0e16, 1e40], "r", linestyle=lstyle)
    txt = ax.text(rhon, 2.0e36, r'$\rho_n$', rotation=90, ha='center', va='center', size=8)
    txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))

    #ax.plot([rhond, rhond], [1.0e16, 1e40], "r", linestyle=lstyle)
    #txt = ax.text(rhond*0.5, 2.0e36, r'$\rho_{\mathrm{ND}}$', rotation=90, ha='center', va='center', size=8)


    if False:
        ax.set_ylim(1.0, 5.0)
        ax.set_yscale('linear')
        ax.set_ylabel(r'Adiabatic index $\gamma$')
        
        ax.fill_between([xmin, xmax], [1.0, 1.0], [10.0, 10.0], facecolor=col, color=None, alpha=1.0, edgecolor=col)
        
        #different ns structures
        y_text = 4.5
        ax.text(5.0e13, y_text, 'Inner\ncrust', rotation=0, ha='center', va='center', size=8)
        ax.text(6.0e14, y_text, 'Core', rotation=0, ha='center', va='center', size=8)
        
        lstyle = 'dotted'
        ax.plot([rhon, rhon], [1.0, 10.0], "r", linestyle=lstyle)
        txt = ax.text(rhon, 4.5, r'$\rho_n$', rotation=90, ha='center', va='center', size=8)
        txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))



    if True:
        i = 0
        for key, value in eosLib.iteritems():
            dense_eos = get_eos(key)
            eos = glue_crust_and_core( SLyCrust, dense_eos )

            rho = np.logspace(13, 18, 100)

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

            if True:
                #l, = ax.plot(rad, mass, color=col, linestyle=linestyle, alpha = 0.9)
                press = eos.pressures(rho)
                l, = ax.plot(rho, press, color=col, linestyle=linestyle, alpha=0.9)

            if False:

                gamma = np.zeros(len(rho))
                for j, r in enumerate(rho):
                    trope = eos._find_interval_given_density(r)
                    gamma[j] = trope.G

                #if key == 'SLy':
                #    col = 'darkorange'

                l, = ax.plot(rho, gamma, color=col, linestyle=linestyle, alpha=0.8)

                # labels for lines
                near_y = None
                near_x = 3.0e15
                rotation_offset=0.0
                offslabels = ['WFF2', 'APR3', 'SLy', 'MS1b']
                if key in offslabels:
                    near_x += 6.0e15
                label_line(l, key, near_y=near_y, near_x=near_x, rotation_offset=rotation_offset)

            i += 1


    #Exact SLy curve
    if False:
        from SLy import SLyGs
        SLyGs = SLyGs[SLyGs[:,0].argsort()] #sorting

        rhoSLy   = 10.0**SLyGs[:,0]
        gammaSLy = SLyGs[:,1]
        ax.plot( rhoSLy, gammaSLy, color='darkorange', alpha=1.0, linewidth=1.0, linestyle='dashed')

    if False:
        from SLy import SLyPs
        rhoSLy   = 10.0**SLyPs[:,0]
        pressSLy = 10.0**SLyPs[:,1]
        ax.plot( rhoSLy, pressSLy, color='darkorange', alpha=1.0, linewidth=1.5, linestyle='dashed')



    dense_eos = get_eos('SLy')
    eos = glue_crust_and_core( SLyCrust, dense_eos )
    rho = np.logspace(13, 18, 100)

    linestyle='solid'
    col = 'darkorange'
    if False:
        press = eos.pressures(rho)
        l, = ax.plot(rho, press, color=col, linestyle=linestyle, alpha=0.9)
    if False:
        gamma = np.zeros(len(rho))
        for j, r in enumerate(rho):
            trope = eos._find_interval_given_density(r)
            gamma[j] = trope.G
        l, = ax.plot(rho, gamma, color=col, linestyle=linestyle, alpha=0.8)






    plt.subplots_adjust(left=0.15, bottom=0.16, right=0.98, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig('core_eos.pdf')
    #plt.savefig('core_gamma.pdf')
