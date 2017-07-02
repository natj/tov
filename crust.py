from polytropes import monotrope
from polytropes import polytrope
import units as cgs


##################################################
#SLy (Skyrme) crust
KSLy = [6.80110e-9, 1.06186e-6, 5.32697e1, 3.99874e-8] #Scaling constants
GSLy = [1.58425, 1.28733, 0.62223, 1.35692] #polytropic indices
RSLy = [1.e4, 2.44034e7, 3.78358e11, 2.62780e12 ] #transition depths

tropes = []
trans = []

pm = None
for (K, G, r) in zip(KSLy, GSLy, RSLy):
    m = monotrope(K*cgs.c**2, G)
    tropes.append( m )

    #correct transition depths to avoid jumps
    if not(pm == None):
        rho_tr = (m.K / pm.K )**( 1.0/( pm.G - m.G ) )
        #print rho_tr, np.log10(rho_tr), r, rho_tr/r
    else:
        rho_tr = r
    pm = m

    trans.append(rho_tr)

#Create crust using polytrope class
SLyCrust = polytrope(tropes, trans)
