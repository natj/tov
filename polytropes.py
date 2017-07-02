import units as cgs

#Monotropic eos
class monotrope:
    
    #transition continuity constant 
    a = 0.0
    cgsunits = cgs.c**2

    def __init__(self, K, G):
        self.K = K / self.cgsunits
        self.G = G
        self.n = 1.0/(G - 1)

    #pressure P(rho)
    def pressure(self, rho):
        return self.cgsunits * self.K * rho**self.G

    #energy density mu(rho)
    def edens(self, rho):
        return (1.0 + self.a)*rho + (self.K/(self.G - 1)) * rho**self.G

    #for inverse functions lets define rho(P)
    def rho(self, press):
        if press < 0.0:
            return 0.0
        return ( press/self.cgsunits/self.K )**(1 / self.G)


# Piecewise polytropes
class polytrope:
    
    def __init__(self, tropes, trans, prev_trope = None ):
        self.tropes      = tropes
        self.transitions = trans

        self.prs  = []
        self.eds  = []
        
        for (trope, transition) in zip(self.tropes, self.transitions):

            if not( prev_trope == None ):
                trope.a = self._ai( prev_trope, trope, transition )
            else:
                transition = 0.0


            ed = trope.edens(transition) 
            pr = trope.pressure(transition)

            self.prs.append( pr )
            self.eds.append( ed )

            prev_ed = ed
            prev_tr = transition
            prev_trope = trope


    def _ai(self, pm, m, tr):
        return pm.a + (pm.K/(pm.G - 1))*tr**(pm.G-1) - (m.K/(m.G - 1))*tr**(m.G-1)

    def _find_interval_given_density(self, rho):
        if rho <= self.transitions[0]:
            return self.tropes[0]

        for q in range( len(self.transitions) - 1 ):
            if self.transitions[q] <= rho < self.transitions[q+1]:
                return self.tropes[q]

        return self.tropes[-1]

    #inverted equations as a function of pressure
    def _find_interval_given_pressure(self, press):
        if press <= self.prs[0]:
            return self.tropes[0]

        for q in range( len(self.prs) - 1):
            if self.prs[q] <= press < self.prs[q+1]:
                return self.tropes[q]

        return self.tropes[-1]


    ################################################## 
    def pressure(self, rho):
        trope = self._find_interval_given_density(rho)
        return trope.pressure(rho)

    #vectorized version
    def pressures(self, rhos):
        press = []
        for rho in rhos:
            pr = self.pressure(rho)
            press.append( pr )
        return press

    def edens_inv(self, press):
        trope = self._find_interval_given_pressure(press)
        rho = trope.rho(press)
        return trope.edens(rho)

    def rho(self, press):
        trope = self._find_interval_given_pressure(press)
        return trope.rho(press)


