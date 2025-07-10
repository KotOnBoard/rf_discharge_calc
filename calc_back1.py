import matplotlib.pyplot as plt
from scipy import optimize, integrate
import numpy as np
import json5
from sympy import Symbol

def AltTe():
    
    
    def LoadConf(cfname):
        with open(f"conf/{cfname}.json5", 'r') as file:
            cf = json5.load(file)
        return cf
    
    def Kiz_novect(a_Te):
        res = 2.34e-14 * a_Te ** 0.59 * np.exp(-17.44 / a_Te)
        return res
    
    def u_Bohm_novect(a_Te):
        return np.sqrt(ct["qe"] * a_Te / ct["Mi"])
    
    def dfr(a_Te):
        return cf["Vp"] * cf["ng"] * Kiz_novect(a_Te) - (cf["Ae"] + cf["Ag"]) * u_Bohm_novect(a_Te)
    
    cfname = "as1"
    cf = LoadConf(cfname)
    ct = {
        "qe": 1.6e-19,  # Заряд электрона [Кл]
        "me": 9.11e-31,  # Масса электрона [кг]
        "Mi": 6.6335209e-26 - 9.1093837e-31,  # Масса иона Ar [кг]
        "eps_0": 8.85e-12,  # Диэлектрическая постоянная [Ф/м]
        "k_B": 1.380649e-23  # Постоянная Больцмана [Дж/К]
    }
    
    cf["ng"] = cf["p0"] / (ct["k_B"] * cf["T0"])  # Концентрация буферного газа [м^-3]
    cf["Vp"] = cf["Ae"] * cf["l_B"]  # Объем плазменного столба
    cf["fE"] = cf["Ae"] / (cf["Ae"] + cf["Ag"])  # Весовой коэф. слоя E
    cf["fG"] = cf["Ag"] / (cf["Ae"] + cf["Ag"])  # Весовой коэф. слоя G
    cf["Tf"] = 1 / cf["f0"]  # Период ВЧ поля [с]
    
    sol = optimize.root_scalar(dfr, bracket=[1, 7], x0=3, x1=5, xtol=1e-5, method='secant')
    return sol

class par():
    """
    Контейнер для сбора и хранения вводных параметров и таблиц.
    """
    def __init__(self):
        self.ro_B = 2.34*(10**3)
        self.M_B = 10.81
        self.ro_Mo = 10.2*(10**3)
        self.M_Mo = 95.95
        self.ro_MoO3 = 4.9*(10**3)
        self.M_MoO3 = 95.95+3.16

        self.e = 1.6e-19
        self.gas = "Ne"
        self.p = 6 #Геометрия ?
        self.l = 0.10 #electrode dist
        self.R = 6.65e-2 #electrode diam
        self.S1 = np.pi*self.R**2
        self.S2 = 4*self.S1 #произвольная геометрическая ассиметрия 4
        self.f = 10**8
        self.Vrf = 22
        self.Pwr = 300 #rf field params
        self.V1 = self.Vrf/2
        self.omega = 2*np.pi*self.f #потенциал плазмы и слоя в симм. случ
        self.sm = 0.001
        self.d = self.l-2*self.sm
        self.k = 1.38e-23
        self.e0 = 8.85e-12
        self.N_A = 6.022e+23
        self.Ti = 300
        self.Tn = 300
        self.m = 9.1e-31
        self.qe = 1.6e-19,  # Заряд электрона [Кл]
        self.me = 9.11e-31,  # Масса электрона [кг]
        self.Mi = 6.6335209e-26 - 9.1093837e-31,  # Масса иона Ar [кг]
        self.eps_0 = 8.85e-12,  # Диэлектрическая постоянная [Ф/м]
        self.k_B = 1.380649e-23  # Постоянная Больцмана [Дж/К]
        self.gas_params = {
               "Ne": {
                   "xi_ex": 18.7,
                   "xi_iz": 21.55,
                   "M": 32e-27,
                   "xi_ex": 18.7,
                   "sig_i": 3.93e-20*(273*300**-1),
                   "a1": 0,
                   "b1": 38.7,
                   "c1": 1,
                   "d1": 276,
                   "e1": 1.64,
                   "a2": 0.31,
                   "b2": 2.99,
                   "c2": 0.5,
                   "d2": 0.2,
                   "e2": 1.93,
                   "aex": 1.5,
                   "bex": 1.98,
                   "cex": 1.85,
                   "aiz": 20.11,
                   "biz": 6.34,
                   "ciz": 2
                   },
               "He": {
                   "xi_ex": 19.77,
                   "xi_iz": 24.58,
                   "M": 6.65e-27,
                   "sig_i": 6.5e-20*(273*300**-1),
                   "a1": 0,
                   "b1": 7.19,
                   "c1": 1,
                   "d1": 3.76,
                   "e1": 2.79,
                   "a2": 5.16,
                   "b2": 6.09,
                   "c2": 0.41,
                   "d2": 15,
                   "e2": 1.91,
                   "aex": 0.99,
                   "bex": 0.63,
                   "cex": 1.75,
                   "aiz": 3.95,
                   "biz": 2.48,
                   "ciz": 1.91
                   },
               "Ar": {
                   "xi_ex": 13.17,
                   "xi_iz": 15.76,
                   "M": 6.63e-26,
                   "sig_i": 3.6e-20*(273*300**-1),
                   "a1": 0,
                   "b1": 24.1,
                   "c1": 1,
                   "d1": 1.03,
                   "e1": 2.83,
                   "a2": 7.76,
                   "b2": -65.5,
                   "c2": 0.455,
                   "d2": 1961,
                   "e2": 1.37,
                   "aex": 6.48,
                   "bex": 1.83,
                   "cex": 1.81,
                   "aiz": 30.1,
                   "biz": 2.51,
                   "ciz": 1.86
                   }
               }
        
        

""" Нахождение энергии электронов 
через баланс ионизации и ухода через дрейф """

class val():
    """
    Контейнер для хранения расчётных параметров.
    """
    def __init__(self):
        self.ng = (2*par.p)/(3*par.k*par.Ti)
        self.lam_i = 1/(self.ng*par.gas_params[par.gas]['sig_i']) 
        self.lam_i_d = self.lam_i/par.d
        
        self.hl = 0.86*(3+par.d/(2*self.lam_i))**(-1/2)
        self.hR = 0.8*(4+par.R/self.lam_i)**(-1/2)
        self.deff = 1/2*(par.R*par.l/(par.R*self.hl+par.l*self.hR))


def TeEquation(Te_val, verbose=True):
    '''Установка интегрального уравнения Kiz и Te 
    для решения SciPy'''
    ub = np.sqrt(par.e*Te_val/par.gas_params[par.gas]['M'])
    if verbose == True: print(f'ub = {ub}')
    Kiz = Kiz_(Te_val, verbose)
    if verbose == True: 
        print(f'Te_val = {Te_val}')
        print(f"error = {Kiz-ub*((val.ng*val.deff)**-1)}")
    return Kiz-ub*((val.ng*val.deff)**-1)#, Kiz

def InitParams():
    """
    Инициализация значений начальных параметров par и рачётных значений val.
    """
    par.__init__(par)
    val.__init__(val)


def TeEval():
    """
    Выводит графики для ручного анализа.
    """
    test1 = [np.sqrt(par.e*Te_val/par.gas_params[par.gas]['M'])*((val.ng*val.deff)**-1) for Te_val in np.arange(0.01, 4, 10**-2)]
    test2 = [TeEquation(Te_val)[1] for Te_val in np.arange(0.01, 4, 10**-2)]
    plt.figure()
    plt.plot(np.arange(0.01, 4, 10**-2), test1, label='ub')
    plt.plot(np.arange(0.01, 4, 10**-2), test2, label='Kiz')
    plt.legend()
    plt.show()
    
def Kiz_alt(Te_val, verbose=False):
    res = (0.01*0.07)*val.ng*(2.34**-14)*Te_val*np.exp(-17.44 / Te_val)
    if verbose == True: print(f'Kiz_alt = {res}')
    return res
def sym_par(Te_val, verbose=False):
    res = (0.01+0.03)*np.sqrt(par.e * Te_val / (6.6335209e-26 - 9.1093837e-31))
    if verbose == True: print(f' sym_par = {res}')
    return res

def AltTeEval():
    """
    Выводит графики для ручного анализа.

    Returns
    -------
    None.

    """
    test3 = [Kiz_alt(Te_val) for Te_val in np.arange(0.01, 10, 10**-1)]
    test4 = [sym_par(Te_val) for Te_val in np.arange(0.01, 10, 10**-1)]
    plt.plot(np.arange(0.01, 10, 10**-1), test3, label='Kiz_Alt')
    plt.plot(np.arange(0.01, 10, 10**-1), test4, label='sym_par')
    plt.ylim(0,200)
    plt.legend()
    plt.show() 
    
def Te_(verbose=False):
    InitParams()
    Te = optimize.root_scalar(TeEquation, bracket=[1, 7], x0=3, x1=5, xtol=1e-30, method='secant')
    if verbose == True: print(f'res = {Te.root}')
    return Te

def Kiz_(Te_val, verbose=False):
    ub = np.sqrt(par.e*Te_val/par.gas_params[par.gas]['M'])
    Kiz_integ = lambda v: (
        (par.gas_params[par.gas]['aiz']*((par.m*v**2/(2*par.e))/par.gas_params[par.gas]['xi_iz']-1))/
        ((par.m*(v**2)/(2*par.e))/par.gas_params[par.gas]['xi_iz']+par.gas_params[par.gas]['biz'])**par.gas_params[par.gas]['ciz']
        *v*np.exp(-(par.m*(v**2)/(2*par.e*Te_val)))*4*np.pi*v**2
        )  
    Kiz = (1e-20 * np.power((par.m/(2*np.pi*par.e*Te_val)), (3/2))*
    integrate.quad(Kiz_integ, (np.sqrt(2*par.e*par.gas_params[par.gas]['xi_iz']/par.m)), 10**7)[0])
    if verbose == True: 
        v = 3*10**6
        print(f"ub+ = {ub}")
        print(f"Kiz = {Kiz}")
        print(f"""TESTING = {
            ub*(val.ng*val.deff)**-1
            } 5.43e+22""")
        print(f"Kiz_integ = {integrate.quad(Kiz_integ, (np.sqrt(2*par.e*par.gas_params[par.gas]['xi_iz']/par.m)), 10**7)}")
    return Kiz
    
def Kel_(Te_val, verbose=False):
    Kel_integ = lambda v: (
    (par.gas_params[par.gas]['a1']+par.gas_params[par.gas]['b1']*(par.m*v**2/(2*par.e))**par.gas_params[par.gas]['c1'])/
    (1+par.gas_params[par.gas]['d1']*(par.m*v**2/(2*par.e))**par.gas_params[par.gas]['e1'])+
    (par.gas_params[par.gas]['a2']+par.gas_params[par.gas]['b2']*(par.m*v**2/(2*par.e))**par.gas_params[par.gas]['c2'])/
    (1+par.gas_params[par.gas]['d2']*(par.m*v**2/(2*par.e))**par.gas_params[par.gas]['e2'])
    *v*np.power(2.79, -par.m*(v**2)/(2*par.e*Te_val))*4*np.pi*v**2
    )
    Kel = (10e-20*np.power((par.m/(2*np.pi*par.e*Te_val)), (3/2))*
    integrate.quad(Kel_integ, (0, 10**7)[0]))
    if verbose == True: print(f"Kel = {Kel}")
    return Kel  

def Kex_(Te_val, verbose=False):
    Kex_integ = lambda v: (
        (par.gas_params[par.gas]['aex']*((par.m*v**2/(2*par.e))/par.gas_params[par.gas]['xi_ex']-1))/
        ((par.m*(v**2)/(2*par.e))/par.gas_params[par.gas]['xi_ex']+par.gas_params[par.gas]['bex'])**par.gas_params[par.gas]['cex']
        *v*np.power(2.79, -par.m*(v**2)/(2*par.e*Te_val))*4*np.pi*v**2
        )    
    Kex = (10e-20*np.power((par.m/(2*np.pi*par.e*Te_val)), (3/2))*
    integrate.quad(Kex_integ, (np.sqrt(2*par.e*par.gas_params[par.gas]['xi_ex']/par.m)), 10**7)[0])
    if verbose == True: print(f"Kex = {Kex}")
    return Kex

InitParams()
print(Te_(verbose=True))
Te = Te_()
print(Kiz_(3.26, verbose=True))
print(val.deff)
"""
Te_val = 1
vm = Kel_(Te_val)
xi_c = 1/Kiz_(Te_val)*(Kiz_(Te_val)*par.gas_params[par.gas]['xi_iz']+
                       Kex_(Te_val)*par.gas_params[par.gas]['xi_ex']+
                       (Kel_(Te_val)*3*par.m*Te_val/par.M))


def dSNormaliser(Te_val):
    while (np.abs(dS) < 0.5):
        Sohm = (1.73*par.m*val.hl/(2*par.e)*par.e0*par.omega**2*vm*
                np.sqrt(1*Te_val)*np.sqrt(1*par.V1)*par.d)
        Sstoc = (0.45*np.sqrt(par.m/par.e)*par.e0*par.omega**2*
                 np.sqrt(1*Te_val)*par.V1)
        ub = np.sqrt(par.e*Te_val/par.gas_params[par.gas]['M'])
        ns = ((Sohm + 2*Sstoc)/(2*par.e*ub*(xi_c + 2*Te_val + Te_val*
                                            np.sqrt(np.log(par.gas_params[par.gas]['M']/
                                                           (2*np.pi*par.m)))
                                + 0.5*Te_val)))
        V = 0.83*par.V1
        Ji = par.e*ns*ub
        sm = np.sqrt(0.82*par.e0*np.power(V, 3/2)/Ji*
                     np.sqrt(2*par.e/par.gas_params[par.gas]['M']))
        J1 = 1.23*par.omega*par.e0/sm*par.V1
        Sabs = 2*par.e*ns*ub*(V+xi_c + 2*Te_val + Te_val*
                                            np.sqrt(np.log(par.gas_params[par.gas]['M']/
                                                           (2*np.pi*par.m)))
                                + 0.5*Te_val)
        dS = Sabs - par.Pwr
        if (dS>10): 
            par.Vrf-=2
        elif (dS<-5):
            par.Vrf+=3
        elif (dS>0.4): 
            par.Vrf-=0.03
        elif (dS<-0.4):
            par.Vrf+=0.03
        par.V1=par.Vrf/2
        dPio = 2*par.e*ns*ub*(par.gas_params[par.gas]['xi_iz']+par.V)/Sabs

Vf = -(Te_val/2)*np.log(par.gas_params[par.gas]['M']/(par.m*2*np.pi/(1.247**2)))        
beta = 1
sm = Symbol('sm')
sm2 = (par.S1+par.S2)/(2*par.S2)*np.power(par.S2/par.S1, beta/2-0.25)*sm
"""