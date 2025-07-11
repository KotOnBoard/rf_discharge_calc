import math
import matplotlib.pyplot as plt
from scipy import optimize, integrate
import numpy as np
import json5
#from sympy import Symbol

ro_B = 2.34*(10**3)
M_B = 10.81
ro_Mo = 10.2*(10**3)
M_Mo = 95.95
ro_MoO3 = 4.9*(10**3)
M_MoO3 = 95.95+3.16

e = 1.6e-19
gas = "Ne"
p = 6 #Геометрия ?
l = 0.10 #electrode dist
R = 6.65e-2 #electrode diam
S1 = np.pi*R**2
S2 = 4*S1 #произвольная геометрическая ассиметрия 4
f = 10**8
Vrf = 22
Pwr = 300 #rf field params
V1 = Vrf/2
omega = 2*np.pi*f #потенциал плазмы и слоя в симм. случ
sm = 0.001
d = l-2*sm
k = 1.38e-23
e0 = 8.85e-12
N_A = 6.022e23
Ti = 300
Tn = 300
m = 9.1e-31
qe = 1.6e-19,  # Заряд электрона [Кл]
me = 9.11e-31,  # Масса электрона [кг]
Mi = 6.6335209e-26 - 9.1093837e-31,  # Масса иона Ar [кг]
eps_0 = 8.85e-12,  # Диэлектрическая постоянная [Ф/м]
k_B = 1.380649e-23  # Постоянная Больцмана [Дж/К]







def LoadGasPar(prname=False):
    global gas_params
    if (prname == False):
        gas_params = {
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
    else: 
        with open(f"conf/{prname}.json5", 'r') as file:
            gas_params = json5.load(file)
    return gas_params


LoadGasPar()
ng = (2*p)/(3*k*Ti)
lam_i = 1/(ng*gas_params[gas]['sig_i']) 
lam_i_d = lam_i/d

hl = 0.86*(3+d/(2*lam_i))**(-1/2)
hR = 0.8*(4+R/lam_i)**(-1/2)
deff = 1/2*(R*l/(R*hl+l*hR))
    

def AltTe():
    
    
    def LoadConf(cfname):
        global cf
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
    
    sol = optimize.root_scalar(dfr, bracket=[1, 7], x0=3, x1=5, xtol=1e-3, method='secant')
    return sol


"""
1. Блок расчёта Температуры электронов и констант скорости взаимодействия.
Нахождение энергии электронов через баланс ионизации и ухода через дрейф
"""
        
def TeCalc(verbose=False):
    """
    

    Returns
    -------
    None.

    """
    LoadGasPar()
    
    def TeEquation(Te_val, verbose=True):
        '''Установка интегрального уравнения Kiz и Te 
        для решения SciPy'''
        ub = np.sqrt(e*Te_val/gas_params[gas]['M'])
        if verbose == True: print(f'ub = {ub}')
        Kiz = Kiz_(Te_val, verbose)
        if verbose == True: 
            print(f'Te_val = {Te_val}')
            print(f"error = {Kiz-ub*((ng*deff)**-1)}")
        return Kiz-ub*((ng*deff)**-1)

    def Te_(verbose=False):
        Te_val = optimize.root_scalar(TeEquation, bracket=[1, 7], x0=3, x1=5, xtol=1e-3, method='secant')
        Te = Te_val.root
        if verbose == True: 
            print(f'Te = {Te}')
            print(Te_val)
        ub = np.sqrt(e*Te/gas_params[gas]['M'])
        return Te, ub

    def Kiz_(Te_val, verbose=False):
        ub = np.sqrt(e*Te_val/gas_params[gas]['M'])
        Kiz_integ = lambda v: (
            (gas_params[gas]['aiz']*((m*v**2/(2*e))/gas_params[gas]['xi_iz']-1))/
            ((m*(v**2)/(2*e))/gas_params[gas]['xi_iz']+gas_params[gas]['biz'])**gas_params[gas]['ciz']
            *v*np.exp(-(m*(v**2)/(2*e*Te_val)))*4*np.pi*v**2
            )  
        Kiz = (1e-20 * np.power((m/(2*np.pi*e*Te_val)), (3/2))*
        integrate.quad(Kiz_integ, (np.sqrt(2*e*gas_params[gas]['xi_iz']/m)), 10**7)[0])
        if verbose == True: 
            v = 10**7
            print(f"ub+ = {ub}")
            print(f"Kiz = {Kiz}")
            print(f"""Kiz_int = {
                (gas_params[gas]['aiz']*((m*v**2/(2*e))/gas_params[gas]['xi_iz']-1))/
                ((m*(v**2)/(2*e))/gas_params[gas]['xi_iz']+gas_params[gas]['biz'])**gas_params[gas]['ciz']
                *v*np.exp(-(m*(v**2)/(2*e*Te_val)))*4*np.pi*v**2
                }""")
            print(f"Kiz_integ = {integrate.quad(Kiz_integ, (np.sqrt(2*e*gas_params[gas]['xi_iz']/m)), 10**7)}")
        return Kiz

    def Kel_(Te_val, verbose=False):
        Kel_integ = lambda v: (
        (gas_params[gas]['a1']+gas_params[gas]['b1']*(m*v**2/(2*e))**gas_params[gas]['c1'])/
        (1+gas_params[gas]['d1']*(m*v**2/(2*e))**gas_params[gas]['e1'])+
        (gas_params[gas]['a2']+gas_params[gas]['b2']*(m*v**2/(2*e))**gas_params[gas]['c2'])/
        (1+gas_params[gas]['d2']*(m*v**2/(2*e))**gas_params[gas]['e2'])
        *v*np.power(2.79, -m*(v**2)/(2*e*Te_val))*4*np.pi*v**2
        )
        Kel = (1e-20*np.power((m/(2*np.pi*e*Te_val)), (3/2))*
        integrate.quad(Kel_integ, 0, 10**7)[0])
        if verbose == True: print(f"Kel = {Kel}")
        return Kel  

    def Kex_(Te_val, verbose=False):
        Kex_integ = lambda v: (
            (gas_params[gas]['aex']*((m*v**2/(2*e))/gas_params[gas]['xi_ex']-1))/
            ((m*(v**2)/(2*e))/gas_params[gas]['xi_ex']+gas_params[gas]['bex'])**gas_params[gas]['cex']
            *v*np.power(2.79, -m*(v**2)/(2*e*Te_val))*4*np.pi*v**2
            )    
        Kex = (1e-20*np.power((m/(2*np.pi*e*Te_val)), (3/2))*
        integrate.quad(Kex_integ, (np.sqrt(2*e*gas_params[gas]['xi_ex']/m)), 10**7)[0])
        if verbose == True: print(f"Kex = {Kex}")
        return Kex
    
    def TeEval():
        """
        Выводит графики для ручного анализа.
        """
        test1 = [np.sqrt(e*Te_val/gas_params[gas]['M'])*((ng*deff)**-1) for Te_val in np.arange(0.01, 4, 10**-2)]
        test2 = [Kiz_(Te_val) for Te_val in np.arange(0.01, 4, 10**-2)]
        plt.figure()
        plt.plot(np.arange(0.01, 4, 10**-2), test1, label='ub')
        plt.plot(np.arange(0.01, 4, 10**-2), test2, label='Kiz')
        plt.ylim(0,1e-15)
        plt.legend()
        plt.show()
    
    def plot_Te():
        en_range = np.arange(1, 7.1, 0.1)

        fig = plt.figure(figsize=(12, 5))
        ax = fig.add_subplot(1, 1, 1)
        minor_ticks = np.arange(0, 101, 4)
        ax.set_xticks(minor_ticks, minor=True)
        plt.plot(en_range, Kiz_(en_range), label='Kiz side')
        plt.plot(en_range, np.sqrt(e*en_range/gas_params[gas]['M'])*((ng*deff)**-1), label='ub side')
        plt.axvline(cf["Te"], color='cyan', linestyle=':')
        plt.grid(which='minor', linestyle='--')
        plt.grid(which='major', linestyle='--')
        plt.minorticks_on()
        plt.yscale('log')
        plt.xlabel('Temperature [eV]')
        plt.ylabel('Particle balance [m^3/s]')
        plt.legend()
        plt.show()
        
    Te, ub = Te_(verbose)
    
    Kiz = Kiz_(Te)
    Kel = Kel_(Te)
    Kex = Kex_(Te)
    TeEval()
    return Te, ub, Kiz, Kel, Kex
    
Te, ub, Kiz, Kel, Kex = TeCalc()
      
"""
Конец Блока 1.
"""

"""
Блок 2: Согласование мощности плазмы с целевой.
"""

def SCalc():
    
    def SabsEq(Vrf):
        global ub
        vm = Kel*ng
        xi_c = 1/Kiz*(Kiz*gas_params[gas]['xi_iz']+
                               Kex*gas_params[gas]['xi_ex']+
                               (Kel*3*m*Te/gas_params[gas]['M']))
        
        Sohm = (1.73*m*hl/(2*e)*e0*omega**2*vm*
                np.sqrt(1*Te)*np.sqrt(1*Vrf/2)*d)
        Sstoc = (0.45*np.sqrt(m/e)*e0*omega**2*
                 np.sqrt(1*Te)*Vrf/2)
        #ub = np.sqrt(e*Te/gas_params[gas]['M'])
        ns = ((Sohm + 2*Sstoc)/(2*e*ub*(xi_c + 2*Te + Te*
                                            np.sqrt(np.log(gas_params[gas]['M']/
                                                           (2*np.pi*m)))
                                + 0.5*Te)))
        V = 0.83*Vrf/2
        Ji = e*ns*ub
        sm = np.sqrt(0.82*e0*np.power(V, 3/2)/Ji*
                     np.sqrt(2*e/gas_params[gas]['M']))
        J1 = 1.23*omega*e0/sm*Vrf/2
        Sabs = 2*e*ns*ub*(V+xi_c + 2*Te + Te*
                                            np.sqrt(np.log(gas_params[gas]['M']/
                                                           (2*np.pi*m)))
                                + 0.5*Te)
        dPio = 2*e*ns*ub*(gas_params[gas]['xi_iz']+V)/Sabs
        return Sabs - Pwr, Ji, sm, dPio

    def Sabs(verbose=False):
        Sabs = optimize.root_scalar(SabsEq, x0=3, x1=5, xtol=1e-3, method='secant')
        if verbose == True: print(f'Sabs = {SabsEq.root}')
        return Sabs
    
    Sabs = Sabs()
    return Sabs




    
def Kiz_alt(Te_val, verbose=False):
    res = (0.01*0.07)*ng*(2.34**-14)*Te_val*np.exp(-17.44 / Te_val)
    if verbose == True: print(f'Kiz_alt = {res}')
    return res
def sym_par(Te_val, verbose=False):
    res = (0.01+0.03)*np.sqrt(e * Te_val / (6.6335209e-26 - 9.1093837e-31))
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
    











TeCalc(verbose=True)
print(f"Scalc = {SCalc()}")

"""
Поиск электрических параметров цепи.
"""

Vf = -(Te/2)*np.log(gas_params[gas]['M']/(m*2*np.pi/(1.247**2)))        
beta = 1
sm2 = (S1+S2)/(2*S2)*np.power(S2/S1, beta/2-0.25)*sm
sm1 = (S2/S1)**(1.5-beta)
ns2 = np.power(S1/S2, beta/2-0.25)
ns1 = np.power(S2/S1, beta-0.5)*ns2
Ji1 = e*ns1*ub
Ji2 = e*ns2*ub
C1 = e0*S1/sm1
C2 = e0*S2/sm2
Vp_amp = C1/(C1+C2)*Vrf
Vp_avg = Vp_amp+np.abs(Vf)-Te/2*np.log(2*np.pi*Vp_amp/Te)
Vp_symm = Vrf/2+np.abs(Vf)-Te/2*np.log(np.pi*Vrf/Te)
Ubias = 2*Vp_avg-Vrf
V1_avg = Vp_avg-Ubias
V2_avg = Vp_avg

"""
Расчёт dE & N(E), s(w), Vdc(w)
"""
dEi1 = 4*np.power(e*V1_avg, 1.5)/(omega*sm1*np.sqrt(2*gas_params[gas]['M']))
dEi2 = 4*np.power(e*V2_avg, 1.5)/(omega*sm2*np.sqrt(2*gas_params[gas]['M']))
dEi_symm = 4*np.power(e*Vp_symm, 1.5)/(omega*sm2*np.sqrt(2*gas_params[gas]['M']))
Ns1 = Ji1/e
#Ns_symm = Ji/e
Ns2 = Ji2/e
