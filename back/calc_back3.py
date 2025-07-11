import time
import matplotlib.pyplot as plt
from scipy import optimize, integrate
import numpy as np
import json5
#from sympy import Symbol
start_time = time.time()


def ParLoader(prname=False, gcname=False):
    """
    Function for automatic parameters filling.

    Parameters
    ----------
    prname : str, optional
        Name of the target json to obtain values of initial parameters from.
        Target file must be within /conf.
        The default is filling dict with preexisting arbitrary parameters.
    gcname : str, optional
        Name of the target json to obtain values of gas parameters from.
        Target file must be within /conf.
        The default is filling dict with preexisting arbitrary parameters.

    Returns
    -------
    par : dict
        Initial parameters and other parameters calculated from initial.
    gas_params : dict
        Parameters for used gas and it's collision approximations.
    out : dict
        Empty dictionary for all calculated parameters for analysis.

    """           
    def LoadGasPar(gcname=False):
        """
        Function for fetching approximation parameters sorted by used gas

        Parameters
        ----------
        gcname : str, optional
            Name of the target json to obtain values of gas parameters from.
            Target file must be within /conf.
            The default is filling dict with preexisting arbitrary parameters.

        Returns
        -------
        gas_params : dict
            Dictionary containing all approximation & gas parameters.

        """
        global gas_params
        if (gcname == False):
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
            with open(f"conf/{gcname}.json5", 'r') as file:
                gas_params = json5.load(file)
        return gas_params
    
    def LoadInitPar(prname=False):
        """
        Function for fetching values of initial physical parameters.

        Parameters
        ----------
        prname : str, optional
            Name of the target json to obtain values of initial parameters from.
            Target file must be within /conf.
            The default is filling dict with preexisting arbitrary parameters.

        Returns
        -------
        par : dict
            Dictionary containing all initial parameters.

        """
        if (prname==False):
            par = {
                'ro_B': 2.34*(10**3),
                'M_B': 10.81,
                'ro_Mo': 10.2*(10**3),
                'M_Mo': 95.95,
                'ro_MoO3': 4.9*(10**3),
                'M_MoO3': 95.95+3.16,
                
                'e': 1.6e-19,
                'gas': "Ne",
                'p': 6, #Геометрия ?
                'l': 0.10, #electrode dist
                'R': 6.65e-2, #electrode diam
                
                'f': 10**8,
                'Vrf': 22,
                'Pwr': 300, #rf field params
                
                 #потенциал плазмы и слоя в симм. случ
                'sm': 0.001,
                
                'k': 1.38e-23,
                'e0': 8.85e-12,
                'N_A': 6.022e23,
                'Ti': 300,
                'Tn': 300,
                'm': 9.1e-31,
                'qe': 1.6e-19,  # Заряд электрона [Кл]
                'me': 9.11e-31,  # Масса электрона [кг]
                'Mi': 6.6335209e-26 - 9.1093837e-31,  # Масса иона Ar [кг]
                'eps_0': 8.85e-12,  # Диэлектрическая постоянная [Ф/м]
                'k_B': 1.380649e-23  # Постоянная Больцмана [Дж/К]
                }
        else: 
            with open(f"conf/{prname}.json5", 'r') as file:
                par = json5.load(file)
        par['S1'] = np.pi*par['R']**2
        par['S2'] = 4*par['S1'] #произвольная геометрическая ассиметрия 4
        par['V1'] = par['Vrf']/2
        par['omega'] = 2*np.pi*par['f']
        par['d'] = par['l']-2*par['sm']
        par['ng'] = (2*par['p'])/(3*par['k']*par['Ti'])
        par['lam_i'] = 1/(par['ng']*gas_params[par['gas']]['sig_i']) 
        par['lam_i_d'] = par['lam_i']/par['d']
        par['hl'] = 0.86*(3+par['d']/(2*par['lam_i']))**(-1/2)
        par['hR'] = 0.8*(4+par['R']/par['lam_i'])**(-1/2)
        par['deff'] = 1/2*(par['R']*par['l']/(par['R']*par['hl']+par['l']*par['hR']))
                    
        return par
    
    global gas_params
    gas_params = LoadGasPar()
    global par
    par = LoadInitPar()
    global out
    out = {}
    
    return par, gas_params, out


def AltTe():
    """
    Backup alternative function for Te evaluation produced by D.S.Samsonov.
    Capable of calculating only in Argon atmosphere.
    Input data is very different and must be checked and modified manually.

    Returns
    -------
    sol : float
        Alternative Electron temperature.

    """
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
    
    def Kiz_alt(Te_val, verbose=False):
        res = (0.01*0.07)*par['ng']*(2.34**-14)*Te_val*np.exp(-17.44 / Te_val)
        if verbose == True: print(f'Kiz_alt = {res}')
        return res
    def sym_par(Te_val, verbose=False):
        res = (0.01+0.03)*np.sqrt(par['e'] * Te_val / (6.6335209e-26 - 9.1093837e-31))
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
    Function for calculating Electron temperature, rate constants and some additional parameters.

    Parameters
    ----------
    verbose : bool, optional
        Prints values of most variables for every iteration.
        The default is False.

    Returns
    -------
    Te : float
        Electron temperature.
    ub : float
        ???
    Kiz : float
        Ionisation rate constant.
    Kel : float
        Elactic collision rate constant.
    Kex : float
        Exitation rate constant.
        DESCRIPTION.

    """    
    def TeEquation(Te_val, verbose=False):
        """
        Inegral equation to be solved for obtaining Electron temperature (Te).

        Parameters
        ----------
        Te_val : float
            Arbittrary Te value.
        verbose : bool, optional
            Prints values of most variables for every iteration.
            The default is False.

        Returns
        -------
        Kiz-ub*((par['ng']*par['deff'])**-1) : float
            Value of equation that must be (near) zero with correct Te.

        """
        ub = np.sqrt(par['e']*Te_val/gas_params[par['gas']]['M'])
        if verbose == True: print(f'ub = {ub}')
        Kiz = Kiz_(Te_val, verbose=verbose)
        if verbose == True: 
            print(f'Te_val = {Te_val}')
            print(f"error = {Kiz-ub*((par['ng']*par['deff'])**-1)}")
        return Kiz-ub*((par['ng']*par['deff'])**-1)

    def Te_(verbose=False):
        """
        Function for solving equation for Electron temperature (Te).

        Parameters
        ----------
        verbose : bool, optional
            Prints values of most variables for every iteration.
            The default is False.

        Returns
        -------
        Te : float
            Calculated Electron temperature.
        ub : float
            ???

        """
        fun = lambda x: TeEquation(x, verbose=verbose)
        Te_val = optimize.root_scalar(fun, bracket=[1, 7], x0=3, x1=5, xtol=1e-3, method='secant')
        Te = Te_val.root
        if verbose == True: 
            print(f'Te = {Te}')
            print(Te_val)
        ub = np.sqrt(par['e']*Te/gas_params[par['gas']]['M'])
        return Te, ub

    def Kiz_(Te_val, verbose=False):
        """
        Calculation of Ionisation rate constant (Kiz).

        Parameters
        ----------
        Te_val : float
            Arbittrary Te value.
        verbose : bool, optional
            Prints values of most variables for every iteration.
            The default is False.

        Returns
        -------
        Kiz : float
            Ionisation rate constant.

        """
        ub = np.sqrt(par['e']*Te_val/gas_params[par['gas']]['M'])
        Kiz_integ = lambda v: (
            (gas_params[par['gas']]['aiz']*((par['m']*v**2/(2*par['e']))/gas_params[par['gas']]['xi_iz']-1))/
            ((par['m']*(v**2)/(2*par['e']))/gas_params[par['gas']]['xi_iz']+gas_params[par['gas']]['biz'])**gas_params[par['gas']]['ciz']
            *v*np.exp(-(par['m']*(v**2)/(2*par['e']*Te_val)))*4*np.pi*v**2
            )  
        Kiz = (1e-20 * np.power((par['m']/(2*np.pi*par['e']*Te_val)), (3/2))*
        integrate.quad(Kiz_integ, (np.sqrt(2*par['e']*gas_params[par['gas']]['xi_iz']/par['m'])), 10**7)[0])
        if verbose == True: 
            v = 10**7
            print(f"ub+ = {ub}")
            print(f"Kiz = {Kiz}")
            print(f"""Kiz_int = {
                (gas_params[par['gas']]['aiz']*((par['m']*v**2/(2*par['e']))/gas_params[par['gas']]['xi_iz']-1))/
                ((par['m']*(v**2)/(2*par['e']))/gas_params[par['gas']]['xi_iz']+gas_params[par['gas']]['biz'])**gas_params[par['gas']]['ciz']
                *v*np.exp(-(par['m']*(v**2)/(2*par['e']*Te_val)))*4*np.pi*v**2
                }""")
            print(f"Kiz_integ = {integrate.quad(Kiz_integ, (np.sqrt(2*par['e']*gas_params[par['gas']]['xi_iz']/par['m'])), 10**7)}")
        return Kiz

    def Kel_(Te_val, verbose=False):
        """
        Calculation of Elastic collision rate constant (Kiz).

        Parameters
        ----------
        Te_val : float
            Arbittrary Te value.
        verbose : bool, optional
            Prints values of most variables for every iteration.
            The default is False.

        Returns
        -------
        Kel : float
            Elastic collision rate constant.

        """
        Kel_integ = lambda v: (
        (gas_params[par['gas']]['a1']+gas_params[par['gas']]['b1']*(par['m']*v**2/(2*par['e']))**gas_params[par['gas']]['c1'])/
        (1+gas_params[par['gas']]['d1']*(par['m']*v**2/(2*par['e']))**gas_params[par['gas']]['e1'])+
        (gas_params[par['gas']]['a2']+gas_params[par['gas']]['b2']*(par['m']*v**2/(2*par['e']))**gas_params[par['gas']]['c2'])/
        (1+gas_params[par['gas']]['d2']*(par['m']*v**2/(2*par['e']))**gas_params[par['gas']]['e2'])
        *v*np.power(2.79, -par['m']*(v**2)/(2*par['e']*Te_val))*4*np.pi*v**2
        )
        Kel = (1e-20*np.power((par['m']/(2*np.pi*par['e']*Te_val)), (3/2))*
        integrate.quad(Kel_integ, 0, 10**7)[0])
        if verbose == True: print(f"Kel = {Kel}")
        return Kel  

    def Kex_(Te_val, verbose=False):
        """
        Calculation of Exitation rate constant (Kiz).

        Parameters
        ----------
        Te_val : float
            Arbittrary Te value.
        verbose : bool, optional
            Prints values of most variables for every iteration.
            The default is False.

        Returns
        -------
        Kex : float
            Exitation rate constant.

        """
        Kex_integ = lambda v: (
            (gas_params[par['gas']]['aex']*((par['m']*v**2/(2*par['e']))/gas_params[par['gas']]['xi_ex']-1))/
            ((par['m']*(v**2)/(2*par['e']))/gas_params[par['gas']]['xi_ex']+gas_params[par['gas']]['bex'])**gas_params[par['gas']]['cex']
            *v*np.power(2.79, -par['m']*(v**2)/(2*par['e']*Te_val))*4*np.pi*v**2
            )    
        Kex = (1e-20*np.power((par['m']/(2*np.pi*par['e']*Te_val)), (3/2))*
        integrate.quad(Kex_integ, (np.sqrt(2*par['e']*gas_params[par['gas']]['xi_ex']/par['m'])), 10**7)[0])
        if verbose == True: print(f"Kex = {Kex}")
        return Kex
    
    def TeEval(a=1e-2, b=10, c=1e-2, d=2e-16):
        """
        Function for manual approximation of Te.
        Plots Kiz & ub side of equation with value of Te on the intersection.
        Intended for testing purposes.

        Parameters
        ----------
        a : float, optional
            Lower range of Te range. The default is 1e-2.
        b : float, optional
            Upper range of Te range. The default is 10.
        c : float, optional
            Step of Te value, resolution of curves. The default is 1e-2.
        d : float, optional
            Upper limit of y axis. The default is 2e-16.

        Returns
        -------
        None.

        """
        test1 = [np.sqrt(par['e']*Te_val/gas_params[par['gas']]['M'])*((par['ng']*par['deff'])**-1) 
                 for Te_val in np.arange(a, b, c)]
        test2 = [Kiz_(Te_val) for Te_val in np.arange(a, b, c)]
        plt.figure()
        plt.plot(np.arange(a, b, c), test1, label='ub side')
        plt.plot(np.arange(a, b, c), test2, label='Kiz side')
        plt.ylim(0,d)
        plt.legend()
        plt.show()
    
    def plot_Te():
        """
        Alternative plotting strategy produced by D.S.Samsonov.
        Intended for check-up purposes.

        Returns
        -------
        None.

        """
        en_range = np.arange(1, 7.1, 0.1)

        fig = plt.figure(figsize=(12, 5))
        ax = fig.add_subplot(1, 1, 1)
        minor_ticks = np.arange(0, 101, 4)
        ax.set_xticks(minor_ticks, minor=True)
        plt.plot(en_range, Kiz_(en_range), label='Kiz side')
        plt.plot(en_range, np.sqrt(par['e']*en_range/gas_params[par['gas']]['M'])*((par['ng']*par['deff'])**-1), label='ub side')
        plt.axvline(cf["Te"], color='cyan', linestyle=':')
        plt.grid(which='minor', linestyle='--')
        plt.grid(which='major', linestyle='--')
        plt.minorticks_on()
        plt.yscale('log')
        plt.xlabel('Temperature [eV]')
        plt.ylabel('Particle balance [m^3/s]')
        plt.legend()
        plt.show()
        
    Te, ub = Te_(verbose=verbose)
    Kiz = Kiz_(Te)
    Kel = Kel_(Te)
    Kex = Kex_(Te)
    if (verbose==True): TeEval()
    return Te, ub, Kiz, Kel, Kex
    

     

"""
Блок 2: Согласование мощности плазмы с целевой.
"""

def VrfCalc(verbose=False):
    """
    Function for calculating absorbed power and matching it by varying RF voltage. 

    Returns
    -------
    Vrf : float
        RF voltage with matching RF power & Power absobed by plasma.

    """
    def VrfEq(Vrf, verbose=False):
        out['vm'] = out['Kel']*par['ng']
        out['xi_c'] = 1/out['Kiz']*(out['Kiz']*gas_params[par['gas']]['xi_iz']+
                               out['Kex']*gas_params[par['gas']]['xi_ex']+
                               (out['Kel']*3*par['m']*out['Te']/gas_params[par['gas']]['M']))
        
        out['Sohm'] = (1.73*par['m']*par['hl']/(2*par['e'])*par['e0']*par['omega']**2*out['vm']*
                np.sqrt(1*out['Te'])*np.sqrt(1*par['Vrf']/2)*par['d'])
        out['Sstoc'] = (0.45*np.sqrt(par['m']/par['e'])*par['e0']*par['omega']**2*
                 np.sqrt(1*out['Te'])*par['Vrf']/2)
        out['ns'] = ((out['Sohm']+2*out['Sstoc'])/
                     (2*par['e']*out['ub']*
                      (out['xi_c']+2*out['Te']+out['Te']*
                       np.sqrt(np.log(gas_params[par['gas']]['M']/
                                      (2*np.pi*par['m'])))+0.5*out['Te'])))
        out['V'] = 0.83*par['Vrf']/2
        out['Ji'] = par['e']*out['ns']*out['ub']
        out['sm'] = np.sqrt(0.82*par['e0']*np.power(out['V'], 3/2)/out['Ji']*
                     np.sqrt(2*par['e']/gas_params[par['gas']]['M']))
        out['J1'] = 1.23*par['omega']*par['e0']/out['sm']*par['Vrf']/2
        out['Sabs'] = (2*par['e']*out['ns']*out['ub']*
        (out['V']+out['xi_c']+2*out['Te']+out['Te']*
         np.sqrt(np.log(gas_params[par['gas']]['M']/
                        (2*np.pi*par['m'])))+0.5*out['Te']))
        out['dPio'] = 2*par['e']*out['ns']*out['ub']*(gas_params[par['gas']]['xi_iz']+out['V'])/out['Sabs']
        return out['Sabs'] - par['Pwr']

    def Vrf_(verbose=False):
        fun = lambda x : VrfEq(x, verbose=verbose)
        Vrf = optimize.root_scalar(fun, x0=3, x1=5, xtol=1e-3, method='secant').root
        if verbose == True: print(f'Vrf = {VrfEq.root}')
        return Vrf
    
    Vrf = Vrf_()
    return Vrf




def UbiasCalc():
    """
    Function for calculating Ubias voltage & other electrical parameters.

    Returns
    -------
    None.
    Automatically collects all obtained parameters in out[] dict.

    """
    out['Vf'] = -(out['Te']/2)*np.log(gas_params[par['gas']]['M']/(par['m']*2*np.pi/(1.247**2)))        
    beta = 1
    out['sm2'] = (par['S1']+par['S2'])/(2*par['S2'])*np.power(par['S2']/par['S1'], beta/2-0.25)*par['sm']
    out['sm1'] = (par['S2']/par['S1'])**(1.5-beta)
    out['ns2'] = np.power(par['S1']/par['S2'], beta/2-0.25)
    out['ns1'] = np.power(par['S2']/par['S1'], beta-0.5)*out['ns2']
    out['Ji1'] = par['e']*out['ns1']*out['ub']
    out['Ji2'] = par['e']*out['ns2']*out['ub']
    out['C1'] = par['e0']*par['S1']/out['sm1']
    out['C2'] = par['e0']*par['S2']/out['sm2']
    out['Vp_amp'] = out['C1']/(out['C1']+out['C2'])*par['Vrf']
    out['Vp_avg'] = out['Vp_amp']+np.abs(out['Vf'])-out['Te']/2*np.log(2*np.pi*out['Vp_amp']/out['Te'])
    out['Vp_symm'] = par['Vrf']/2+np.abs(out['Vf'])-out['Te']/2*np.log(np.pi*par['Vrf']/out['Te'])
    out['Ubias'] = 2*out['Vp_avg']-par['Vrf']
    out['V1_avg'] = out['Vp_avg']-out['Ubias']
    out['V2_avg'] = out['Vp_avg']
    return None

def dECalc():
    """
    Function for calculating dE & other ion  parameters.

    Returns
    -------
    None.
    Automatically collects all obtained parameters in out[] dict.

    """
    out['dEi1'] = (4*np.sign(par['e']*out['V1_avg'])
                   *np.power(np.abs(par['e']*out['V1_avg']), (3/2))
                   /(par['omega']*out['sm1']*np.sqrt(2*gas_params[par['gas']]['M'])))
    out['dEi2'] = (4*np.sign(par['e']*out['V2_avg'])
                   *np.power(np.abs(par['e']*out['V2_avg']), (3/2))
                   /(par['omega']*out['sm2']*np.sqrt(2*gas_params[par['gas']]['M'])))
    out['dEi_symm'] = (4*np.sign(par['e']*out['Vp_symm'])
                       *np.power(np.abs(par['e']*out['Vp_symm']), (3/2))
                       /(par['omega']*out['sm2']*np.sqrt(2*gas_params[par['gas']]['M'])))
    out['Ns1'] = out['Ji1']/par['e']
    out['Ns_symm'] = out['Ji']/par['e']
    out['Ns2'] = out['Ji2']/par['e']
    return None

def main():
    ParLoader()
    out['Te'], out['ub'], out['Kiz'], out['Kel'], out['Kex'] = TeCalc(verbose=False)
    par['Vrf'] = VrfCalc(verbose=False)
    out['Vrf'] = par['Vrf']
    UbiasCalc()
    dECalc()

main()

    
print('Execution time: %s seconds' % (time.time() - start_time))