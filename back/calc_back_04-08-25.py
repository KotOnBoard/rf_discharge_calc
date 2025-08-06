import time
import typing
import matplotlib.pyplot as plt
from scipy import optimize, integrate
import numpy as np
import json5
import sys
start_time = time.time()


def ParLoader(prname=False, vrname=False, gcname=False):
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
                       "ciz": 2,
                       "xi_Mo": 25,
                       "xi_0": 25,
                       "y00": 0.336,
                       "xc0": 25,
                       "w0": 344.90083,
                       "A0": -182.02909,
                       "xi_B": 52,
                       "xi_1": 52,
                       "y01": 0.85,
                       "xc1": 52,
                       "w1": 475.94718,
                       "A1": -635.40798,
                       "xi_Si": 1,
                       "xi_2": 1,
                       "y02": 1,
                       "xc2": 1,
                       "w2": 1,
                       "A2": 1,
                       "xi_MoO3": 25,
                       "xi_3": 40,
                       "y03": 0.2262,
                       "xc3": 40,
                       "w3": 318.32173,
                       "A3": -113.01079,
                       "gam_Mo_hi": lambda x: 0.8236-0.80571*np.exp(-x/(1.6e-19*1039.32852)),
                       "gam_B_hi": lambda x: 1.85289-1.9649*np.exp(-x/(1.6e-19*952.03776)),
                       "gam_MoO3_hi": lambda x: 0.15682-0.16238*np.exp(-x/(1.6e-19*1310.76845))
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
                       "ciz": 1.91,
                       "xi_Mo": 61,
                       "xi_0": 61,
                       "y00": 0.029,
                       "xc0": 61,
                       "w0": 344.15473,
                       "A0": -15.67888,
                       "xi_B": 20,
                       "xi_1": 20,
                       "y01": 0.17246,
                       "xc1": 20,
                       "w1": 177.45386,
                       "A1": -47.93846,
                       "gam_Mo_hi": lambda x: 0.05636-0.06908*np.exp(-x/(1.6e-19*451.8792)),
                       "gam_B_hi": lambda x: (0.63346-2*11205.55135*13921.51887
                                              /(np.pi*(4*np.power(x/1.6e-19-4440, 2)
                                                       +np.power(13921.51887, 2)))
                                              -2*160.27248*348.90715
                                              /(np.pi*(4*np.power(x/1.6e-19+28, 2)
                                                       +np.power(348.90715, 2)))),
                       "gam_MoO3_hi": lambda x: 0.03212-0.03677*np.exp(-x/(1.6e-19*767.13872)),
                       "xi_MoO3": 80,
                       "xi_3": 80
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
                       "ciz": 1.86,
                       "xi_Mo": 25,
                       "xi_0": 25,
                       "y00": 0.313,
                       "xc0": 25,
                       "w0": 251.67992,
                       "A0": -123.72542,
                       "xi_B": 90,
                       "xi_1": 90,
                       "y01": 0.85,
                       "xc1": 90,
                       "w1": 748.99493,
                       "A1": -1000,
                       "xi_Si": 31,
                       "xi_2": 31,
                       "y02": 0.35,
                       "xc2": 31,
                       "w2": 333.74389,
                       "A2": -183.47724,
                       "xi_MoO3": 25,
                       "xi_3": 25,
                       "y03": 0.2262,
                       "xc3": 40,
                       "w3": 318.32173,
                       "A3": -113.01079,
                       "gam_Mo_hi": lambda x: 2.74606-2.97354*np.exp(-x/(1.6e-19*1687.79318)),
                       "gam_B_hi": lambda x: 1.64807-1.62468*np.exp(-x/(1.6e-19*1368.54118)),
                       "gam_MoO3_hi": lambda x: 0.15682-0.16238*np.exp(-x/(1.6e-19*1310.76845))
                       }
                   }
        else: 
            with open(f"conf/{gcname}.json5", 'r') as file:
                gas_params = json5.load(file)
        return gas_params
    
    def LoadInitPar(prname=False, vrname=False):
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
                'ro_B': 2.34e3,
                'M_B': 10.81,
                'ro_Mo': 10.2e3,
                'M_Mo': 95.95,
                'ro_MoO3': 4.9e3,
                'M_MoO3': 95.95+3*16,
                'e': 1.6e-19, # Заряд электрона [Кл]
                'sm': 0.009,
                'k': 1.380649e-23,   # Постоянная Больцмана [Дж/К]
                'e0': 8.85E-12,  # Диэлектрическая постоянная [Ф/м]
                'N_A': 6.022e23,
                'Ti': 300,
                'Tn': 300,
                'm': 9.11e-31,  # Масса электрона [кг]
                }
        else:
            with open(f"conf/{prname}.json5", 'r') as filepr:
                par = json5.load(filepr)
        if (vrname==False):
            par['gas'] = "Ar"
            par['p'] = 1 #pressure
            par['l'] = 0.07 #electrode dist
            par['R'] = 0.11767 #electrode diam
            par['f'] = 80
            #par['Vrf'] = 159.38543719750095
            par['Pwr'] = 75 #rf field params
            par['Assy'] = 2
        elif (type(vrname)==dict):
            par.update(vrname)
        elif (type(vrname)==str): 
            with open(f"conf/{prname}.json5", 'r') as filevr:
                par.update(json5.load(filevr))
        par['f']*=1e6   
        try:
            par['S1']=par['S1']
        except:
            par['S1'] = np.pi*par['R']**2
        try:
            par['S2']=par['S2']
        except:
            par['S2'] = par['Assy']*par['S1'] #произвольная геометрическая ассиметрия 4
        try:
            par['omega']=par['omega']
        except:
            par['omega'] = 2*np.pi*par['f']
        try:
            par['d']=par['d']
        except:
            par['d'] = par['l']-2*par['sm'] 
        try:
            par['ng']=par['ng']
        except:
            par['ng'] = (2*par['p'])/(3*par['k']*par['Ti'])
        try:
            par['lam_i']=par['lam_i']
        except:
            par['lam_i'] = 1/(par['ng']*gas_params[par['gas']]['sig_i']) 
        try:
            par['lam_i_d']=par['lam_i_d']
        except:
            par['lam_i_d'] = par['lam_i']/par['d']
        try:
            par['hl']=par['hl']
        except:
            par['hl'] = 0.86*(3+par['d']/(2*par['lam_i']))**(-1/2)
        try:
            par['hR']=par['hR']
        except:
            par['hR'] = 0.8*(4+par['R']/par['lam_i'])**(-1/2)
        try:
            par['deff']=par['deff']
        except:
            par['deff'] = par['l']/(2*par['hl']) #1/2*par['R']*par['l']/(par['R']*par['hl']+par['l']*par['hR'])
                    
        return par
    
    gas_params = LoadGasPar(gcname=gcname)
    par = LoadInitPar(prname=prname, vrname=vrname)
    out = {}
    
    return par, gas_params, out


"""
1. Блок расчёта Температуры электронов и констант скорости взаимодействия.
Нахождение энергии электронов через баланс ионизации и ухода через дрейф
"""
        
def TeCalc(par, gas_params, out, verbose=False):
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
        Bohm velocity.
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
            Bohm velocity.

        """
        fun = lambda x: TeEquation(x, verbose=verbose)
        Te_val = optimize.root_scalar(fun, bracket=[1e-3, 30], x0=3, x1=5, xtol=1e-3, method='brentq')
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
        (((gas_params[par['gas']]['a1']+gas_params[par['gas']]['b1']*np.power(par['m']*v**2/(2*par['e']), gas_params[par['gas']]['c1']))/
        (1+gas_params[par['gas']]['d1']*np.power(par['m']*v**2/(2*par['e']),gas_params[par['gas']]['e1'])))+
        ((gas_params[par['gas']]['a2']+gas_params[par['gas']]['b2']*np.power(par['m']*v**2/(2*par['e']), gas_params[par['gas']]['c2']))/
        (1+gas_params[par['gas']]['d2']*np.power(par['m']*v**2/(2*par['e']), gas_params[par['gas']]['e2']))))
        *v*np.exp(-(par['m']*(v**2)/(2*par['e'])/Te_val))*4*np.pi*v**2
        )
        Kel = (1e-20*np.power(par['m']/(2*np.pi*par['e']*Te_val), 3/2)*
        integrate.quad(Kel_integ, 0, 1e7)[0])
        if verbose == True: 
            print(f"Kel = {Kel}")
            print(f'Kel_Integ = {integrate.quad(Kel_integ, 0, 1e7)}')
            v=1e7
            print(f"InInteg = {(1+gas_params[par['gas']]['d2']*np.power(par['m']*v**2/(2*par['e']), gas_params[par['gas']]['e2']))}")
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
            *v*np.exp(-(par['m']*(v**2)/(2*par['e']*Te_val)))*4*np.pi*v**2
            )    
        Kex = (1e-20*np.power((par['m']/(2*np.pi*par['e']*Te_val)), (3/2))*
        integrate.quad(Kex_integ, (np.sqrt(2*par['e']*gas_params[par['gas']]['xi_ex']/par['m'])), 10**7)[0])
        if verbose == True: print(f"Kex = {Kex}")
        return Kex
    
    def TeEval(a=1e-2, b=10, c=1e-2, d=2e-14):
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
        plt.yscale('log')
        plt.ylim(1e-18,1e-14)
        plt.legend()
        plt.show()
    
    out['Te'], out['ub'] = Te_(verbose=verbose)
    out['Kiz'] = Kiz_(out['Te'], verbose=verbose)
    out['Kel'] = Kel_(out['Te'], verbose=verbose)
    out['Kex'] = Kex_(out['Te'], verbose=verbose)
    if (verbose==True): TeEval()
    return out
    

     

"""
Блок 2: Согласование мощности плазмы с целевой.
"""

def VrfCalc(par, gas_params, out, verbose=False):
    """
    Function for calculating absorbed power and matching it by varying RF voltage. 

    Returns
    -------
    Vrf : float
        RF voltage with matching RF power & Power absobed by plasma.

    """
    if (verbose==True): 
        check = {}
        check['d'] = []
        check['hl'] = []
        check['Vrf'] = []
        check['V1'] = []
        check['Sohm'] = []
        check['Sstoc'] = []
        check['ns'] = []
        check['n0'] = []
        check['V'] = []
        check['Ji'] = []
        check['sm'] = []
        check['J1'] = []
        check['Sabs'] = []
    
             
    def VrfEq(Vrf, verbose=False):
        #par['d'] = par['l']-2*par['sm']
        #par['hl'] = 0.86*np.power(np.abs(3+par['d'])/(2*par['lam_i']),(-1/2))*np.sign(3+par['d'])
        
        out['V1'] = np.abs(Vrf)/2
        out['Sohm'] = (1.73*par['m']*par['hl']/(2*par['e'])*par['e0']*par['omega']**2*out['vm']*
                np.sqrt(out['Te'])*np.sqrt(out['V1'])*par['d'])
        out['Sstoc'] = 0.45*np.sqrt(par['m']/par['e'])*par['e0']*par['omega']**2*np.sqrt(out['Te'])*out['V1']
        out['ns'] = ((out['Sohm']+2*out['Sstoc'])
        /(2*par['e']*out['ub']
          *(out['xi_c']+2*out['Te']+out['Te']
            *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
            +0.5*out['Te'])))
        out['n0'] = out['ns']/par['hl']
        out['V'] = 0.83*out['V1']
        out['Ji_symm'] = par['e']*out['ns']*out['ub']
        par['sm'] = np.sqrt(0.82*par['e0']*np.power(out['V'], 3/2)/out['Ji_symm']
                            *np.sqrt(2*par['e']/gas_params[par['gas']]['M']))
        out['J1'] = 1.23*par['omega']*par['e0']/par['sm']*out['V1']
        out['Sabs'] = (2*par['e']*out['ns']*out['ub']
                *(out['V']+out['xi_c']+2*out['Te']+out['Te']
                  *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
                  +0.5*out['Te'])*np.pi*par['R']**2)
        out['dPio'] = 2*par['e']*out['ns']*out['ub']*(gas_params[par['gas']]['xi_iz']+out['V'])/out['Sabs']
        
        if (verbose==True) :
            check['d'].append(par['d'])
            check['hl'].append(par['hl'])
            check['Vrf'].append(Vrf)
            check['V1'].append(out['V1'])
            check['Sohm'].append(out['Sohm'])
            check['Sstoc'].append(out['Sstoc'])
            check['ns'].append(out['ns'])
            check['n0'].append(out['n0'])
            check['V'].append(out['V'])
            check['Ji'].append(out['Ji_symm'])
            check['sm'].append(par['sm'])
            check['J1'].append(out['J1'])
            check['Sabs'].append(out['Sabs'])
            print(f'Vrf = {Vrf}')
            print(f'Sabs = {out["Sabs"]}')
            print(f'ns = {out["ns"]}')
            print(f'Ji = {out["Ji_symm"]}')
        return out['Sabs'] - par['Pwr']

    def Vrf_(verbose=False):
        out["vm"] = out["Kel"]*par["ng"]
        out["xi_c"] = (1/out["Kiz"]*(out["Kiz"]*gas_params[par['gas']]['xi_iz']+
                              out["Kex"]*gas_params[par['gas']]['xi_ex']+
                              out["Kel"]*3*par["m"]*out["Te"]/gas_params[par['gas']]['M']))
        
        fun = lambda x : VrfEq(x, verbose=verbose)
        VrfSolve = optimize.root_scalar(fun, bracket=[1e0, 1e10], x0=500, x1=5000, xtol=1e-3, method='brentq')
        
        if (verbose==True): print(f'VrfSolve = {VrfSolve}')
        return VrfSolve.root
    
    def VrfSample(Vrf, verbose=False, sm=0.001, d=par['l']-2*par['sm'], hl=0.86*(3+par['d']/(2*par['lam_i']))**(-1/2)):
        #d = par['l']-2*sm
        #hl = 0.86*np.power(np.abs(3+d)/(2*par['lam_i']),(-1/2))*np.sign(3+d)
        
        V1 = Vrf/2
        Sohm = (1.73*par['m']*hl/(2*par['e'])*par['e0']*par['omega']**2*out['vm']*
                np.sqrt(out['Te'])*np.sign(V1)*np.sqrt(np.abs(V1))*d)
        Sstoc = 0.45*np.sqrt(par['m']/par['e'])*par['e0']*par['omega']**2*np.sqrt(out['Te'])*V1
        ns = ((Sohm+2*Sstoc)
        /(2*par['e']*out['ub']
          *(out['xi_c']+2*out['Te']+out['Te']
            *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
            +0.5*out['Te'])))
        n0 = ns/hl
        V = 0.83*V1
        Ji = par['e']*ns*out['ub']
        sm = np.sqrt(0.82*par['e0']*np.power(np.abs(V), 3/2)/np.abs(Ji)
                     *np.sqrt(2*par['e']/gas_params[par['gas']]['M']))*np.sign(V)*np.sign(Ji)
        J1 = 1.23*par['omega']*par['e0']/sm*V1
        Sabs = (2*par['e']*ns*out['ub']
                *(V+out['xi_c']+2*out['Te']+out['Te']
                  *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
                  +0.5*out['Te'])*np.pi*par['R']**2)
        return d, hl, Sohm, Sstoc, ns, Ji, sm, J1, Sabs, n0
   
    def VrfEval(a=1, b=1e4+1, c=10, d=500):
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
        range_ = np.arange(a, b, c)
        plt.figure()
        plt.plot(range_, [par['Pwr']]*(len(range_)), label='Pwr')
        plt.plot(range_, [VrfSample(Vrf, verbose=verbose)[0] for Vrf in range_], label='d')
        plt.plot(range_, [VrfSample(Vrf, verbose=verbose)[1] for Vrf in range_], label='hl')
        plt.plot(range_, [VrfSample(Vrf, verbose=verbose)[2] for Vrf in range_], label='Sohm')
        plt.plot(range_, [VrfSample(Vrf, verbose=verbose)[3] for Vrf in range_], label='Sstoc')
        plt.plot(range_, [VrfSample(Vrf, verbose=verbose)[4] for Vrf in range_], label='ns')
        plt.plot(range_, [VrfSample(Vrf, verbose=verbose)[5] for Vrf in range_], label='Ji')
        plt.plot(range_, [VrfSample(Vrf, verbose=verbose)[6] for Vrf in range_], label='sm')
        plt.plot(range_, [VrfSample(Vrf, verbose=verbose)[7] for Vrf in range_], label='J1')
        plt.plot(range_, [VrfSample(Vrf, verbose=verbose)[8] for Vrf in range_], label='Sabs')
        #plt.yscale('log')
        plt.xlim(Vrf-20, Vrf+20)
        plt.ylim(par['Pwr']-0.5,par['Pwr']+0.5)
        plt.legend()
        plt.show()
        
    def VrfDynEval():
        range_ = [i for i in range(len(check['d']))]
        plt.figure()
        plt.plot(range_, check['d'], label='d', color='tab:blue')
        plt.plot(range_, check['hl'], label='hl', color='tab:orange')
        plt.plot(range_, check['Vrf'], label='Vrf', color='tab:green')
        plt.plot(range_, check['V1'], label='V1', color='tab:red')
        plt.plot(range_, check['Sohm'], label='Sohm', color='tab:purple')
        plt.plot(range_, check['Sstoc'], label='Sstoc', color='tab:brown')
        plt.plot(range_, check['ns'], label='ns', color='tab:pink')
        plt.plot(range_, check['n0'], label='n0', color='tab:gray')
        plt.plot(range_, check['V'], label='V', color='tab:olive')
        plt.plot(range_, check['Ji'], label='Ji', color='tab:cyan')
        plt.plot(range_, check['sm'], label='sm', color='darkblue')
        plt.plot(range_, check['J1'], label='J1', color='orangered')
        plt.plot(range_, check['Sabs'], label='Sabs', color='darkred')
        plt.yscale('log')
        #plt.xlim(0,5)
        #plt.ylim(0,0.2)
        plt.legend()
        plt.show()  
        
    Vrf = Vrf_(verbose=verbose) 
    par['Vrf'] = Vrf
    out['Vrf'] = Vrf
    if (verbose==True): 
        VrfEval()
        VrfDynEval()

    return par, out

def SabsCalc(par, gas_params, out, verbose=False):
    
    out["vm"] = out["Kel"]*par["ng"]
    out["xi_c"] = (1/out["Kiz"]*(out["Kiz"]*gas_params[par['gas']]['xi_iz']+
                          out["Kex"]*gas_params[par['gas']]['xi_ex']+
                          out["Kel"]*3*par["m"]*out["Te"]/gas_params[par['gas']]['M']))
    out['V1'] = par['Vrf']/2
    out['Sohm'] = (1.73*par['m']*par['hl']/(2*par['e'])*par['e0']*par['omega']**2*out['vm']*
            np.sqrt(out['Te'])*np.sqrt(out['V1'])*par['d'])
    out['Sstoc'] = 0.45*np.sqrt(par['m']/par['e'])*par['e0']*par['omega']**2*np.sqrt(out['Te'])*out['V1']
    out['ns'] = ((out['Sohm']+2*out['Sstoc'])
    /(2*par['e']*out['ub']
      *(out['xi_c']+2*out['Te']+out['Te']
        *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
        +0.5*out['Te'])))
    out['n0'] = out['ns']/par['hl']
    out['V'] = 0.83*out['V1']
    out['Ji_symm'] = par['e']*out['ns']*out['ub']
    par['sm'] = np.sqrt(0.82*par['e0']*np.power(out['V'], 3/2)/out['Ji_symm']
                        *np.sqrt(2*par['e']/gas_params[par['gas']]['M']))
    out['J1'] = 1.23*par['omega']*par['e0']/par['sm']*out['V1']
    out['Sabs'] = (2*par['e']*out['ns']*out['ub']
            *(out['V']+out['xi_c']+2*out['Te']+out['Te']
              *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
              +0.5*out['Te'])*np.pi*par['R']**2)
    out['dPio'] = 2*par['e']*out['ns']*out['ub']*(gas_params[par['gas']]['xi_iz']+out['V'])/out['Sabs']

    return par, out



def UbiasCalc(par, gas_params, out, ):
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
    out['sm1'] = np.power(par['S2']/par['S1'], 1.5-beta)*out['sm2']
    out['ns2'] = np.power(par['S1']/par['S2'], beta/2-0.25)*out['ns']
    out['ns1'] = np.power(par['S2']/par['S1'], beta-0.5)*out['ns2']
    out['Ji1'] = par['e']*out['ns1']*out['ub']
    out['Ji2'] = par['e']*out['ns2']*out['ub']
    out['C1'] = par['e0']*par['S1']/out['sm1']
    out['C2'] = par['e0']*par['S2']/out['sm2']
    out['Vp_amp'] = out['C1']/(out['C1']+out['C2'])*out['Vrf']
    out['Vp_avg'] = out['Vp_amp']+np.abs(out['Vf'])-out['Te']/2*np.log(2*np.pi*out['Vp_amp']/out['Te'])
    out['Vp_symm'] = par['Vrf']/2+np.abs(out['Vf'])-out['Te']/2*np.log(np.pi*par['Vrf']/out['Te'])
    out['Ubias'] = 2*out['Vp_avg']-par['Vrf']
    out['V1_avg'] = out['Vp_avg']-out['Ubias'] #для симметричного случая
    out['V2_avg'] = out['Vp_avg']
    #out['V1_avg'] = np.power(par['S2']/par['S1'], 2.5-beta)*out['V2_avg'] #устаревший
    return out

def dECalc(par, gas_params, out):
    """
    Function for calculating dE & other ion  parameters.

    Returns
    -------
    None.
    Automatically collects all obtained parameters in out[] dict.

    """

    out['dEi1'] = (4*np.power(par['e']*np.abs(out['V1_avg']), 1.5)*np.sign(out['V1_avg']) #194.37009
                   /(par['omega']*out['sm1']*np.sqrt(2*gas_params[par['gas']]['M'])))
    out['dEi1 [eV]'] = out['dEi1']/par['e']
    out['dEi2'] = (4*np.power(par['e']*np.abs(out['V2_avg']), 1.5)*np.sign(out['V2_avg']) #18.75306
                   /(par['omega']*out['sm2']*np.sqrt(2*gas_params[par['gas']]['M'])))
    out['dEi2 [eV]'] = out['dEi2']/par['e']
    out['dEi_symm'] = (4*np.power(par['e']*np.abs(out['Vp_symm']), 1.5)*np.sign(out['Vp_symm']) #167.38546
                       /(par['omega']*out['sm2']*np.sqrt(2*gas_params[par['gas']]['M'])))
    out['dEi_symm [eV]'] = out['dEi_symm']/par['e']
    out['Ns1'] = out['Ji1']/par['e']
    out['Ns_symm'] = out['Ji_symm']/par['e']
    out['Ns2'] = out['Ji2']/par['e']
    return out

def FiENorm(par, gas_params, out, verbose=False):
    
    j_ = ['1','2','_symm','1_avg','2_avg','p_symm']
    el =  ['Mo', 'B', 'Si', 'MoO3']
    xi_hilo = 330
    for j in [0,1,2]:
        FiE = lambda x: (2*out[f'Ns{j_[j]}']/(par['omega']*out[f'dEi{j_[j]}'])
                            *np.power(1-4*np.power(x-par['e']*out[f'V{j_[j+3]}'], 2)
                                             /np.power(out[f'dEi{j_[j]}'], 2), -0.5))

        out[f'NN{j_[j]}'] = integrate.quad(FiE, par['e']*out[f'V{j_[j+3]}']-(out[f'dEi{j_[j]}']/2),
                                           par['e']*out[f'V{j_[j+3]}']+(out[f'dEi{j_[j]}']/2))[0]
                                                                                              
        for i in [0,1,3]:
            gam = lambda x: (gas_params[par['gas']][f'y0{i}']+2*gas_params[par['gas']][f'A{i}']/np.pi
                                *(gas_params[par['gas']][f'w{i}']/
                                  (4*np.power(x/par['e']-gas_params[par['gas']][f'xc{i}'], 2)
                                   +np.power(gas_params[par['gas']][f'w{i}'], 2))))
            
            if ((out[f'V{j_[j+3]}'] < xi_hilo) and (i < 3)): 
                gam1_eff = lambda x: (gam(x)*FiE(x)/out[f'NN{j_[j]}']*np.heaviside(x-gas_params[par['gas']][f'xi_{i}']*par['e'], 0))
            else: 
                gam1_eff = lambda x: (gas_params[par['gas']][f'gam_{el[i]}_hi'](x)*FiE(x)/out[f'NN{j_[j]}']*np.heaviside(x-gas_params[par['gas']][f'xi_{i}']*par['e'], 0))
            out[f'Г{j_[j]}_{el[i]}_eff'] = (integrate.quad(gam1_eff, par['e']*out[f'V{j_[j+3]}']-(out[f'dEi{j_[j]}']/2),
                                                      par['e']*out[f'V{j_[j+3]}']+(out[f'dEi{j_[j]}']/2))[0])
            if (verbose==True):
                print(f"mode {j_[j]}_{el[i]}  = {(out[f'V{j_[j+3]}'] > xi_hilo)}") 
                try: 
                    print(f"gamma low {j_[j]}_{el[i]} = {gam(par['e']*out[f'V{j_[j+3]}'])}")
                except: ...
                print(f"heaviside_{j_[j]}_{el[i]} = {(par['e']*out[f'V{j_[j+3]}']+(out[f'dEi{j_[j]}']/2))-gas_params[par['gas']][f'xi_{i}']*par['e']}")
                print(f"gamma high {j_[j]}_{el[i]} = {(gas_params[par['gas']][f'gam_{el[i]}_hi'](par['e']*out[f'V{j_[j+3]}'])*FiE(par['e']*out[f'V{j_[j+3]}'])/out[f'NN{j_[j]}']*np.heaviside(par['e']*out[f'V{j_[j+3]}']-gas_params[par['gas']][f'xi_{i}']*par['e'], 0))}")
                print(f"gamma{j_[j]}_{el[i]}_eff_Integ = {gam1_eff(par['e']*out[f'V{j_[j+3]}'])}")
                print(f"V{j_[j+3]} = {out[f'V{j_[j+3]}']}")
                print(f"FiE{j_[j]} = {FiE(par['e']*out[f'V{j_[j+3]}'])}")
                print(f"{'='*30}")
    for j in [0,1,2]:        
        out[f'd{j_[j]}_MoB'] = out[f'Г{j_[j]}_Mo_eff']/(out[f'Г{j_[j]}_B_eff']+out[f'Г{j_[j]}_Mo_eff']+1e-10)
        out[f'd{j_[j]}_MoO3'] = out[f'Г{j_[j]}_Mo_eff']/(out[f'Г{j_[j]}_MoO3_eff']+out[f'Г{j_[j]}_Mo_eff']+1e-10)
        
        for i in [0,1,3]:
            out[f'd{j_[j]}M_{el[i]}'] = out[f'Ji{j_[j]}']/(par['e']*par['N_A']*1e3)*par[f'M_{el[i]}']*out[f'Г{j_[j]}_{el[i]}_eff']
            out[f'd{j_[j]}H_{el[i]}'] = out[f'd{j_[j]}M_{el[i]}']/par[f'ro_{el[i]}']*6e10
    
    return out
    
            


def calc(prname=False, vrname=False, gcname=False, verbose=False):
    err = False
    try: 
        try: par, gas_params, out = ParLoader(prname=prname, vrname=vrname, gcname=gcname)
        except: 
            err = f"{'⚠'*20}\nUnable to load initial parameters\n{'⚠'*20}"
            print(err)
            sys.exit(1)
        try: out.update(TeCalc(par, gas_params, out, verbose=verbose))
        except: 
            err = f"{'⚠'*20}\nElectron temperature and/or rate coefficients can't be calculated\n{'⚠'*20}"
            print(err)
            sys.exit(2)
        if (par['Pwr']):
            try:
                par_, out_ = VrfCalc(par, gas_params, out, verbose=verbose)
                par.update(par_)
                out.update(out_)
            except: 
                err = f"{'⚠'*20}\nCan't calculate Vrf from given Pwr\n{'⚠'*20}"
                print(err)
                sys.exit(31)
        elif (par['Vrf']):
            try:
                par_, out_ = SabsCalc(par, gas_params, out)
                par.update(par_)
                out.update(out_)
            except:
                err = f"{'⚠'*20}\nCan't calculate Sabs from given Vrf\n{'⚠'*20}"
                print(err)
                sys.exit(32)
        try:
            out.update(UbiasCalc(par, gas_params, out))
        except:
            err = f"{'⚠'*20}\nCan't calculate Ubias and/or other circuit parameters\n{'⚠'*20}"
            print(err)
            sys.exit(4)
        try:
            out.update(dECalc(par, gas_params, out))
        except:
            err = f"{'⚠'*20}\nCan't calculate dEi\n{'⚠'*20}"
            print(err)
            sys.exit(5)
        try:
            out.update(FiENorm(par, gas_params, out, verbose=verbose))
        except:
            err = f"{'⚠'*20}\nUnable to find etching params\n{'⚠'*20}"
            print(err)
            sys.exit(6)
    finally: 
        return par, out, err

#par, out, err = calc(verbose=False)
#print('Execution time: %s seconds' % (time.time() - start_time))
#sys.exit(0)