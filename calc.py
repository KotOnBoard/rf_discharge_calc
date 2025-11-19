#import time
import matplotlib.pyplot as plt
from scipy import optimize, integrate
import numpy as np
import json5
import sys
import csv
from scipy.interpolate import Akima1DInterpolator#interp1d
#start_time = time.time()


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
    par : dict
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
                'sm': 0.002,
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
            par['p'] = 2 #pressure
            par['l'] = 0.1 #electrode dist
            par['R'] = 0.314 #electrode diam
            par['f'] = 80
            #par['Vrf'] = 159.38543719750095
            par['Pwr'] = 75 #rf field params
            par['Assy'] = 2
            par['S1'] = 0.04625
        elif (type(vrname)==dict):
            par.update(vrname)
        elif (type(vrname)==str): 
            with open(f"conf/{prname}.json5", 'r') as filevr:
                par.update(json5.load(filevr))
        par['f']*=1e6 
        if 'S1' not in par:
            par['S1'] = np.pi*par['R']**2
        else: 
            par['R'] = np.sqrt(par['S1']/np.pi)
        if 'S2' not in par:
            par['S2'] = par['Assy']*par['S1'] #произвольная геометрическая ассиметрия 4
        if 'omega' not in par:
            par['omega'] = 2*np.pi*par['f']
        if 'd' not in par:
            par['d'] = par['l']-2*par['sm'] 
        if 'ng' not in par:
            par['ng'] = (2*par['p'])/(3*par['k']*par['Ti'])
        if 'lam_i' not in par:
            par['lam_i'] = 1/(par['ng']*gas_params[par['gas']]['sig_i']) 
        if 'lam_i_d' not in par:
            par['lam_i_d'] = par['lam_i']/par['d']
        if 'hl' not in par:
            par['hl'] = 0.86*(3+par['d']/(2*par['lam_i']))**(-1/2)
        if 'hR' not in par:
            par['hR'] = 0.8*(4+par['R']/par['lam_i'])**(-1/2)
        if 'deff' not in par:
            par['deff'] = par['l']/(2*par['hl']) #1/2*par['R']*par['l']/(par['R']*par['hl']+par['l']*par['hR'])
                    
        return par
    
    gas_params = LoadGasPar(gcname=gcname)
    par = LoadInitPar(prname=prname, vrname=vrname)
    
    return par, gas_params


"""
1. Блок расчёта Температуры электронов и констант скорости взаимодействия.
Нахождение энергии электронов через баланс ионизации и ухода через дрейф
"""
        
def TeCalc(par, gas_params, verbose=False):
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
    """
    def CrossectionInit(par):
        cont = [[],[]]
        with open(f'lib/cross_sections/Ciz;{par["gas"]}.csv', newline='') as f:
            reader = csv.reader(f)
            for k,row in enumerate(reader):
                if k<14: ...
                else: 
                    cont[0].append(float(row[0]))
                    cont[1].append(float(row[1]))
        cross = Akima1DInterpolator(cont[0], cont[1], extrapolate=True) #interp1d(cont[0], cont[1], fill_value=0)
        return cross
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
        return Kiz-ub*((par['ng']*par['deff'])**-1) #???
        #return par['d']*par['S1'] * par["ng"] * Kiz - (par['S1']+par['S2']) * ub #D.S. model

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
        """
        ???
        Kiz = 2.34e-14 * Te_val ** 0.59 * np.exp(-17.44 / Te_val) #D.S. model
        """
        #Legacy
        
        Kiz_integ = lambda v: (
            (gas_params[par['gas']]['aiz']*((par['m']*v**2/(2*par['e']))/gas_params[par['gas']]['xi_iz']-1))/
            ((par['m']*(v**2)/(2*par['e']))/gas_params[par['gas']]['xi_iz']+gas_params[par['gas']]['biz'])**gas_params[par['gas']]['ciz']
            *v*np.exp(-(par['m']*(v**2)/(2*par['e']*Te_val)))*4*np.pi*v**2
            )
        Kiz = (1e-20 * np.power((par['m']/(2*np.pi*par['e']*Te_val)), (3/2))*
        integrate.quad(Kiz_integ, (np.sqrt(2*par['e']*gas_params[par['gas']]['xi_iz']/par['m'])), 10**7)[0])
        
        """Kiz_integ = lambda v: (max(Csection(par['m']*np.power(v,2)/(2*par['e'])),0)*np.power(v,3)* #сечение должно быть от 
                               np.power(par['m']/(2*np.pi*par['e']*Te_val),3/2)*
                               np.exp(-par['m']*v**2/(2*par['e']*Te_val)))
        Kiz = (4*np.pi*integrate.quad(Kiz_integ, 0, 10**7)[0])"""
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
    
    #Csection = CrossectionInit(par)
    par['Te'], par['ub'] = Te_(verbose=verbose)
    par['Kiz'] = Kiz_(par['Te'], verbose=verbose)
    par['Kel'] = Kel_(par['Te'], verbose=verbose)
    par['Kex'] = Kex_(par['Te'], verbose=verbose)
    if (verbose==True): TeEval()
    return par
    

     

"""
Блок 2: Согласование мощности плазмы с целевой.
"""

def VrfCalc(par, gas_params, verbose=False):
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
        par['d'] = par['l']-2*par['sm']
        par['hl'] = 0.86*np.power(3+par['d']/(2*par['lam_i']),(-1/2))
        par['V1'] = np.abs(Vrf)/2
        par['Sohm'] = (1.73*par['m']*par['hl']/(2*par['e'])*par['e0']*par['omega']**2*par['vm']*
                np.sqrt(par['Te'])*np.sqrt(par['V1'])*par['d'])
        par['Sstoc'] = 0.45*np.sqrt(par['m']/par['e'])*par['e0']*par['omega']**2*np.sqrt(par['Te'])*par['V1']
        par['ns'] = ((par['Sohm']+2*par['Sstoc'])
        /(2*par['e']*par['ub']
          *(par['xi_c']+2*par['Te']+par['Te']
            *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
            +0.5*par['Te'])))
        par['n0'] = par['ns']/par['hl']
        par['V'] = 0.83*par['V1']
        par['Ji_symm'] = par['e']*par['ns']*par['ub']
        #print('Ji_symm =',par['Ji_symm'],' V^3/2 =',np.power(par['V'], 3/2))
        par['sm'] = np.sqrt(0.82*par['e0']*np.power(par['V'], 3/2)/par['Ji_symm']
                            *np.sqrt(2*par['e']/gas_params[par['gas']]['M']))
        #par['sm'] = 0.002
        par['J1'] = 1.23*par['omega']*par['e0']/par['sm']*par['V1']
        par['Sabs'] = (2*par['e']*par['ns']*par['ub']
                *(par['V']+par['xi_c']+2*par['Te']+par['Te']
                  *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
                  +0.5*par['Te'])*np.pi*par['R']**2)
        par['dPio'] = 2*par['e']*par['ns']*par['ub']*(gas_params[par['gas']]['xi_iz']+par['V'])/par['Sabs']
        #par['d'] = par['l']-2*par['sm']
        
        if (verbose==True) :
            check['d'].append(par['d'])
            check['hl'].append(par['hl'])
            check['Vrf'].append(Vrf)
            check['V1'].append(par['V1'])
            check['Sohm'].append(par['Sohm'])
            check['Sstoc'].append(par['Sstoc'])
            check['ns'].append(par['ns'])
            check['n0'].append(par['n0'])
            check['V'].append(par['V'])
            check['Ji'].append(par['Ji_symm'])
            check['sm'].append(par['sm'])
            check['J1'].append(par['J1'])
            check['Sabs'].append(par['Sabs'])
            print(f'Vrf = {Vrf}')
            print(f'Sabs = {par["Sabs"]}')
            print(f'ns = {par["ns"]}')
            print(f'Ji = {par["Ji_symm"]}')
        return par['Sabs'] - par['Pwr']

    def Vrf_(verbose=False):
        par["vm"] = par["Kel"]*par["ng"]
        par["xi_c"] = (1/par["Kiz"]*(par["Kiz"]*gas_params[par['gas']]['xi_iz']+
                              par["Kex"]*gas_params[par['gas']]['xi_ex']+
                              par["Kel"]*3*par["m"]*par["Te"]/gas_params[par['gas']]['M']))
        
        fun = lambda x : VrfEq(x, verbose=verbose)
        VrfSolve = optimize.root_scalar(fun, x0=500, x1=5000, xtol=1e-3, method='secant')
        
        if (verbose==True): print(f'VrfSolve = {VrfSolve}')
        return VrfSolve.root
    
    def VrfSample(Vrf, verbose=False, sm=0.001, d=par['l']-2*par['sm'], hl=0.86*(3+par['d']/(2*par['lam_i']))**(-1/2)):
        #d = par['l']-2*sm
        #hl = 0.86*np.power(np.abs(3+d)/(2*par['lam_i']),(-1/2))*np.sign(3+d)
        
        V1 = Vrf/2
        Sohm = (1.73*par['m']*hl/(2*par['e'])*par['e0']*par['omega']**2*par['vm']*
                np.sqrt(par['Te'])*np.sign(V1)*np.sqrt(np.abs(V1))*d)
        Sstoc = 0.45*np.sqrt(par['m']/par['e'])*par['e0']*par['omega']**2*np.sqrt(par['Te'])*V1
        ns = ((Sohm+2*Sstoc)
        /(2*par['e']*par['ub']
          *(par['xi_c']+2*par['Te']+par['Te']
            *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
            +0.5*par['Te'])))
        n0 = ns/hl
        V = 0.83*V1
        Ji = par['e']*ns*par['ub']
        sm = np.sqrt(0.82*par['e0']*np.power(np.abs(V), 3/2)/np.abs(Ji)
                     *np.sqrt(2*par['e']/gas_params[par['gas']]['M']))*np.sign(V)*np.sign(Ji)
        J1 = 1.23*par['omega']*par['e0']/sm*V1
        Sabs = (2*par['e']*ns*par['ub']
                *(V+par['xi_c']+2*par['Te']+par['Te']
                  *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
                  +0.5*par['Te'])*np.pi*par['R']**2)
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
    if (verbose==True): 
        VrfEval()
        VrfDynEval()

    return par

def SabsCalc(par, gas_params, verbose=True):
    par["vm"] = par["Kel"]*par["ng"]
    par["xi_c"] = (1/par["Kiz"]*(par["Kiz"]*gas_params[par['gas']]['xi_iz']+
                                 par["Kex"]*gas_params[par['gas']]['xi_ex']+
                                 par["Kel"]*3*par["m"]*par["Te"]/gas_params[par['gas']]['M']))
    par['d'] = par['l']-2*par['sm']
    par['hl'] = 0.86*np.power(3+par['d']/(2*par['lam_i']),(-1/2))
    par['V1'] = np.abs(par['Vrf'])/2
    par['Sohm'] = (1.73*par['m']*par['hl']/(2*par['e'])*par['e0']*par['omega']**2*par['vm']*
            np.sqrt(par['Te'])*np.sqrt(par['V1'])*par['d'])
    par['Sstoc'] = 0.45*np.sqrt(par['m']/par['e'])*par['e0']*par['omega']**2*np.sqrt(par['Te'])*par['V1']
    par['ns'] = ((par['Sohm']+2*par['Sstoc'])
    /(2*par['e']*par['ub']
      *(par['xi_c']+2*par['Te']+par['Te']
        *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
        +0.5*par['Te'])))
    par['n0'] = par['ns']/par['hl']
    par['V'] = 0.83*par['V1']
    par['Ji_symm'] = par['e']*par['ns']*par['ub']
    #print('Ji_symm =',par['Ji_symm'],' V^3/2 =',np.power(par['V'], 3/2))
    par['sm'] = np.sqrt(0.82*par['e0']*np.power(par['V'], 3/2)/par['Ji_symm']
                        *np.sqrt(2*par['e']/gas_params[par['gas']]['M']))
    #par['sm'] = 0.002
    par['J1'] = 1.23*par['omega']*par['e0']/par['sm']*par['V1']
    par['Sabs'] = (2*par['e']*par['ns']*par['ub']
            *(par['V']+par['xi_c']+2*par['Te']+par['Te']
              *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
              +0.5*par['Te'])*par['S1'])
    par['Pwr'] = par['Sabs']
    par['dPio'] = 2*par['e']*par['ns']*par['ub']*(gas_params[par['gas']]['xi_iz']+par['V'])/par['Sabs']
    
    return par



def UbiasCalc(par, gas_params,):
    """
    Function for calculating Ubias voltage & other electrical parameters.

    Returns
    -------
    None.
    Automatically collects all obtained parameters in par[] dict.

    """
    par['Vf'] = -(par['Te']/2)*np.log(gas_params[par['gas']]['M']/(par['m']*2*np.pi/(1.247**2)))        
    beta = 1
    par['sm2'] = (par['S1']+par['S2'])/(2*par['S2'])*np.power(par['S2']/par['S1'], beta/2-0.25)*par['sm']
    par['sm1'] = np.power(par['S2']/par['S1'], 1.5-beta)*par['sm2']
    par['ns2'] = np.power(par['S1']/par['S2'], beta/2-0.25)*par['ns']
    par['ns1'] = np.power(par['S2']/par['S1'], beta-0.5)*par['ns2']
    par['Ji1'] = par['e']*par['ns1']*par['ub']
    par['Ji2'] = par['e']*par['ns2']*par['ub']
    par['C1'] = par['e0']*par['S1']/par['sm1']
    par['C2'] = par['e0']*par['S2']/par['sm2']
    par['Vp_amp'] = par['C1']/(par['C1']+par['C2'])*par['Vrf']
    par['Vp_avg'] = par['Vp_amp']+np.abs(par['Vf'])-par['Te']/2*np.log(2*np.pi*par['Vp_amp']/par['Te'])
    par['Vp_symm'] = par['Vrf']/2+np.abs(par['Vf'])-par['Te']/2*np.log(np.pi*par['Vrf']/par['Te'])
    par['Ubias'] = 2*par['Vp_avg']-par['Vrf']
    par['V1_avg'] = par['Vp_avg']-par['Ubias'] #для симметричного случая
    par['V2_avg'] = par['Vp_avg']
    #par['V1_avg'] = np.power(par['S2']/par['S1'], 2.5-beta)*par['V2_avg'] #устаревший
    return par

def dECalc(par, gas_params):
    """
    Function for calculating dE & other ion  parameters.

    Returns
    -------
    None.
    Automatically collects all obtained parameters in par[] dict.

    """
    par['dEi1'] = np.abs(4*np.power(par['e']*np.abs(par['V1_avg']), 1.5)*np.sign(par['V1_avg']) #194.37009
                       /(par['omega']*par['sm1']*np.sqrt(2*gas_params[par['gas']]['M'])))
    par['dEi1 [eV]'] = par['dEi1']/par['e']
    par['dEi2'] = np.abs(4*np.power(par['e']*np.abs(par['V2_avg']), 1.5)*np.sign(par['V2_avg']) #18.75306
                       /(par['omega']*par['sm2']*np.sqrt(2*gas_params[par['gas']]['M'])))
    par['dEi2 [eV]'] = par['dEi2']/par['e']
    par['dEi_symm'] = np.abs(4*np.power(par['e']*np.abs(par['Vp_symm']), 1.5)*np.sign(par['Vp_symm']) #167.38546
                           /(par['omega']*par['sm2']*np.sqrt(2*gas_params[par['gas']]['M'])))
    par['dEi_symm [eV]'] = par['dEi_symm']/par['e']
    par['Ns1'] = par['Ji1']/par['e']
    par['Ns_symm'] = par['Ji_symm']/par['e']
    par['Ns2'] = par['Ji2']/par['e']
    return par

def FiENorm(par, gas_params, verbose=False):
    j_ = ['1','2','_symm','1_avg','2_avg','p_symm']
    el = ['Mo', 'B', 'Si']
    xi_hilo = 330
    for j in [0,1,2]:
        FiE = lambda x: (2*par[f'Ns{j_[j]}']/(par['omega']*par[f'dEi{j_[j]} [eV]'])
                            *np.power(1-4*np.power(x-par[f'V{j_[j+3]}'], 2)
                                             /np.power(par[f'dEi{j_[j]} [eV]'], 2), -0.5))

        par[f'NN{j_[j]}'] = integrate.quad(FiE, par[f'V{j_[j+3]}']-(par[f'dEi{j_[j]} [eV]']/2),
                                           par[f'V{j_[j+3]}']+(par[f'dEi{j_[j]} [eV]']/2))[0]
        
        for ind,i in enumerate(el):
            cont = [[],[]]
            with open(f'lib/sputtering_yeild/{i};{par["gas"]}.csv', newline='') as f:
                reader = csv.reader(f)
                for k,row in enumerate(reader):
                    if k<9: ...
                    else: 
                        cont[0].append(np.log10(float(row[0])))
                        cont[1].append(np.log10(float(row[1])))
            log_interp = Akima1DInterpolator(cont[0], cont[1], extrapolate=True) #interp1d(cont[0], cont[1], fill_value=0)
            interp = lambda y: np.power(10.0, log_interp(np.log10(y)))                                                                          
            """gam = lambda x: (gas_params[par['gas']][f'y0{i}']+2*gas_params[par['gas']][f'A{i}']/np.pi
                                *(gas_params[par['gas']][f'w{i}']/
                                  (4*np.power(x/par['e']-gas_params[par['gas']][f'xc{i}'], 2)
                                   +np.power(gas_params[par['gas']][f'w{i}'], 2))))"""
            
                    
            gam1_eff = lambda x: (interp(x)*FiE(x)/par[f'NN{j_[j]}'])

            if par[f'V{j_[j+3]}']>0 and par[f'V{j_[j+3]}']-(par[f'dEi{j_[j]} [eV]']/2)>0: 
                par[f'G{j_[j]}_{i}_eff'] = (integrate.quad(gam1_eff, par[f'V{j_[j+3]}']-(par[f'dEi{j_[j]} [eV]']/2),
                                                           par[f'V{j_[j+3]}']+(par[f'dEi{j_[j]} [eV]']/2))[0])
            else: par[f'G{j_[j]}_{i}_eff'] = 0
            """Проверка адекватности показаний. Всё более 1e6 и менее 10^-6 приравнивается к нулю,  т.к. не явялется адекватным"""
            if par[f'G{j_[j]}_{i}_eff']>1e6 or par[f'G{j_[j]}_{i}_eff']<1e-6: par[f'G{j_[j]}_{i}_eff'] = 0
            if (verbose==True):
                print(f"mode {j_[j]}_{i}  = {(par[f'V{j_[j+3]}'] > xi_hilo)}") 
                print(f"heaviside_{j_[j]}_{i} = {(par['e']*par[f'V{j_[j+3]}']+(par[f'dEi{j_[j]}']/2))-gas_params[par['gas']][f'xi_{i}']*par['e']}")
                print(f"gamma high {j_[j]}_{i} = {(gas_params[par['gas']][f'gam_{i}_hi'](par['e']*par[f'V{j_[j+3]}'])*FiE(par['e']*par[f'V{j_[j+3]}'])/par[f'NN{j_[j]}']*np.heaviside(par['e']*par[f'V{j_[j+3]}']-gas_params[par['gas']][f'xi_{i}']*par['e'], 0))}")
                print(f"gamma{j_[j]}_{i}_eff_Integ = {gam1_eff(par['e']*par[f'V{j_[j+3]}'])}")
                print(f"V{j_[j+3]} = {par[f'V{j_[j+3]}']}")
                print(f"FiE{j_[j]} = {FiE(par['e']*par[f'V{j_[j+3]}'])}")
                print(f"{'='*30}")
    for j in [0,1,2]:        
        par[f'd{j_[j]}_MoB'] = par[f'G{j_[j]}_Mo_eff']/(par[f'G{j_[j]}_B_eff']+par[f'G{j_[j]}_Mo_eff']+1e-10)
        
        for i in el:
            par[f'd{j_[j]}M_{i}'] = par[f'Ji{j_[j]}']/(par['e']*par['N_A']*1e3)*par[f'M_{i}']*par[f'G{j_[j]}_{i}_eff']
            par[f'd{j_[j]}H_{i}'] = par[f'd{j_[j]}M_{i}']/par[f'ro_{i}']*6e10
    
    return par
    
            


def main(vrname=False, prname=False, gcname=False, verbose=False):
        err = False
        try: 
            try: par, gas_params = ParLoader(prname=prname, vrname=vrname, gcname=gcname)
            except: 
                err = 1#"Unable to load initial parameters"
                #print(err)
                sys.exit(1)
            error_lock = True
            recur_count = 0
            while error_lock:
                temp_d = par['d']
                try: par.update(TeCalc(par, gas_params, verbose=verbose))
                except: 
                    err = 2#"Electron temperature and/or rate coefficients can't be calculated"
                    #print(err)
                    sys.exit(2)
                if 'Pwr' in par:
                    try:
                        par_ = VrfCalc(par, gas_params, verbose=verbose)
                        par.update(par_)
                    except: 
                        err = 3#"Can't calculate Vrf from given Pwr"
                        #print(err)
                        sys.exit(31)
                elif 'Vrf' in par:
                    try:
                        par_ = SabsCalc(par, gas_params)
                        par.update(par_)
                    except:
                        err = 3#"Can't calculate Sabs from given Vrf"
                        #print(err)
                        sys.exit(32)
                if not(par['d']>temp_d*1.1 or par['d']<temp_d*0.9):
                    error_lock = False
                    #print(f'+++d error is achieved. d = {par["d"]}.+++')
                elif recur_count>20:
                    error_lock = False
                else:
                    #print(f'---d error is {par["d"]-temp_d}. Recalculating Te & Vrf...---')
                    recur_count+=1
            try:
                par.update(UbiasCalc(par, gas_params))
            except:
                err = 4#"Can't calculate Ubias and/or other circuit parameters"
                #print(err)
                sys.exit(4)
            
            try:
                par.update(dECalc(par, gas_params))
            except:
                err = 5#"Can't calculate dEi"
                #print(err)
                sys.exit(5)
            try:
                par.update(FiENorm(par, gas_params, verbose=verbose))
            except:
                err = 6#"Unable to find etching params"
                #print(err)
                sys.exit(6)
               
        finally: 
            par['mi']=gas_params[par['gas']]['M']
            par['f']=par['f']/1e6
            return par

#par = main(verbose=False)
#print('Execution time: %s seconds' % (time.time() - start_time))
#sys.exit(0)