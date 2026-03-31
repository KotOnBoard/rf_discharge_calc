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
                'M_Si': 28.02,
                'ro_Si': 2330,
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
            par['gas'] = 'He'
            par['l'] = 0.1
            par['f'] = 40
            par['Pwr'] = 200
            par['Assy'] = 8
            par['p'] = 3
            par['S1'] = 0.01
            """
            par['gas'] = "Ne"
            par['p'] = 6 #pressure
            par['l'] = 0.1 #electrode dist
            par['R'] = 0.0665 #electrode diam
            par['f'] = 10e7/1e6
            #par['Vrf'] = 159.38543719750095
            par['Pwr'] = 300 #rf field params
            par['Assy'] = 4
            #par['S1'] = 0.04625"""
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
            par['omega'] = 2*np.pi*par['f'] #циклическая частота тока
        if 'd' not in par:
            par['d'] = par['l']-2*par['sm'] #толщина булки
        if 'ng' not in par:
            par['ng'] = (2*par['p'])/(3*par['k']*par['Ti']) #Концентрация атомов газа
        if 'lam_i' not in par:
            par['lam_i'] = 1/(par['ng']*gas_params[par['gas']]['sig_i']) #mean free path
        if 'lam_i_d' not in par:
            par['lam_i_d'] = par['lam_i']/par['d'] #критерий существования столкновений
        if 'hl' not in par:
            par['hl'] = 0.86*(3+par['d']/(2*par['lam_i']))**(-1/2) #axial sheath edge pg. 331 
        if 'hR' not in par:
            par['hR'] = 0.8*(4+par['R']/par['lam_i'])**(-1/2) #radial sheath edge
        if 'deff' not in par:
            par['deff'] = par['l']/(2*par['hl']) #1/2*par['R']*par['l']/(par['R']*par['hl']+par['l']*par['hR'])
            #эффективная площадь плазмы 
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
        return Kiz-ub*((par['ng']*par['deff'])**-1) #pg 408 10.2.10
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

def VrfCalc(par, gas_params, mode='Pwr', verbose=False):
    """
    Function for calculating absorbed power and matching it by varying RF voltage. 

    Returns
    -------
    Vrf : float
        RF voltage with matching RF power & Power absobed by plasma.

    """

    def VrfEq(Vrf, verbose=False):
        #print("==============================\n par['Vrf'] = ",Vrf)
        par["vm"] = par["Kel"]*par["ng"]
        par["xi_c"] = (1/par["Kiz"]*(par["Kiz"]*gas_params[par['gas']]['xi_iz']+
                              par["Kex"]*gas_params[par['gas']]['xi_ex']+
                              par["Kel"]*3*par["m"]*par["Te"]/gas_params[par['gas']]['M']))  
        par['Vf'] = -(par['Te']/2)*np.log(gas_params[par['gas']]['M']/(par['m']*2*np.pi
                                                                       /(np.pow(1.247,2))))
        par['sm2'] = ((par['S1']+par['S2'])/(2*par['S2'])
                      *np.pow(par['S2']/par['S1'],1/2-0.25)*par['sm'])
        par['sm1'] = np.pow(par['S2']/par['S1'],1.5-1)*par['sm2']
        par['d'] = par['l']-par['sm1']-par['sm2']
        par['hl'] = 0.86*np.power(3+par['d']/(2*par['lam_i']),(-1/2))
        par['C1'] = par['e0']*par['S1']/par['sm1']
        par['C2'] = par['e0']*par['S2']/par['sm2']
        par['Vp_amp'] = par['C1']/(par['C1']+par['C2'])*Vrf
        #print("par['Vp_amp'] = ",par['Vp_amp'])
        par['Vp_avg'] = (par['Vp_amp']+np.abs(par['Vf'])-par['Te']/2
                         *np.log(2*np.pi*par['Vp_amp']/par['Te']))
        par['Vp_symm'] = Vrf/2+np.abs(par['Vf'])-par['Te']/2*np.log(np.pi*Vrf/par['Te'])
        par['Ubias'] = 2*par['Vp_avg']-Vrf
        #print("par['2*Vp_avg'] = ",2*par['Vp_avg'])
        #print("par['Ubias'] = ",par['Ubias'])
        par['V1_avg'] = par['Vp_avg']-par['Ubias']
        #print("par['V1_avg'] = ",par['V1_avg'])
        par['V2_avg'] = par['Vp_avg']
        par['Sohm'] = (1.73*par['m']*par['hl']/(2*par['e'])*par['e0']
                       *np.pow(par['omega'],2)*par['vm']*np.sqrt(par['Te'])
                       *np.sqrt(np.abs(par['V1_avg']))*par['d'])
        for i in 1,2: #Sstoc[1,2] = f(V[1,2]_avg)
            par[f'Sstoc{i}'] = (0.45*np.sqrt(par['m']/par['e'])*par['e0']
                                *np.pow(par['omega'],2)*np.sqrt(par['Te'])
                                *np.abs(par[f'V{i}_avg']))
        par['ns1'] = ((par['Sohm']+par['Sstoc1']+par['Sstoc2'])
                      /(2*par['e']*par['ub']
                               *(par['xi_c']+2*par['Te']+par['Te']
                                      *np.sqrt(np.log(gas_params[par['gas']]['M']
                                                      /(2*np.pi*par['m'])))+par['Te']/2)))  
        par['n0'] = par['ns1']/par['hl']
        par['V'] = 0.83*par['V1_avg']
        par['Ji_symm'] = par['e']*par['ns1']*par['ub']
        par['Ji1'] = par['e']*par['ns1']*par['ub']
        #print("par['V'] = ",par['V'])
        par['sm1'] = np.sqrt(0.82*par['e0']*np.pow(par['V'],3/2)/par['Ji1']
                             *np.sqrt(2*par['e']/gas_params[par['gas']]['M']))
        par['J1'] = 1.23*par['omega']*par['e0']/par['sm1']*par['V1_avg']
        #Cycle
        par['sm2'] = np.pow(par['S2']/par['S1'],1-1.5)*par['sm1']
        par['d'] = par['l']-par['sm1']-par['sm2']
        par['hl'] = 0.86*np.power(3+par['d']/(2*par['lam_i']),(-1/2))
        #print("par['sm1'] = ",par['sm1'])
        #print("par['sm2'] = ",par['sm2'])
        par['ns2'] = np.sqrt(par['S1']/par['S2'])*par['ns1']
        par['Ji1'] = par['e']*par['ns1']*par['ub']
        par['Ji2'] = par['e']*par['ns2']*par['ub']
        par['C1'] = par['e0']*par['S1']/par['sm1']
        #print("par['C1'] = ",par['C1'])
        par['C2'] = par['e0']*par['S2']/par['sm2']
        #print("par['C2'] = ",par['C2'])
        par['Vp_amp'] = par['C1']/(par['C1']+par['C2'])*Vrf
        #print("par['Vp_amp'] = ",par['Vp_amp'])
        par['Vp_avg'] = (par['Vp_amp']+np.abs(par['Vf'])-par['Te']/2
                         *np.log(2*np.pi*par['Vp_amp']/par['Te']))
        par['Vp_symm'] = Vrf/2+np.abs(par['Vf'])-par['Te']/2*np.log(np.pi*Vrf/par['Te'])
        par['Ubias'] = 2*par['Vp_avg']-Vrf
        par['V1_avg'] = par['Vp_avg']-par['Ubias']
        par['V2_avg'] = par['Vp_avg']
        par['Sohm'] = (1.73*par['m']*par['hl']/(2*par['e'])*par['e0']
                       *np.pow(par['omega'],2)*par['vm']*np.sqrt(par['Te'])
                       *np.sqrt(np.abs(par['V1_avg']))*par['d']) 
        for i in 1,2: #Sstoc[1,2] = f(V[1,2]_avg)
            par[f'Sstoc{i}'] = (0.45*np.sqrt(par['m']/par['e'])*par['e0']
                                *np.pow(par['omega'],2)*np.sqrt(par['Te'])
                                *np.abs(par[f'V{i}_avg']))
        par['ns1'] = ((par['Sohm']+par['Sstoc1']+par['Sstoc2'])
                      /(2*par['e']*par['ub']
                               *(par['xi_c']+2*par['Te']+par['Te']
                                      *np.sqrt(np.log(gas_params[par['gas']]['M']
                                                      /(2*np.pi*par['m'])))+par['Te']/2)))  
        par['n0'] = par['ns1']/par['hl']
        par['V'] = 0.83*par['V1_avg']
        par['Ji1'] = par['e']*par['ns1']*par['ub']
        par['sm1'] = np.sqrt(0.82*par['e0']*np.pow(par['V'],3/2)/par['Ji1']*np.sqrt(2*par['e']/gas_params[par['gas']]['M']))
        par['J1'] = 1.23*par['omega']*par['e0']/par['sm1']*par['V1_avg']
        par['Sabs'] = (par['e']*par['ns1']*par['ub']
                       *(par['V1_avg']+par['xi_c']+2*par['Te']+par['Te']
                         *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
                         +0.5*par['Te'])*par['S1']
                    + (par['e']*par['ns2']*par['ub']
                       *(par['V2_avg']+par['xi_c']+2*par['Te']+par['Te']
                         *np.sqrt(np.log(gas_params[par['gas']]['M']/(2*np.pi*par['m'])))
                         +0.5*par['Te'])*par['S2']))
            
        if (verbose==True) :
            print(f'Vrf = {Vrf}')
            print(f'Sabs = {par["Sabs"]}')
            print(f'ns = {par["ns1"]}')
            print(f'Ji = {par["Ji1"]}')
        return par['Sabs'] - par['Pwr']
         
    def Vrf_(verbose=False):

        fun = lambda x : VrfEq(x, verbose=verbose)
        VrfSolve = optimize.root_scalar(fun, x0=50, x1=5000, xtol=1e-3, method='secant')
        #VrfSolve = optimize.root_scalar(fun, bracket=(10000,10), rtol=1e-3, method='brentq')
        if (verbose==True): print(f'VrfSolve = {VrfSolve}')
        return VrfSolve.root
    
    if mode=='Pwr':
        Vrf = Vrf_(verbose=verbose) 
        par['Vrf'] = Vrf
    elif mode=='Vrf':
        VrfEq(par['Vrf'])
        par['Pwr'] = par['Sabs']

    return par

def dECalc(par, gas_params, verbose=False):
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
    par['Ns2'] = par['Ji2']/par['e']
    j_ = ['1','2','1_avg','2_avg']
    if verbose == True:
        for j in [0,1]:
            print(f"Integ_floor = {par[f'V{j_[j+2]}']-(par[f'dEi{j_[j]} [eV]']/2)}")
            print(f"First_power = {par[f'V{j_[j+2]}']-(par[f'dEi{j_[j]} [eV]']/2)-par[f'V{j_[j+2]}']}")
            print(f"Second_power = {1-4*np.power(par[f'V{j_[j+2]}']-(par[f'dEi{j_[j]} [eV]']/2)-par[f'V{j_[j+2]}'],2)}")
    return par

def FiENorm(par, gas_params, verbose=False):
    j_ = ['1','2','1_avg','2_avg']
    el = ['Mo', 'B', 'Si', 'MoO3']
    for j in [0,1]:
        if verbose == True:
            print(f"dEi{j_[j]} [eV] = {par[f'dEi{j_[j]} [eV]']}ddd")
            print(f"V{j_[j+2]} = {par[f'V{j_[j+2]}']}ddd")
        FiE = lambda x: (2*par[f'Ns{j_[j]}']/(par['omega']*par[f'dEi{j_[j]} [eV]'])
                            *np.power(1-4*np.power(x-par[f'V{j_[j+2]}'], 2)
                                             /np.power(par[f'dEi{j_[j]} [eV]'], 2), -0.5))

        par[f'NN{j_[j]}'] = integrate.quad(FiE, par[f'V{j_[j+2]}']-(par[f'dEi{j_[j]} [eV]']/2),
                                           par[f'V{j_[j+2]}']+(par[f'dEi{j_[j]} [eV]']/2))[0]
        
        for ind,i in enumerate(el):
            cont = [[],[]]
            with open(f'lib/sputtering_yeild/{i};{par["gas"]}.csv', newline='') as f:
                reader = csv.reader(f)
                for k,row in enumerate(reader):
                    if k<9: ...
                    else: 
                        cont[0].append(np.log10(float(row[0])))
                        cont[1].append(np.log10(float(row[1])))
                f.close()
            """⚠Main problem №1⚠"""
            log_interp = Akima1DInterpolator(cont[0], cont[1], extrapolate=True) #interp1d(cont[0], cont[1], fill_value=0)
            interp = lambda y: np.power(10.0, log_interp(np.log10(y)))                                                                          
            """gam = lambda x: (gas_params[par['gas']][f'y0{i}']+2*gas_params[par['gas']][f'A{i}']/np.pi
                                *(gas_params[par['gas']][f'w{i}']/
                                  (4*np.power(x/par['e']-gas_params[par['gas']][f'xc{i}'], 2)
                                   +np.power(gas_params[par['gas']][f'w{i}'], 2))))"""
            """⚠Main problem №2⚠"""       
            gam1_eff = lambda x: (interp(x)*FiE(x)/par[f'NN{j_[j]}'])
            if par[f'V{j_[j+2]}']>0 and par[f'V{j_[j+2]}']-(par[f'dEi{j_[j]} [eV]']/2)>0: 
                par[f'G{j_[j]}_{i}_eff'] = (integrate.quad(gam1_eff, par[f'V{j_[j+2]}']-(par[f'dEi{j_[j]} [eV]']/2),
                                                           par[f'V{j_[j+2]}']+(par[f'dEi{j_[j]} [eV]']/2))[0])
            else: par[f'G{j_[j]}_{i}_eff'] = 0
            """Проверка адекватности показаний. Всё более 1e6 и менее 10^-6 приравнивается к нулю,  т.к. не явялется адекватным"""
            if par[f'G{j_[j]}_{i}_eff']>1e6 or par[f'G{j_[j]}_{i}_eff']<1e-10: par[f'G{j_[j]}_{i}_eff'] = 0
            if (verbose==True):
                print(f"{j_[j]}_{i} = {par[f'G{j_[j]}_{i}_eff']}")
                print(f"V{j_[j+2]} = {par[f'V{j_[j+2]}']}")
                print(f"{'='*30}")
        
        """Designated MoO3 block (reverted)
        try:
            MoO3 = {'xi': 40, 'y0': 0.2262, 'xc': 40, 'w': 318.32173, 'A': -113.01079}
            gamMoO3 = lambda y: (MoO3['y0']+(2*MoO3['A']/np.pi)*(MoO3['w']/(4*np.power(y-MoO3['xc'],2)+np.power(MoO3['w'],2))))
            gamMoO3_eff = lambda x: (gamMoO3(x)*FiE(x)*np.heaviside(x-MoO3['xi'],0)/par[f'NN{j_[j]}'])
            if par[f'V{j_[j+2]}']>0 and par[f'V{j_[j+2]}']-(par[f'dEi{j_[j]} [eV]']/2)>0: 
                par[f'G{j_[j]}_MoO3_eff'] = (integrate.quad(gamMoO3_eff, par[f'V{j_[j+2]}']-(par[f'dEi{j_[j]} [eV]']/2),
                                                            par[f'V{j_[j+2]}']+(par[f'dEi{j_[j]} [eV]']/2))[0])
            else: par[f'G{j_[j]}_MoO3_eff'] = 0
            '''Проверка адекватности показаний. Всё более 1e6 и менее 10^-6 приравнивается к нулю,  т.к. не явялется адекватным'''
            if par[f'G{j_[j]}_MoO3_eff']>1e6 or par[f'G{j_[j]}_MoO3_eff']<1e-6: par[f'G{j_[j]}_MoO3_eff'] = 0
        except: print('⚠⚠⚠⚠⚠MoO3 failed⚠⚠⚠⚠⚠')
        ======================"""
        
    for j in [0,1]:        
        par[f'd{j_[j]}_MoB'] = par[f'G{j_[j]}_Mo_eff']/(par[f'G{j_[j]}_B_eff']+1e-10)
        par[f'd{j_[j]}_MoO3'] = par[f'G{j_[j]}_Mo_eff']/(par[f'G{j_[j]}_MoO3_eff']+1e-10)
        
        for i in el:
            par[f'd{j_[j]}M_{i}'] = par[f'Ji{j_[j]}']/(par['e']*par['N_A']*1e3)*par[f'M_{i}']*par[f'G{j_[j]}_{i}_eff']
            par[f'd{j_[j]}H_{i}'] = par[f'd{j_[j]}M_{i}']/par[f'ro_{i}']*6e10
            if i!='Mo':
                try:
                    par[f's{j_[j]}_{i}'] = (par[f'G{j_[j]}_{i}_eff']-par[f'G{j_[j]}_Mo_eff'])/(par[f'G{j_[j]}_{i}_eff']+par[f'G{j_[j]}_Mo_eff'])
                except: 
                    par[f's{j_[j]}_{i}'] = 'nan'
    
    return par
    
            


def main(vrname=False, prname=False, gcname=False, verbose=False):
        err = False
        try: 
            try: par, gas_params = ParLoader(prname=prname, vrname=vrname, gcname=gcname)
            except: 
                err = 1#"Unable to load initial parameters"
                if verbose == True: 
                    print(err)
                sys.exit(1)
            try: par.update(TeCalc(par, gas_params, verbose=verbose))
            except: 
                err = 2#"Electron temperature and/or rate coefficients can't be calculated"
                if verbose == True: 
                    print(err)
                sys.exit(2)
            if 'Pwr' in par:
                try:
                    par_ = VrfCalc(par, gas_params, mode='Pwr', verbose=verbose)
                    par.update(par_)
                except: 
                    err = 31#"Can't calculate Vrf from given Pwr"
                    if verbose == True: 
                        print(err)
                    sys.exit(31)
            elif 'Vrf' in par:
                try:
                    par_ = VrfCalc(par, gas_params, mode='Vrf', verbose=verbose)
                    par.update(par_)
                except:
                    err = 32#"Can't calculate Sabs from given Vrf"
                    if verbose == True: 
                        print(err)
                    sys.exit(32)           
            try:
                par.update(dECalc(par, gas_params, verbose=verbose))
            except:
                err = 5#"Can't calculate dEi"
                if verbose == True: 
                    print(err)
                sys.exit(5)
            try:
                par.update(FiENorm(par, gas_params, verbose=verbose))
            except:
                err = 6#"Unable to find etching params"
                if verbose == True: 
                    print(err)
                sys.exit(6)
               
        finally: 
            par['f']=par['f']/1e6
            return par
        
#par = main(verbose=False)
#print(par, len(par))
#print('Execution time: %s seconds' % (time.time() - start_time))
#sys.exit(0)