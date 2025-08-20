import calc
import pandas as pd
import time
import pyjson5
import os
import argparse
import numpy as np
import string
start_time = time.time()

def Loader(confname='dE_f.txt'):
    with open(f"conf/{confname}", 'r') as file:
        par = pyjson5.load(file)
    return par

def ArgParcer():
    parse = argparse.ArgumentParser(description="A script that accepts keyword arguments.")
    parse.add_argument('confname')
    parse.add_argument('filename')
    parse.add_argument('--plot', '--plot', action='store_true')
    parse.add_argument('--plotfilename', '--plotfilename', action='store_const')
    parse.add_argument('--filtername', '--filtername', action='store_const')
    args = parse.parse_args()
    #args.confname = 'conf.json5'
    #args.filename = 'test'
    
    
    return args

def ParameterParser(par):
    flush = list(par.keys())
    iter_par = {}
    for i in flush:
        match par[f'{i}']:
            case list():
                iter_par[f'{i}']=par[f'{i}']
                par.pop(f'{i}')
            case str():
                if ',' in par[f'{i}']:
                    unfold = par[f'{i}'].lstrip(string.ascii_lowercase).lstrip('([{<').rstrip('>)]}').split(',')
                    #print('===',unfold)
                    iter_par[f'{i}']=np.arange(float(unfold[0]),float(unfold[1]),float(unfold[2]))
                    par.pop(f'{i}')
    return iter_par

def iterator(par, iter_par):
    output = pd.DataFrame({})
    if (len(iter_par) > 0):
        output = evol(_par=par, iter_par=iter_par, output=output)
    else:
        calc.calc(vrname=par)   
    return output

def evol(_par: dict, iter_par:dict, output, recur=0):
    if (recur<len(iter_par.keys())):
        for i in list(iter_par[f'{list(iter_par.keys())[recur]}']):
            _par[f'{list(iter_par.keys())[recur]}'] = i
            output = evol(_par=_par, iter_par=iter_par, output=output, recur=recur+1)
    else:
        if (len(output)==0):
            __par, __err = calc.calc(vrname=_par)
            par_ = pd.DataFrame(__par, index=[0])
            output = par_
        else:
            __par, __err = calc.calc(vrname=_par)
            par_ = pd.DataFrame(__par, index=[len(output)])
            output = pd.concat([output, par_], axis=0)
    return output        
    
def main(confname='conf.txt', filename='Test', filtername=False):   
    args = ArgParcer()
    par = Loader(args.confname)
    iter_par = ParameterParser(par)
    output = iterator(par, iter_par)
    if args.filtername==None:
        if not os.path.isfile(f'output/{args.filename}.xlsx'):
            output.to_excel(f'output/{args.filename}.xlsx')
        else:
            filenum = 1
            while (os.path.isfile(f'output/{args.filename}_{filenum}.xlsx')==True):
                filenum+=1
            else: output.to_excel(f'output/{args.filename}_{filenum}.xlsx')
    else:
        with open(f"conf/{args.filtername}", 'r') as file:
            filt = pyjson5.load(file)
        if not os.path.isfile(f'output/{args.filename}.xlsx'):
            output.loc[:,filt].to_excel(f'output/{args.filename}.xlsx')
        else:
            filenum = 1
            while (os.path.isfile(f'output/{args.filename}_{filenum}.xlsx')==True):
                filenum+=1
            else: output.loc[:,filt].to_excel(f'output/{args.filename}_{filenum}.xlsx')
    return output
    

main()
print('Execution time: %s seconds' % (time.time() - start_time))
