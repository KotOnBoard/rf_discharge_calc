import calc
import pandas as pd
import time
import pyjson5
import os
import argparse
import numpy as np
import string

def Loader(confname):
    with open(f"conf/{confname}", 'r') as file:
        _par = pyjson5.load(file)
        if len(_par)>1:
            par = _par[0]
            kwargs = list(_par[1:])
        else:
            par = _par[0]
            kwargs = None
    return par, kwargs

def ArgParcer():
    
    parse = argparse.ArgumentParser(description="A script that accepts keyword arguments.")
    parse.add_argument('--confname', '--confname', default=None)
    parse.add_argument('--filename', '--filename', default=None)
    #parse.add_argument('--plot', '--plot', action='store_true')
    #parse.add_argument('--plotfilename', '--plotfilename', action='store_const')
    #parse.add_argument('--filtername', '--filtername', action='store_const')
    parse.add_argument('-tutor', '-tutor', action='store_false', default=True)
    #parse.add_argument('args', nargs=argparse.REMAINDER)
    arg = parse.parse_args()
    #args.confname = 'conf.json5'
    #args.filename = 'test'
    
    return arg

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
    pairs = PairFinder(iter_par)
    return iter_par, pairs

def PairFinder(iter_par):
    keys = list(iter_par.keys())
    pairs={}
    for i, val in enumerate(iter_par.values()):
        length = len(val)
        if length not in pairs and length>1:
            pairs[length] = []
        if length>1:
            pairs[length].append(keys[i])
    for key in list(pairs):
        if len(pairs[key])<2: pairs.pop(key)    
    return pairs

def PairRevolver(pair, par, iter_par, progress):
    output = pd.DataFrame({})
    pair_par={}
    for key in pair:
        if key in iter_par: 
            pair_par[key] = iter_par[key]
            iter_par.pop(key)
    for index in range(len(pair_par[pair[0]])):
        for key in pair_par:
            par[key]=pair_par[key][index]
        out, err = RecurIterator(par, iter_par, progress)
        if (len(output)==0): 
            output = pd.DataFrame(out, index=[0])
        else: 
            _out = pd.DataFrame(out, index=[len(output)])
            output = pd.concat([output, _out], axis=0)
    return output, progress

def PairChooser(pairs, par):
    pair = input(f'Multiple ranges with equal length are detected: {pairs}\nPlease specify which pair to iterate over or press Enter to iterate over combinations.\n(e.g. p,Pwr)\n>>> ').split(',')
    match pair:
        case '':
            return None
        case list():
            print(pair, type(pair))
            return pair
        case _: raise ValueError("'pair' kwarg is provided but only one key is specified.")

def printProgressBar (progress, err, printEnd = "\r"):
    progress[0]+=1    
    match err:
        case 1: progress[1]+=1
        case 2: progress[2]+=1
        case 3: progress[3]+=1
        case 4: progress[4]+=1
        case 5: progress[5]+=1
        case 6: progress[6]+=1          
    print(f'\rIteration: {progress[0]}|Te: {progress[0]-np.sum(progress[1:3])}/{progress[0]} |Vrf: {progress[0]-np.sum(progress[1:4])}/{progress[0]} |Ubias: {progress[0]-np.sum(progress[1:5])}/{progress[0]} |dEi: {progress[0]-np.sum(progress[1:6])}/{progress[0]} |Selectivity: {progress[0]-np.sum(progress[1:7])}/{progress[0]}', end = printEnd)
    return progress
def RecurIterator(par, iter_par, progress):
    output = pd.DataFrame({})
    if (len(iter_par) > 0):
        output, progress = evol(_par=par, iter_par=iter_par, output=output, progress=progress)
    else:
        output, err = calc.calc(vrname=par) 
        progress = printProgressBar(progress, err)
    return output, progress

def evol(_par: dict, iter_par:dict, output, progress, recur=0):
    if (recur<len(iter_par.keys())):
        for i in list(iter_par[f'{list(iter_par.keys())[recur]}']):
            _par[f'{list(iter_par.keys())[recur]}'] = i
            output, progress = evol(_par=_par, iter_par=iter_par, output=output, progress=progress, recur=recur+1)
    else:
        if (len(output)==0):
            __par, err = calc.calc(vrname=_par)
            par_ = pd.DataFrame(__par, index=[0])
            output = par_
            progress = printProgressBar(progress, err)
        else:
            __par, err = calc.calc(vrname=_par)
            par_ = pd.DataFrame(__par, index=[len(output)])
            output = pd.concat([output, par_], axis=0)
            progress = printProgressBar(progress, err)
    return output, progress        

def Iterator(kwargs, iter_par, pairs, par, progress):
    match kwargs:
        case None:
            output, progress = RecurIterator(par, iter_par, progress)
        case list():
            if 'pair' in kwargs:
                if len(pairs)==1:
                    pair = pairs[list(pairs.keys())[0]]
                    if len(pair)==2:
                        output, progress = PairRevolver(pair, par, iter_par, progress)
                    elif len(pair)>2:
                        print(f'Warning: more than two parameters are paired: {pair}\nIteration will be completed over all of them.')
                        output, progress = PairRevolver(pair, par, iter_par, progress)
                elif len(pairs)>1:
                    pair = PairChooser(pairs, par)
                    if pair:
                        timer()
                        output, progress = PairRevolver(pair, par, iter_par, progress)
                    else:
                        timer()
                        output, progress = RecurIterator(par, iter_par, progress)
                else: raise IndexError("Can't iterate in pairs because parameter ranges are different lengths.")
    return output, progress

def timer():
    global start_time
    start_time = time.time()
    
def mainHelper(args, progress):
    if not(args.confname):
        args.confname = input('Enter coniguration file name with extension (e.g. confexample.json5)\n>>> ')
    par, kwargs = Loader(args.confname)
    if not(args.filename):
        filename = input(f'Specify output file name without extension and subdirectory (if needed) or press Enter for {args.confname.split('.')[0]}.xlsx\n(e.g. Example)\n>>> ')
        match filename:
            case '':
                args.filename = args.confname.split('.')[0]
            case _:
                args.filename = filename
    iter_par, pairs = ParameterParser(par)
    match kwargs:
        case list():
            if 'pair' in kwargs:
                if len(pairs)==1:
                    switch = input(f'Found one combination of equal-length lists: {pairs}.\n[Enter - confirm| N - for regular iteration| C - choose specific combination|Other - exit]\n>>> ')
                    match switch:
                        case '': 
                            pair = pairs[list(pairs.keys())[0]]
                            timer()
                            output, progress = PairRevolver(pair, par, iter_par, progress)
                        case 'N':
                            timer()
                            output, progress = RecurIterator(par, iter_par, progress)
                        case 'C':
                            pair = PairChooser(pairs, par)
                            if pair:
                                timer()
                                output, progress = PairRevolver(pair, par, iter_par, progress)
                            else:
                                timer()
                                output, progress = RecurIterator(par, iter_par, progress)
                        case _: raise KeyboardInterrupt('Program exit')
                elif len(pairs)>1:
                    pair = PairChooser(pairs, par)
                    if pair:
                        timer()
                        output, progress = PairRevolver(pair, par, iter_par, progress)
                    else:
                        timer()
                        output, progress = RecurIterator(par, iter_par, progress)
                else: raise IndexError("Can't iterate in pairs because parameter ranges are different lengths.")
        case _:
            match len(pairs):
                case 1:
                    switch = input(f'No kwargs given bit found one combination of equal-length lists: {pairs}.\n[Enter - confirm| N - for regular iteration| C - choose specific combination|Other - exit]\n>>> ')
                    match switch:
                        case '': 
                            pair = pairs[list(pairs.keys())[0]]
                            timer()
                            output, progress = PairRevolver(pair, par, iter_par, progress)
                        case 'N': 
                            timer()
                            output, progress = RecurIterator(par, iter_par, progress)
                        case 'C':
                            timer()
                            pair = PairChooser(pairs, par)
                            if pair:
                                timer()
                                output, progress = PairRevolver(pair, par, iter_par, progress)
                            else:
                                timer()
                                output, progress = RecurIterator(par, iter_par, progress)
                        case _: raise KeyboardInterrupt('Program exit')
                case 2:
                    pair = PairChooser(pairs, par)
                    if pair:
                        timer()
                        output, progress = PairRevolver(pair, par, iter_par, progress)
                    else:
                        timer()
                        output, progress = RecurIterator(par, iter_par, progress)
                case _: 
                    switch = input('No kwargs and no pairs found.\nPress Enter for combination iteration or type exit for exit.\n')
                    match switch:
                        case '': 
                            timer()
                            output, progress = RecurIterator(par, iter_par, progress)
                        case _: raise KeyboardInterrupt('Program exit')
    return output, progress, args
    
def main(confname='conf.txt', filename='Test', filtername=False):   
    args = ArgParcer()
    progress = [0,0,0,0,0,0,0]
    if args.tutor:
        output, progress, args = mainHelper(args, progress)
    else:
        #if args.confname:
        par, kwargs = Loader(args.confname)
        #elif args.tutor:
        #    args.confname = input('Define config name with valid parameters and subdirectory, if needed.\nExample: confexample.json5\n')
        #    par = Loader(args.confname)
        iter_par, pairs = ParameterParser(par)
        timer()
        output, progress = Iterator(kwargs, iter_par, pairs, par, progress)
    if True:#args.filtername==None:
        if '/' in args.filename:
            subdir = 'output/' + args.filename.rstrip(string.ascii_letters).rstrip('/')
            if not(os.path.isdir(subdir)):
                os.makedirs(subdir)
        if not os.path.isfile(f'output/{args.filename}.xlsx'):
            output.to_excel(f'output/{args.filename}.xlsx')
        else:
            filenum = 1
            while (os.path.isfile(f'output/{args.filename}_{filenum}.xlsx')==True):
                filenum+=1
            else: output.to_excel(f'output/{args.filename}_{filenum}.xlsx')
        '''else:
            with open(f"conf/{args.filtername}", 'r') as file:
                filt = pyjson5.load(file)
            if not os.path.isfile(f'output/{args.filename}.xlsx'):
                output.loc[:,filt].to_excel(f'output/{args.filename}.xlsx')
            else:
                filenum = 1
                while (os.path.isfile(f'output/{args.filename}_{filenum}.xlsx')==True):
                    filenum+=1
                else: output.loc[:,filt].to_excel(f'output/{args.filename}_{filenum}.xlsx')'''
        return output
    

main()
print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
