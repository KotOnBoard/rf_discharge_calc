import calc as calc
import pandas as pd
import time
import pyjson5
import os
import argparse
import numpy as np
import string
from tqdm import tqdm
from functools import reduce
from numpy import cumprod
from multiprocessing import Pool
import csv
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
    #arg.confname = 'PLA-55.json5'
    #args.filename = 'test'
    
    return arg

def ParameterParser(par): #Не очищать словарь от постоянных значений
    flush = list(par.keys())
    for i in flush:
        match par[f'{i}']:
            case str():
                if ',' in par[f'{i}']:
                    unfold = par[f'{i}'].lstrip(string.ascii_lowercase).lstrip('([{<').rstrip('>)]}').split(',')
                    par[f'{i}']=np.arange(float(unfold[0]),float(unfold[1]),float(unfold[2]))
    pairs = PairFinder(par)
    return par, pairs

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


def PairChooser(pairs, par):
    pair = input(f'Multiple ranges with equal length are detected: {pairs}\nPlease specify which pair to iterate over or press Enter to iterate over combinations.\n(e.g. p,Pwr)\n>>> ').split(',')
    match pair:
        case '':
            return None
        case list():
            print(pair, type(pair))
            return pair
        case _: raise ValueError("'pair' kwarg is provided but only one key is specified.")

def Iterator(args, par, pair=None):
    #pair_lenght = 
    _par = par
    pair_par={}
    if pair:
        for key in pair:
            if key in _par: 
                pair_par[key] = _par[key]
                _par.pop(key)
        ky = list(pair_par.keys())
        pl = len(list(pair_par.values())[0])
    else: pl = 1
    kx = list(_par.keys())
    ln_a = [len(i) for i in _par.values()] + [1, ]
    ln_a.reverse()
    x = cumprod(ln_a).tolist()
    x.reverse()
    ml = reduce(lambda x, y: x*y, ln_a)
    
    
    
    pool = Pool() #defaults to number of available CPU's

    prev_i=0
    arr = []
    for i in range(ml*pl):
        inter = []
        inter_d = {}
        if pair:
            for k in ky:
                inter_d[k] = pair_par[k][i//ml]
        _i = i%ml
        for j in range(len(x)-1):
            inter.append(_par[kx[j]][_i%x[j]//x[j+1]])
     
        for k,v in zip(kx, inter):
            inter_d[k] = v
        arr.append(inter_d)
        """if (i%10000 == 0) or (i == ml*pl-1):
            #print(arr)
            ret = pool.map(calc.main, arr) # тут как раз все заморачиваться должно
            #print(ret)
            dump_out(ret, args, prev_i+2)
            arr = []
            prev_i = i"""
    ret = pool.map(calc.main, tqdm(arr))
    #ret = [calc.main(i) for i in tqdm(arr)]
    print("Calculation complete! Exporting... (this migh take a while)")
    #dump_out(ret, args, prev_i+2)  
    dump_back(ret, args)
    pass

def dump_back(out, args):
    keys = out[0].keys()
    with open(args.path, 'w', newline='') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(out)

def dump_out(output, args, i):
    out = pd.DataFrame(output)
    if os.path.isfile(args.path): 
        with pd.ExcelWriter(args.path,mode='a',engine="openpyxl",if_sheet_exists="overlay") as writer:
            out.to_excel(writer, startrow=i, header=False, index=False)
    else: 
        with pd.ExcelWriter(args.path, engine='xlsxwriter') as writer:
            out.to_excel(writer, index=False)
    
    pass

      

def ForkLoader(args, kwargs, pairs, par):
    match kwargs:
        case None:
            Iterator(args, par)
        case list():
            if 'pair' in kwargs:
                if len(pairs)==1:
                    pair = pairs[list(pairs.keys())[0]]
                    if len(pair)==2:
                        start_time = time.time()
                        Iterator(args, par, pair)
                        print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                    elif len(pair)>2:
                        print(f'Warning: more than two parameters are paired: {pair}\nIteration will be completed over all of them.')
                        start_time = time.time()
                        Iterator(args, par, pair)
                        print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                elif len(pairs)>1:
                    pair = PairChooser(pairs, par)
                    if pair:
                        start_time = time.time()
                        Iterator(args, par, pair)
                        print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                    else:
                        start_time = time.time()
                        Iterator(args, par)
                        print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                else: raise IndexError("Can't iterate in pairs because parameter ranges are different lengths.")
    return
    
def FrameworkHelper(args):
    if not(args.confname):
        args.confname = input('Enter coniguration file name with extension (e.g. confexample.json5)\n>>> ')
    par, kwargs = Loader(args.confname)
    if not(args.filename):
        filename = input(f"Specify output file name without extension and subdirectory (if needed) or press Enter for {args.confname.split('.')[0]}.xlsx\n(e.g. Example)\n>>> ")
        match filename:
            case '':
                args.filename = args.confname.split('.')[0]
            case _:
                args.filename = filename
    args.path = path_setter(args)    
    par, pairs = ParameterParser(par)
    match kwargs:
        case list():
            if 'pair' in kwargs:
                if len(pairs)==1:
                    switch = input(f'Found one combination of equal-length lists: {pairs}.\n[Enter - confirm| N - for regular iteration| C - choose specific combination|Other - exit]\n>>> ')
                    match switch:
                        case '': 
                            pair = pairs[list(pairs.keys())[0]]
                            start_time = time.time()
                            Iterator(args, par, pair)
                            print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                        case 'N':
                            start_time = time.time()
                            Iterator(args, par)
                            print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                        case 'C':
                            pair = PairChooser(pairs, par)
                            if pair:
                                start_time = time.time()
                                Iterator(args, par, pair)
                                print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                            else:
                                start_time = time.time()
                                Iterator(args, par)
                                print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                        case _: raise KeyboardInterrupt('Program exit')
                elif len(pairs)>1:
                    pair = PairChooser(pairs, par)
                    if pair:
                        start_time = time.time()
                        Iterator(args, par, pair)
                        print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                    else:
                        start_time = time.time()
                        Iterator(args, par)
                        print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                else: raise IndexError("Can't iterate in pairs because parameter ranges are different lengths.")
        case _:
            match len(pairs):
                case 1:
                    switch = input(f'No kwargs given bit found one combination of equal-length lists: {pairs}.\n[Enter - confirm| N - for regular iteration| C - choose specific combination|Other - exit]\n>>> ')
                    match switch:
                        case '': 
                            pair = pairs[list(pairs.keys())[0]]
                            start_time = time.time()
                            Iterator(args, par, pair)
                            print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                        case 'N': 
                            start_time = time.time()
                            Iterator(args, par)
                            print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                        case 'C':
                            pair = PairChooser(pairs, par)
                            if pair:
                                start_time = time.time()
                                Iterator(args, par, pair)
                                print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                            else:
                                start_time = time.time()
                                Iterator(args, par)
                                print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                        case _: raise KeyboardInterrupt('Program exit')
                case 2:
                    pair = PairChooser(pairs, par)
                    if pair:
                        start_time = time.time()
                        Iterator(args, par, pair)
                        print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                    else:
                        start_time = time.time()
                        Iterator(args, par)
                        print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                case _: 
                    switch = input('No kwargs and no pairs found.\nPress Enter for combination iteration or type exit for exit.\n')
                    match switch:
                        case '': 
                            start_time = time.time()
                            Iterator(args, par)
                            print('\nCompleted successfully.\nExecution time: %s seconds' % (time.time() - start_time))
                        case _: raise KeyboardInterrupt('Program exit')
    return

def path_setter(args, ext='csv'):
    if '/' in args.filename:
        subdir = 'output/' + args.filename.rstrip(string.ascii_letters).rstrip('/')
        if not(os.path.isdir(subdir)):
            os.makedirs(subdir)
    if not os.path.isfile(f'output/{args.filename}.{ext}'):
        path = f'output/{args.filename}.{ext}'
    else:
        filenum = 1
        while (os.path.isfile(f'output/{args.filename}_{filenum}.{ext}')==True):
            filenum+=1
        else: path = f'output/{args.filename}_{filenum}.{ext}'
    return path
    
def Framework(confname='conf.txt', filename='Test', filtername=False):   
    args = ArgParcer()
    if args.tutor:
        FrameworkHelper(args)
    else:
        par, kwargs = Loader(args.confname)
        par, pairs = ParameterParser(par)
        args.path = path_setter(args)  
            
        ForkLoader(args, kwargs, pairs, par)
    return
    

if __name__ == '__main__':
    Framework()

