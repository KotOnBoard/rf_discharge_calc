import pandas as pd
import argparse
import matplotlib.gridspec as grd
import matplotlib.pyplot as plt
import json5
import os
import string
from mpl_toolkits import axisartist
from mpl_toolkits.axes_grid1 import host_subplot

def ArgParcer():
    par = argparse.ArgumentParser(description="A script that accepts keyword arguments.")
    par.add_argument('--tablename', '--tablename', required=False)
    par.add_argument('--filename', '--filename', required=False)
    par.add_argument('--construct', '--construct')
    par.add_argument('--dim', '--dim', help='Dimension of figures', default='2')
    par.add_argument('--shared', '--shared', help='Defines if any axes are shared', default=None)
    par.add_argument('--layout', '--layout', help='Figure layout', default='tight')
    par.add_argument('-show', '-show', action='store_true')
    par.add_argument('-filterby', '-filterby')
    par.add_argument('-trans', '-trans', action='store_true')
    par.add_argument('--out', '--out', help='Additional path for output.')
    par.add_argument('-log', '-log')
    args = par.parse_args()
    if not args.out:
        args.out = args.tablename
    return args

def RequestFormatter(args):
    if args.filename:
        with open(f"conf/{args.filename}", 'r') as file:
            ar = pd.Series(json5.load(file))
    elif args.construct:
        ar = pd.Series(json5.loads(args.construct))
    d = pd.read_excel(f'output/{args.tablename}.xlsx')
    for ind,arg in enumerate(ar):
        match arg:
            case str():
                if ',' in arg:
                    if args.filterby: print("⚠AmbiguousInputWarning: Both Terminal and File filterby arguments are given. In cases like that Terminal input takes priority. It's advised to use only one filterby origin.")
                    else:
                        filtr = arg.lstrip(string.ascii_lowercase).lstrip('([{<').rstrip('>)]}').split(',')
                        if len(filtr)%2: raise IndexError('Some parameters of "filterby" are missing. Number of parameters must be divisible by 2.')
                        else: 
                            a = list(range(0,len(filtr),2))
                            b = list(range(1,len(filtr),2))
                            for i,j in enumerate(b):
                                try:
                                    filtr[j] = float(filtr[j])
                                except: pass
                            for i,j in enumerate(a):
                                d = d.loc[d[filtr[j]]==filtr[b[i]]] #d.loc[d[ar[i][2]]==key[j],[ar[i][0],ar[i][1]]]
                        ar.pop(ind)
                #else: raise ValueError('Encountered invalid "filterby" formatting or other miscellanious string parameter.')
    if args.filterby:
        filtr = args.filterby.lstrip('([{').rstrip(')]}').split(',')
        if len(filtr)%2: raise IndexError('Some parameters of "filterby" are missing. Number of parameters must be divisible by 2.')
        else: 
            a = list(range(0,len(filtr),2))
            b = list(range(1,len(filtr),2))
            for i,j in enumerate(b):
                try:
                    filtr[j] = float(filtr[j])
                except: pass
            for i,j in enumerate(a):
                d = d.loc[d[filtr[j]]==filtr[b[i]]] #d.loc[d[ar[i][2]]==key[j],[ar[i][0],ar[i][1]]]
                
                    
        
    return d, ar
    
def GridSpecConstructor(num):   
    if num in range(1,3):
        nrows = 1
        ncols = num
    elif num==4:
        nrows = 2
        ncols = 2
    elif num in range(8,40,4):
        nrows = num//4
        ncols = 4
    elif num in range(6,90,3):
        nrows = num//3
        ncols = 3
    else:
        nrows = num//3+1
        ncols = 3           
    return nrows, ncols

def Plot(args, d, ar):
    #args.tablename = 'Test_1'
    #args.filename = 'filter2.json5'
    params = ['gas','p','l','R','f','Pwr','Assy']
    

    fig = plt.figure(layout=args.layout)
    ax = fig.add_subplot(visible=False)
    if args.trans:
        ncols, nrows = GridSpecConstructor(len(ar))
    else:
        nrows, ncols = GridSpecConstructor(len(ar))
    if type(ar[0])==type(list()):
        n = len(ar)
    else: 
        n = 1
        place = [[ar[i] for i in range(len(ar))]]
        ar = pd.Series(place)
    match args.dim:
        case '2':
            for i in range(n):
                try: params.pop(params.index(ar[i][0]))
                except: pass
                match args.shared:
                    case None: 
                        ax = [host_subplot(grd.SubplotSpec(grd.GridSpec(nrows,ncols), i), axes_class=axisartist.Axes)]
                        for j in range(1,len(ar[i])-1):
                            ax.append(ax[0].twinx())
                            ax[j].axis['right'] = ax[j].new_fixed_axis(loc="right", offset=(60*(j-1), 0))
                        for j in range(1,len(ar[i])):
                            p,=ax[j-1].plot(d[ar[i][0]], d[ar[i][j]], label=ar[i][j], marker='.')
                            ax[j-1].set_ylabel(ar[i][j])
                            if j==1: ax[j-1].axis['left'].label.set_color(p.get_color())
                            else: ax[j-1].axis['right'].label.set_color(p.get_color())
                    case 'x'|'y': 
                        ax = fig.add_subplot(grd.SubplotSpec(grd.GridSpec(nrows,ncols), i), sharex=ax)
                        for j in range(1,len(ar[i])):
                            ax.plot(d[ar[i][0]], d[ar[i][j]], label=ar[i][j], marker='.')
                    case 'y': ax = fig.add_subplot(grd.SubplotSpec(grd.GridSpec(nrows,ncols), i), sharey=ax)
                
                plt.grid(which='minor', linestyle='--')
                plt.grid(which='major', linestyle='--')
                plt.minorticks_on()
                if args.log: 
                    match args.log:
                        case 'y':
                            plt.yscale('log')
                        case 'x':
                            plt.xscale('log')
                        case 'xy'|'yx':
                            plt.yscale('log')
                            plt.xscale('log')
                plt.xlabel(ar[i][0])
                if args.shared:
                    plt.ylabel(str([ar[i][j] for j in range(1,len(ar[i]))]).lstrip('[').rstrip(']').replace("'",""))
        case '3':
            for i in range(len(ar)):
                try: params.pop(params.index(ar[i][0]))
                except: pass
                key = pd.unique(d[ar[i][2]])
                match args.shared:
                    case None: ax = fig.add_subplot(grd.SubplotSpec(grd.GridSpec(nrows,ncols), i))
                    case 'x': ax = fig.add_subplot(grd.SubplotSpec(grd.GridSpec(nrows,ncols), i), sharex=ax)
                    case 'y': ax = fig.add_subplot(grd.SubplotSpec(grd.GridSpec(nrows,ncols), i), sharey=ax)        
                for j in range(0,len(key)):
                    fil = d.loc[d[ar[i][2]]==key[j],[ar[i][0],ar[i][1]]]
                    plt.plot(fil[ar[i][0]], fil[ar[i][1]], label=f'{key[j]}', marker='.')
                plt.grid(which='minor', linestyle='--')
                plt.grid(which='major', linestyle='--')
                plt.minorticks_on()
                if args.log: 
                    match args.log:
                        case 'y':
                            plt.yscale('log')
                        case 'x':
                            plt.xscale('log')
                        case 'xy'|'yx':
                            plt.yscale('log')
                            plt.xscale('log')
                plt.xlabel(ar[i][0])
                plt.ylabel(ar[i][1])
    plt.legend()
    title = ''
    for i in params:
        title+=f"{i} = "
        title+=f"{str(pd.unique(d[i])).replace(' ',', ').strip('[({})]')}, "
    fig.suptitle(title.rstrip(", "))
    
    if not os.path.isfile(f'output/{args.out}.png'):
        plt.savefig(f'output/{args.out}.png', dpi=300)
    else:
        filenum = 1
        while (os.path.isfile(f'output/{args.out}_{filenum}.png')==True):
            filenum+=1
        else: plt.savefig(f'output/{args.out}_{filenum}.png', dpi=300)
    if args.show:
        plt.show()

def main(tablename=False):
    match tablename:
        case pd.DataFrame():
            ...
        case False:
            args = ArgParcer()
            d, ar = RequestFormatter(args)
            Plot(args, d, ar)
    
main()
