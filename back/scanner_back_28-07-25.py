import calc
import pandas as pd
from datetime import datetime
import time
import pyjson5
import sys
start_time = time.time()

def Loader(confname):
    with open(f"conf/{confname}", 'r') as file:
        par = pyjson5.load(file)
    return par

def Parcer(par):
    flush = list(par.keys())
    iter_par = {}
    for i in flush:
        if (type(par[f'{i}'])==list):
            iter_par[f'{i}']=par[f'{i}']
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
            __par, __out, __err = calc.calc(vrname=_par)
            par_ = pd.DataFrame(__par, index=[0])
            out_ = pd.DataFrame(__out, index=[0])
            output = pd.merge(par_, out_, left_index=True, right_index=True)
        else:
            __par, __out, __err = calc.calc(vrname=_par)
            par_ = pd.DataFrame(__par, index=[len(output)])
            out_ = pd.DataFrame(__out, index=[len(output)])
            _out = pd.merge(par_, out_, left_index=True, right_index=True)
            output = pd.concat([output, _out], axis=0)
    return output        
    
def main(confname='conf.json5', filename='parameters'):
    if len(sys.argv)>1 :
        if len(sys.argv)>2:
            confname, filename = sys.argv[1], sys.argv[2]
        else: 
            confname = sys.argv[1]
    par = Loader(confname)
    iter_par = Parcer(par)
    output = iterator(par, iter_par)
    output.to_excel(f'output/{filename}_{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}.xlsx')
    return output
    

main()
print('Execution time: %s seconds' % (time.time() - start_time))
