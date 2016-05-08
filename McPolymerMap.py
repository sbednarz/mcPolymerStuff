#
# Mapping of mcPolymer simulation results to a Python dictionary
#


import numpy as np
import re
import glob


def processMwd(case):
    mwd = {}
    files = glob.glob(case+'/*.HlogM.dat')
    for f in files:
        s = re.search('\/(.*?)\.', f).group(1)
        t = re.search('(\d\..*\d)', f).group(0)
        time = float(t)
        if mwd.has_key(s):
            mwd[s][time] = np.genfromtxt(f)    
        else:
            mwd[s] = {}
    
    return mwd



def processCld(case):
    cld = {}
    files = glob.glob(case+'/*.cld')
    for f in files:
        
        with open(f) as tmpf:
            data = tmpf.readlines()
        
        data0 = data[0].split()
        data1 = data[1].split()
        #print data0
        #print data1
        s = re.search('\/(.*?)\.', f).group(1)
        t = re.search('(\d\..*\d)', f).group(0)
        time = float(t)
        if not cld.has_key(s):
            cld[s] = {}
        if not cld[s].has_key(time):
            cld[s][time] = {}
        cld[s][time]['cld'] = np.genfromtxt(f, skip_header=3)
        cld[s][time]['nchains'] = float(data0[1])
        cld[s][time]['nmonomers'] = float(data0[7])
        cld[s][time]['pn'] = float(data1[1])
        cld[s][time]['pw'] = float(data1[3])
        cld[s][time]['pd'] = float(data1[5])
    return cld



def processConc(case):
    conc = {}
    files = glob.glob(case+'/c*dat')
    for f in files:
        s = re.search('\/c(.*)\.', f).group(1)
        conc[s] = np.genfromtxt(f)
    return conc



def processModel(case):
    filename = glob.glob(case+'/sim*tcl')[0]
    with open(filename) as f:
        model = f.readlines()

    params = {'rates':{},
              'conc0':{},
              'monomers':{},
              'initiator':{}, 
              'species':{}, 
              'speciesmacro':{},
              'var':{}
             }
    
    for i in range(0,len(model)):
        line = model[i]
        k = line.split()
        if len(k)>0:
            
            if k[0]=='#var!':
                line = model[i+1]
                k = line.split()
                params['var']['param'] = k[1]
                params['var']['val'] = float(k[2])
                
            elif k[0]=='#desc!':
                line = model[i+1].rstrip()
                params['desc'] = line
                
            elif k[0]=='Monomer':
                params['monomers'][k[1]] = float(k[2])
                
            elif k[0]=='Initiator':
                params['initiator'][k[1]] = float(k[2])
                
            elif k[0]=='Species':
                params['species'][k[1]] = ''
                params['conc0'][ k[1] ] = 0.0
                
            elif k[0]=='SpeciesMacro':
                params['speciesmacro'][k[1]] = ''
                params['conc0'][ k[1] ] = 0.0
                
            elif k[0]=='Concentration':
                params['conc0'][ k[1] ] = float(k[2])
                
            elif k[0]=='RateConstant':
                params['rates'][ k[1] ] = float(k[2])
                
            elif k[0]=='Simulation':
                params['time'] = float(k[1])
            
            elif k[0]=='InitSimulation':
                params['nmc']= float(k[1])

    params['model'] = filename
    return params

#------------------

sim = {}

with open('desc.txt', 'r') as f:
    desc = f.read()

sim['desc'] = desc

cases = glob.glob('case*')
for case in cases:
    params = processModel(case)
    mwd = processMwd(case)
    cld = processCld(case)
    conc = processConc(case)
    
    key = "{}={}".format(p['var']['param'], p['var']['val'])
    sim[key] = {}
    sim[key]['params'] = params
    sim[key]['conc'] = conc
    sim[key]['cld'] = cld
    sim[key]['mwd'] = mwd
