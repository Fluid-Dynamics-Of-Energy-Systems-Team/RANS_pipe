import f90nml
from subprocess import Popen, PIPE
import numpy as np
import os

#parameter cases hassan
cases = {
    "caseA60": {
        "Re" : 360,
        "Pr" : 3.19, 
        "Fr_1" : 0.0,
        "fprefix":'cA60',
        "Qwall":2.4,
        "LOD": 60*1.02
    },
    'caseB': {
        "Re" : 360,
        "Pr" : 3.19, 
        "Fr_1" : -9.96,
        "fprefix":'cC',
        "Qwall":2.4,
        "LOD": 30*1.02
    },
    "caseC60": {
        "Re" : 360,
        "Pr" : 3.19, 
        "Fr_1" : -79.67,
        "fprefix":'cC60',
        "Qwall":2.4,
        "LOD": 60*1.02
    },
    "caseD": {
        "Re" : 360,
        "Pr" : 3.19, 
        "Fr_1" : -268.89,
        "fprefix":'cD',
        "Qwall":2.4,
        "LOD": 30*1.02
    },
    "caseE": {
        "Re" : 360,
        "Pr" : 3.19, 
        "Fr_1" : 99.59,
        "fprefix":'cE',
        "Qwall":2.4,
        "LOD": 30*1.02
    }
}


turbmodels = {
 "SA": 1,
 "MK": 2,
 "SST": 3,
 "V2F": 4,
 "ABE": 5        
}

turbdiffmodels = {
    "CPrt": 0,
    "IF"  : 1,
    "Tang": 2,
    "KC"  : 3,
    "Kays": 4,
    "Bae" : 5,
    "DWX" : 6
}

modifications = {
    "std": 0,
    "otero": 1,
    "aupoix": 2
}






def case_name(case,turbmod,modification,turbdiffmodel):
    return "_".join([case['fprefix'],turbmod,turbdiffmodel,modification])


def set_parameters(case,turbmod,modification,turbdiffmodel):
    adjust=dict()
    adjust['input']=dict()
    adjust['input']['turbmod'] = turbmod
    adjust['input']['EOSmode'] = 2
    adjust['input']['Re']=case['Re']
    adjust['input']['Pr']=case['Pr']
    adjust['input']['Qwall']=case['Qwall']
    adjust['input']['LoD']=case['LOD']
    adjust['input']['turbdiffmod']=turbdiffmodel
    adjust['input']['modifDiffTerm']=modification
    return adjust

def write_file(casefolder, inputfile, outputfile, parameters):
    f90nml.patch(inputfile, parameters, os.path.join(casefolder,outputfile))
    return

if __name__== "__main__":
    periodic_template = "per_input.nml"
    developing_template = 'input.nml'


    #mesh
    kelem = 384
    imax = 96

    #relax parameters
    cfl     = 0.3
    alphac  = 0.1
    alphak  = 0.1
    alphae  = 0.1
    alphav2 = 0.1
    i=0
    for cname, c in cases.items():
        for tm_name, tm in turbmodels.items():
            for tdm_name,tdm  in turbdiffmodels.items():
                for mod_name, mod in modifications.items():
                    if not ((tm_name == "SA" or tm_name == "SST") and ( tdm_name == "DWX")):
                        name =case_name(c,tm_name, mod_name, tdm_name)
                        os.system("mkdir "+ name)
                        parameters = set_parameters(c,tm,mod,tdm)
                        write_file(name,periodic_template,"periodic.nml",parameters )
                        write_file(name,developing_template,"developing.nml" ,parameters)
                    
#     return 
# template_f = 'input.nml'
# input_f    = 'input_new.nml'


# adjust=dict()
# adjust['input']=dict()


# #turbulence model (1. 
# turbmodels = [1,2,3,4,5]
# turbdiffmodels = [0]
# modifdiffterms = [0,1]
# systemsolve =1,3)
# for tm in turbmodels:
#     for tdm in turbdiffmodels:
#         for mod in modifdiffterms:
#             adjust['input']['systemSolve']=1
#             adjust['input']['turbmod']=tm
#             adjust['input']['periodic']=1
#             adjust['input']['EOSmode']=2
#             adjust['input']['Re']=5300
#             adjust['input']['bulkmod']=1
#             adjust['input']['select_init']=0
#             adjust['input']['imax']=96
#             adjust['input']['kmax']=32
#             adjust['input']['isothermalBC']=0
#             adjust['input']['Qwall']=0
#             adjust['input']['modifDiffTerm'] = mod
#             adjust['input']['turbdiffmod']=tdm
#             f90nml.patch(template_f, adjust, input_f)
#             with Popen(['mpirun', '-np', '4', './run',input_f ], stdout=PIPE, bufsize=1 ) as p:
#                 for line in p.stdout:
#                     print(line.decode().strip())
