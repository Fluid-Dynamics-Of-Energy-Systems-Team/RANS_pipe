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
    adjust['input']['LoD']=case['LOD']
    adjust['input']['turbdiffmod']=turbdiffmodel
    adjust['input']['modifDiffTerm']=modification
    return adjust

def write_file(casefolder, inputfile, outputfile, parameters):
    f90nml.patch(inputfile, parameters, os.path.join(casefolder,outputfile))
    return


def write_jobfile(name):
    f = open("job.template","r")
    string = f.read()
    f.close()
    string = string.replace("REPLACE_JOB", name)
    string = string.replace("REPLACE_JOB_LOG"+".log", name)
    f = open(os.path.join(name, "job"),"w")
    f.write(string)
    f.close()

if __name__== "__main__":
    submit = True
    periodic_template = "per_input.nml"
    developing_template = 'input.nml'

    #table
    table_loc="/home/azureuser/Documents/clean_repos/RANS_pipe/tables/co2_table.dat"
    nTab=2000
    fluid_name="co2"


    #mesh
    kelem = 384
    imax = 96
    K_start_heat=5

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
                        os.system("mkdir "+ name+'/pipe')
                        parameters = set_parameters(c,tm,mod,tdm)
                        parameters['input']['table_loc']=table_loc
                        parameters['input']['fluid_name']=fluid_name
                        parameters['input']['nTab']=nTab
                        parameters['input']['imax']=imax
                        parameters['input']['K_start_heat']=K_start_heat
                        parameters['input']['cfl']= cfl
                        parameters['input']['alphac']= alphac
                        parameters['input']['alphak']= alphak
                        parameters['input']['alphae']= alphae
                        parameters['input']['alphav2']= alphav2
                        parameters['input']['output_fname']=name+'.csv'
                        write_file(name,periodic_template,"periodic.nml",parameters)
                        parameters['input']['Qwall']=c['Qwall']
                        parameters['input']['kelem']=kelem
                        write_file(name,developing_template,"developing.nml" ,parameters)
                        write_jobfile(name)
                        #if submit:
                        #    p = Popen(['sbatch', 'job'], cwd=os.path.join(os.getcwd(),name))

