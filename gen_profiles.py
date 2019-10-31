import f90nml
from subprocess import Popen, PIPE
import numpy as np

template_f = 'input.nml'
input_f    = 'input_new.nml'


adjust=dict()
adjust['input']=dict()

turbmodels = np.arange(1,5)
systemsolve = np.arange(1,2)
for tm in turbmodels:
    for ss in systemsolve:
        adjust['input']['systemSolve']=ss
        adjust['input']['turbmod']=tm
        adjust['input']['periodic']=1
        adjust['input']['EOSmode']=0
        adjust['input']['select_init']=1
        adjust['input']['imax']=96
        adjust['input']['isothermalBC']=0
        adjust['input']['Qwall']=0
        f90nml.patch(template_f, adjust, input_f)
        with Popen(['mpirun', '-np', '4', './run',input_f ], stdout=PIPE, bufsize=1 ) as p:
            for line in p.stdout:
                print(line.decode().strip())
