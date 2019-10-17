import f90nml
from subprocess import Popen, PIPE


template_f = 'input.nml'
input_f    = 'input_new.nml'


adjust=dict()
adjust['input']=dict()
adjust['input']['eosmode']=0


f90nml.patch(template_f, adjust, input_f)
with Popen(['mpirun', '-np', '4', './forced_real_VF',input_f ], stdout=PIPE, bufsize=1 ) as p:
    for line in p.stdout:
        print(line.decode().strip())
