import f90nml
from subprocess import Popen, PIPE


inputfile = 'new_sample.nml'


patch_nml=dict()
patch_nml['input']=dict()


patch_nml['input']['eosmode']=0



f90nml.patch('input.nml', patch_nml, inputfile)
with Popen(['mpirun', '-np', '4', './forced_real_VF',inputfile ], stdout=PIPE, bufsize=1 ) as p:
    for line in p.stdout:
        print(line.decode().strip())
