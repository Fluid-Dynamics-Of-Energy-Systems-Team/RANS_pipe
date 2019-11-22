#python script to generate the periodic inflow profiles

import f90nml
from subprocess import Popen, PIPE
import numpy as np

template_f = 'input.nml'
input_f    = 'input_new.nml'


adjust=dict()
adjust['input']=dict()

turbmodels = np.arange(1)
imaxes = [96,192,192,192,192]
kmaxes = [100,200,400,800,1600]
for index,(imax,kmax) in enumerate(zip(imaxes,kmaxes)):
   adjust['input']['turbmod']=0
   adjust['input']['CFL']=max(0.7,10/kmax)
   adjust['input']['K_start_heat']=int(kmax/10)
   adjust['input']['systemSolve']=4
   adjust['input']['periodic']=2
   adjust['input']['EOSmode']=0
   adjust['input']['Re']=1e6
   adjust['input']['nstep']=int(1e8)
   adjust['input']['LoD']=200
   adjust['input']['imax']=int(imax)
   adjust['input']['kelem']=int(kmax)
   adjust['input']['isothermalBC']=0
   adjust['input']['Qwall']=0
   adjust['input']['output_fname']   ="laminarbl_"   +"_".join([str(imax),str(kmax)])+".csv"
   adjust['input']['output_fname_bl']="laminarbl_pp_"+"_".join([str(imax),str(kmax)])+".csv"
   if (index > 0):
      adjust['input']['select_init']=2
      imax_old = imaxes[index-1]
      kmax_old = kmaxes[index-1]
      adjust['input']['kelem_old']=kmax_old
      adjust['input']['imax_old']=imax_old
      adjust['input']['read_fname']="laminarbl_"+"_".join([str(imax_old),str(kmax_old)])+".csv"
   else:
      adjust['input']['select_init']=0
   adjust['input']['dtmax']= max(0.001,1./kmax)
   f90nml.patch(template_f, adjust, input_f)
   with Popen(['mpirun', '-np', '4', './run',input_f ], stdout=PIPE, bufsize=1 ) as p:
      for line in p.stdout:
         print(line.decode().strip())
