
clc


plot_var = 4;   % VARIABLES ="X", "Y", (1,2) 
% "U", "W", "C", "T", (3-6) 
% "k", "eps", "v2", "om", "nuSA", (7-11)
% "yplus", "RHO","Pe", mu, mut (12-16)
rmax = 0.50005; rmin = 0.0;
cas  = 'cBL';      % cAL, cBL, cCL, cDL, CEL
mod  = 'OM';      % MK, VF, OM, SA, 
vers = 'noMod';   % noMod, modNew, Aupoix

filename  = sprintf('clusterAll/%s_%s',cas,mod);
filename2 = sprintf('%s_%s',filename,vers);

%data = readTecplot('/Users/renep/Dropbox/Research/RANS_pipe/caseA/MK/',4, 64,192);
%data = readTecplot('caseA/OM/',4, 96, 192); %64, 192);
%data = readTecplot(sprintf('%s/%s/',filename2,mod),4, 96, 192); %64, 192);
data = readTecplot('SA/',1, 96, 32,17); %64, 192);


figure(2)
contourf(data(:,:,1),data(:,:,2),data(:,:,plot_var),20)
colorbar

%caxis([36000. 38000])        
axis([0 20 rmin rmax])  % 0 x/L and 0 r
