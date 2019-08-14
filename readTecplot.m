
% 'VARIABLES ="X","Y","U","W",                 (1-4)
%             "C","T","k","eps",               (5-8)
%             "v2","omega","nuSA","yplus",     (9-12)
%             "RHO","Pe","mu","mut",           (13-16)
%             "lamcp","cp","alphat","kt",      (17-20)
%             "epst","Pk","Gk"                 (21)


function [ data ] = readTecplot(folder,ncpu, imax,kmax,nvar)

    data=[];
    if nargin<5
       nvar= 21; 
    end
    
    for i=0:ncpu-1
        fname = sprintf('%stecp.%05d',folder,i);
        
        if (i==0)
            dread = dlmread(fname, '', 2, 0); 
        else
            dread = dlmread(fname); 
        end
        
        data = [data;dread];
    end
    
    
    data = reshape(data,imax+2,kmax+2*ncpu,nvar);
    
end

