
% VARIABLES ="X", "Y", 
% "U", "W", "C", "T", (3-6) 
% "k", "eps", "v2",  "om", "nuSA", (7-11)
% "yplus", "RHO","Pe", mu, mut, ReTauStar (12-16)

function [ data ] = readTecplot(folder,ncpu, imax,kmax,nvar)

    data=[];
    if nargin<5
       nvar= 16; 
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

