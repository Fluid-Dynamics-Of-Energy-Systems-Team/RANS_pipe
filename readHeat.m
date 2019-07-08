

function [ data ] = readHeat(folder,ncpu, imax,kmax,nvar)

    data=[];
    if nargin<5
       nvar= 16; 
    end
    
    for i=0:ncpu-1
        fname = sprintf('%sheat.%05d',folder,i);
        
        if (i==0)
            dread = dlmread(fname, '', 2, 0); 
        else
            dread = dlmread(fname); 
        end
        
        data = [data;dread];
    end
    
    
    data = reshape(data,imax+2,kmax+2*ncpu,nvar);
    
end

