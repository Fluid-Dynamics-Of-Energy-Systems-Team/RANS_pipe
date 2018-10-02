

function [ data ] = ReadRansX(folder,ncpu)

    data=[];
    
    for i=0:ncpu-1
        fname = sprintf('%sprofX.%05d',folder,i);
        dread = dlmread(fname);
        if (i > 0) dread = dread(1:end,:); end
        data = [data;dread];
    end
    
end

