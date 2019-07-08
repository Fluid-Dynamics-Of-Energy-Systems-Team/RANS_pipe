% close all 

cd /Users/rpecnik/Dropbox/Research/PIPE

[data] = readTecplot("0/",1, 24,48,17);

figure(1); hold off 
contourf(data(:,:,1),data(:,:,2),data(:,:,4)); 
% contourf(data(:,:,4)); 
colorbar
% axis equal



figure(2); hold off 
y = 0.5-data(:,2,2);
plot(y, data(:,2,4), 'bx-', 'LineWidth', 2); hold on 
plot(y, 360*y.*(1-y), 'ro', 'LineWidth', 2); hold on 




