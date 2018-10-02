

clear all 
% close all


set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


FSZ = 24;
LW = 3;
TickLength = [0.025 0.05];
MarkerSize = 10;
LineWidth = 4;

% widthPos = 0.41;
% heightPos = 0.38;
% plotmovey = -0.05;


colDat = 2;


skip = 3;
% tmp = dlmread('DNS_Bae/A.txt');  dataA = tmp(1:skip:end,:);
% tmp = dlmread('DNS_Bae/B.txt');  dataB = tmp(1:skip:end,:);
% tmp = dlmread('DNS_Bae/C.txt');  dataC = tmp(1:skip:end,:);
% tmp = dlmread('DNS_Bae/D.txt');  dataD = tmp(1:skip:end,:);
% tmp = dlmread('DNS_Bae/E.txt');  dataE = tmp(1:1:end,:);
% tmp = dlmread('DNS_Bae/F.txt');  dataA = tmp(1:skip:end,:);
% tmp = dlmread('DNS_Bae/G.txt');  dataA = tmp(1:skip:end,:);


skip = 1;
tmp = dlmread('DNS_Hassan/Wall_Temp/Jbulk');   dataA = [tmp(1:skip:end,1) tmp(1:skip:end,4)]; % A 60 long
tmp = dlmread('DNS_Hassan/Wall_Temp/Bbulk');   dataB = [tmp(1:skip:end,1) tmp(1:skip:end,4)];
tmp = dlmread('DNS_Hassan/Wall_Temp/Hbulk');   dataC = [tmp(1:skip:end,1) tmp(1:skip:end,4)];   % C 60 long
tmp = dlmread('DNS_Hassan/Wall_Temp/Dbulk');   dataD = [tmp(1:skip:end,1) tmp(1:skip:end,4)];
tmp = dlmread('DNS_Hassan/Wall_Temp/Ebulk');   dataE = [tmp(1:skip:end,1) tmp(1:skip:end,4)];   % downward
% tmp = dlmread('DNS_Hassan/Wall_Temp/Fbulk');   dataF = [tmp(1:skip:end,1) tmp(1:skip:end,4)];   % forced high Q
% tmp = dlmread('DNS_Hassan/Wall_Temp/Gbulk');   dataG = [tmp(1:skip:end,1) tmp(1:skip:end,4)];   % mixed high Q, case G

figure(1); 
% alimits = [0 50 1 1.2];

mod = 'MK';
mod2 = 'YP';
dir = '/Users/renep/Documents/Research/RANS_pipe/cluster/';

% subplot(3,2,1); 
hold off
% h1=plot(dataA(:,1), dataA(:,colDat),'k-','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
% h2=plot(dataB(:,1), dataB(:,colDat),'k--','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
% h3=plot(dataC(:,1), dataC(:,colDat),'k-.','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
% h4=plot(dataE(:,1), dataE(:,colDat),'k:','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;

h1=plot(dataA(:,1), dataA(:,colDat),'r-','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
h2=plot(dataB(:,1), dataB(:,colDat),'--','Color',[0.0 0.75 0],'MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
h3=plot(dataC(:,1), dataC(:,colDat),'b-.','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
% h4=plot(dataE(:,1), dataE(:,colDat),'k:','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;


axis([0 60 1 1.2]);
xlabel('Axial position $x/D$');
ylabel('Wall temperature $T/T_0$');
set(gca,'XMinorTick','on'); %,'XLabel',[],'XTickLabel',[])
set(gca,'YMinorTick','on'); %,'YLabel',[],'YTickLabel',[])
set(gca,'fontsize', FSZ,'ticklength',TickLength,'linewidth',1);
set(gca,'TickDir','out');

% set(gca,'GridAlpha',0.1);     

legend('Forced','Location','southeast')
legend('Forced', 'Low buoyancy, up','Location','southeast')
legend('Forced', 'Low buoyancy, up','High buoyancy, up','Location','southeast')
% legend('Forced', 'Low buoyancy, up','High buoyancy, up', 'Mod. buoyancy, down','Location','southeast')
% legend('Direct numerical simulation (Nemati et al. 2015)', 'Standard turbulence model', 'New model based on Pecnik and Patel, 2016', 'Location', 'northwest')
% set(gca,'TickDir','out');

set(gcf,'PaperUnits','points');
scaling = 1.1;
set(gcf,'PaperPosition',[0 0 750 450]*scaling);
print(strcat('DNS3.eps'),'-depsc')




