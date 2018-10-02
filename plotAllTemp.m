

clear all 
% close all


set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


FSZ = 24;
LW = 3;
MSZ = 10;
TickLength = [0.025 0.05];
MarkerSize = 10;
MarkerLineWidth = 1;
LineWidth = 2;

widthPos = 0.265;
heightPos = 0.275;

plotmovey = -0.0;





%dlmread('/Users/renep/Dropbox/Research/RANS_pipe/DNS_Hassan/Wall_Temp/...)
dataA = dlmread('DNS_Hassan/Wall_Temp/Abulk');
dataB = dlmread('DNS_Hassan/Wall_Temp/Bbulk');
dataC = dlmread('DNS_Hassan/Wall_Temp/Cbulk');
dataD = dlmread('DNS_Hassan/Wall_Temp/Dbulk');
dataE = dlmread('DNS_Hassan/Wall_Temp/Ebulk');   % downward
dataF = dlmread('DNS_Hassan/Wall_Temp/Fbulk');   % forced high Q
dataG = dlmread('DNS_Hassan/Wall_Temp/Gbulk');   % mixed high Q, case G
dataH = dlmread('DNS_Hassan/Wall_Temp/Hbulk');   % C 60 long
dataJ = dlmread('DNS_Hassan/Wall_Temp/Jbulk');   % A 60 long


skip = 100;
cas = 'cAL';
widthPos = 0.38;
heightPos = 0.38;

figure(1); 


% ------------------------ SA -------------------------
subplot(2,2,1); hold off
h1=plot(dataJ(1:skip:end,1), dataJ(1:skip:end,4),'ko','MarkerSize',MarkerSize,'LineWidth',MarkerLineWidth); hold on;
mod = 'SA'; 

vers = 'noMod'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
h2 = plot(dataRans(1:1:end,1)-0.15, T,'r--','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

vers = 'Aupoix'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
h3 = plot(dataRans(1:1:end/4,1)-0.15, T(1:1:end/4),'r:','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;
plot(dataRans(end/2:1:3*end/4,1)-0.15, T(end/2:1:3*end/4),'r:','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

vers = 'ModNew'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
h4 = plot(dataRans(end/4:1:end/2,1)-0.15, T(end/4:1:end/2),'b-','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;
plot(dataRans(3*end/4:1:end,1)-0.15, T(3*end/4:1:end),'b-','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

set(gca,'XTickLabel',[]);
ylabel('$T/T_0$');
set(gca, 'fontsize', 24, 'XMinorTick', 'on', 'YMinorTick', 'on')
set(gca, 'TickDir','out','TickLength',TickLength);

xlim([0 60]);
ylim([1 1.4]);
legend([h1,h2,h3,h4],'DNS', 'Conventional', 'Aupoix', 'Our mods')
title('Spalart \& Almaras', 'fontsize', 24)

pos = get(gca,'Position');
set(gca,'Position',[pos(1)    pos(2)-0.02    widthPos    heightPos]);



% ------------------------ MK -------------------------
subplot(2,2,2); hold off
h1=plot(dataJ(1:skip:end,1), dataJ(1:skip:end,4),'ko','MarkerSize',MarkerSize,'LineWidth',MarkerLineWidth); hold on;
mod = 'MK'; 

vers = 'noMod_yplus'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
plot(dataRans(1:1:end,1)-0.15, T,'r--','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

vers = 'Aupoix_yplus'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
plot(dataRans(1:1:end,1)-0.15, T,'r:','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

vers = 'ModNew'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
plot(dataRans(1:1:end,1)-0.15, T,'b-','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca, 'fontsize', 24, 'XMinorTick', 'on', 'YMinorTick', 'on')
set(gca, 'TickDir','out','TickLength',TickLength);

xlim([0 60]);
ylim([1 1.4]);
title('MK', 'fontsize', 24)

pos = get(gca,'Position');
set(gca,'Position',[pos(1)    pos(2)-0.02    widthPos    heightPos]);



% ------------------------ SST -------------------------
subplot(2,2,3); hold off
h1=plot(dataJ(1:skip:end,1), dataJ(1:skip:end,4),'ko','MarkerSize',MarkerSize,'LineWidth',MarkerLineWidth); hold on;
mod = 'OM'; 

vers = 'noMod'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
plot(dataRans(1:1:end,1)-0.15, T,'r--','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

vers = 'Aupoix'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
plot(dataRans(1:1:end,1)-0.15, T,'r:','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

vers = 'ModNew'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
plot(dataRans(1:1:end,1)-0.15, T,'b-','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

xlabel('$x/D$');
ylabel('$T/T_0$');
set(gca, 'fontsize', 24, 'XMinorTick', 'on', 'YMinorTick', 'on')
set(gca, 'TickDir','out','TickLength',TickLength);

xlim([0 60]);
ylim([1 1.4]);
title('SST', 'fontsize', 24)

pos = get(gca,'Position');
set(gca,'Position',[pos(1)    pos(2)    widthPos    heightPos]);



% ------------------------ V2F -------------------------
subplot(2,2,4); hold off
h1=plot(dataJ(1:skip:end,1), dataJ(1:skip:end,4),'ko','MarkerSize',MarkerSize,'LineWidth',MarkerLineWidth); hold on;
mod = 'VF'; 

vers = 'noMod'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
plot(dataRans(1:1:end,1)-0.15, T,'r--','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

vers = 'Aupoix'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
plot(dataRans(1:1:end,1)-0.15, T,'r:','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

vers = 'ModNew'; filename = sprintf('caseA/%s_%s',cas,mod);
dataRans = ReadRansX(sprintf('%s_%s/',filename,vers),4); T = 1.0*(dataRans(:,5)-1)+1;
plot(dataRans(1:1:end,1)-0.15, T,'b-','MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;

set(gca,'YTickLabel',[]);
xlabel('$x/D$','Interpreter','latex');
set(gca, 'fontsize', 24, 'XMinorTick', 'on', 'YMinorTick', 'on')
set(gca, 'TickDir','out','TickLength',TickLength);

xlim([0 60]);
ylim([1 1.4]);

title('V2F', 'fontsize', 24)

pos = get(gca,'Position');
set(gca,'Position',[pos(1)    pos(2)    widthPos    heightPos]);


set(gcf,'PaperUnits','points');
set(gcf,'PaperPosition',[0 0 1000 700]);
print(strcat('TempAllModelsCaseA.eps'),'-depsc')

