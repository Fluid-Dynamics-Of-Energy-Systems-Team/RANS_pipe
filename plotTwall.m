clear all 
%close all


set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

ncore = 1;
FSZ = 24;
LW = 3;
MSZ = 10;
TickLength = [0.025 0.05];
MarkerSize = 8;
LineWidth = 2;

widthPos = 0.265;
heightPos = 0.275;

plotmovey = -0.0;

titleSpec  = '%s solving %s';                              % title of the plot
nameSpec   = '%s';                              % name of the output file
fileSpec   = '%s_%s_%s_%s';     % e.g. cCL_MK_noMod_DWX
fileSpec2   = '%s_%s_%s_%s_%s'; % e.g. cCL_MK_noMod_yplus_DWX

%Reading DNS data
% hassan
dataA = dlmread('DNS_Hassan/Wall_Temp/Abulk');
dataB = dlmread('DNS_Hassan/Wall_Temp/Bbulk');
dataC = dlmread('DNS_Hassan/Wall_Temp/Cbulk');
dataD = dlmread('DNS_Hassan/Wall_Temp/Dbulk');
dataE = dlmread('DNS_Hassan/Wall_Temp/Ebulk');   % downward
dataF = dlmread('DNS_Hassan/Wall_Temp/Fbulk');   % forced high Q
dataG = dlmread('DNS_Hassan/Wall_Temp/Gbulk');   % mixed high Q, case G
dataH = dlmread('DNS_Hassan/Wall_Temp/Hbulk');   % C 60 long
dataJ = dlmread('DNS_Hassan/Wall_Temp/Jbulk');   % A 60 long
% bae
datA = dlmread('DNS_Bae/A.txt');   
datB = dlmread('DNS_Bae/B.txt');   
datC = dlmread('DNS_Bae/C.txt');  
datD = dlmread('DNS_Bae/D.txt'); 
datE = dlmread('DNS_Bae/E.txt'); 





skipH = 100;
skipB = 2;
 

turbMod             = {'MK','VF'};
turbMod_name        = {'MK','V2F'};
diffMod_name        = {'noMod','modNew','Aupoix'};
tempturbMod_name    = {'const','Prt','DWX','DWXZh'};
nturb =  2;
cases_vect    = {'cAL','cBL', 'cCL','cDL', 'cEL'};
ncases= 1;
y_plus =0;  %1: plot noMod and Aupox with yplus, 0: with ystar
figure(1); hold off

mod =  turbMod{nturb};
cas = cases_vect{ncases};
data = 'H'; %H or B, hassan and bae, respectively
if(data=='H') %hassan data
    switch cas
        case {'cA','cAL'}
            h1=plot(dataJ(1:skipH:end,1), dataJ(1:skipH:end,4),'ko','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
            tmax = max(dataJ(1:end,4));
        case {'cB','cBL'} 
            h1=plot(dataB(1:skipH:end,1), dataB(1:skipH:end,4),'ko','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
            tmax = max(dataB(1:end,4));
        case {'cC','cCL'}
            h1=plot(dataH(1:skipH:end,1), dataH(1:skipH:end,4),'ko','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
            tmax = max(dataH(1:end,4));
        case {'cD','cDL'}
            h1=plot(dataD(1:skipH:end,1), dataD(1:skipH:end,4),'ko','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
            tmax = max(dataD(1:end,4));
        case {'cE','cEL'} 
            h1=plot(dataE(1:skipH:end,1), dataE(1:skipH:end,4),'ko','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
            tmax = max(dataE(1:end,4));
    end
else %bae data
    switch cas
        case {'cA','cAL'}
            h11=plot(datA(1:skipB:end,1), datA(1:skipB:end,2),'k^','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
            tmax = max(datA(1:end,2));
        case {'cB','cBL'} 
            h11=plot(datB(1:skipB:end,1), datB(1:skipB:end,2),'k^','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
            tmax = max(datB(1:end,2));
        case {'cC','cCL'}
            h11=plot(datC(1:skipB:end,1), datC(1:skipB:end,2),'k^','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
            tmax = max(datC(1:end,2));
        case {'cD','cDL'}
            h11=plot(datD(1:skipB:end,1), datD(1:skipB:end,2),'k^','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
            tmax = max(datD(1:end,2));
        case {'cE','cEL'} 
            h11=plot(datE(1:end,1), datE(1:end,2),'k^','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
            tmax = max(datE(1:end,2));
    end
end

%%mod New
filename2 = sprintf('%s/',turbMod{nturb}); %sprintf('Results/Cluster/%s_%s_modNew_DWX/',cas,mod);
dataRans = ReadRansX(filename2,4);
T = 1.0*(dataRans(1:1:end,5)-1)+1;
h2=plot(dataRans(1:1:end,1)-0.15, T,'bx','MarkerSize',MarkerSize,'LineWidth',LineWidth);


filename2 = sprintf('Results/%s_%s_noMod_const/',cas,mod);
dataRans = ReadRansX(filename2,4);
T = 1.0*(dataRans(1:1:end,5)-1)+1;
h3=plot(dataRans(1:1:end,1)-0.15, T,'r-','MarkerSize',MarkerSize,'LineWidth',LineWidth);

filename2 = sprintf('Results/%s_%s_modNew_const/',cas,mod);
dataRans = ReadRansX(filename2,4);
T = 1.0*(dataRans(1:1:end,5)-1)+1;
h4=plot(dataRans(1:1:end,1)-0.15, T,'g-.','MarkerSize',MarkerSize,'LineWidth',LineWidth);

% filename2 = sprintf('Results/%s_%s_modNew_const_Sochi/',cas,mod);
% dataRans = ReadRansX(filename2,4);
% T = 1.0*(dataRans(1:1:end,5)-1)+1;
% h5=plot(dataRans(1:1:end,1)-0.15, T,'m.','MarkerSize',MarkerSize,'LineWidth',LineWidth);


xlabel('Streamwise distance, $x/D$','Interpreter','latex');
ylabel('Wall temperature, $T/T_0$','Interpreter','latex');


title(sprintf(titleSpec,turbMod_name{nturb},cases_vect{ncases})); %,'FontWeight','bold');
titletext = get(gca,'title');
set(gca,'fontsize', 24)

set(gca, 'XMinorTick', 'on')
% set(gca, 'XTick', [0:5:30])

set(gca, 'YMinorTick', 'on')
% set(gca, 'YTick', [0:5:30])

ymax = 1.0*tmax + (tmax-1);
if(ymax<1.2)
    ymax=1.2;
end
xlim([0 60]);
ylim([1 ymax]);


if(exist('h5'))
    if(data=='H')
        legend([h2 h3 h4 h5 h1],'Standard TM + Sochi','Standard TM','Modified TM', 'Modified TM + Sochi','DNS','Location','southeast');
    else
        legend([h2 h3 h4 h5 h11],'Standard TM + Sochi','Standard TM','Modified TM', 'Modified TM + Sochi','DNS','Location','southeast');
    end
else
    if(data=='H')
        legend([h2 h3 h4 h1],'Simulating','Standard TM','Modified TM', 'DNS','Location','southeast');
    else
        legend([h2 h3 h4 h11],'Simulating','Standard TM','Modified TM', 'DNS','Location','southeast');
    end
end
