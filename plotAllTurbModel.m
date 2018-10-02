

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

titleSpec  = '%s';                              % title of the plot
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
ncases= 3;
y_plus =0;  %1: plot noMod and Aupox with yplus, 0: with ystar
figure(1); hold off
for j=ncases:ncases
    %figure(j);
    for i=nturb:nturb
        %subplot(2,nturb/2,i)
        mod = turbMod{i};
        cas = cases_vect{j};
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
        filename2 = sprintf('Results/Cluster/%s_%s_modNew_DWX/',cas,mod);
        dataRans = ReadRansX(filename2,4);
        T = 1.0*(dataRans(1:1:end,5)-1)+1;
        h2=plot(dataRans(1:1:end,1)-0.15, T,'bx','MarkerSize',MarkerSize,'LineWidth',LineWidth);
        if(ncases>1)
            filename2 = sprintf('Results/Cluster/%s_%s_modNew_DWXZh/',cas,mod); %filename2 = sprintf('Results/Cluster/%s_%s_modNew_NK/',cas,mod);
            dataRans = ReadRansX(filename2,4);
            T = 1.0*(dataRans(1:1:end,5)-1)+1;
            h3=plot(dataRans(1:1:end,1)-0.15, T,'b-.','MarkerSize',MarkerSize,'LineWidth',LineWidth);
        end
        
        filename2 = sprintf('Results/Cluster/%s_%s_modNew_Prt/',cas,mod);
        dataRans = ReadRansX(filename2,4);
        T = 1.0*(dataRans(1:1:end,5)-1)+1;
        h4=plot(dataRans(1:1:end,1)-0.15, T,'b:','MarkerSize',MarkerSize,'LineWidth',LineWidth);
        
        
        filename2 = sprintf('Results/%s_%s_modNew_const/',cas,mod);
        dataRans = ReadRansX(filename2,4);
        T = 1.0*(dataRans(1:1:end,5)-1)+1;
        h5=plot(dataRans(1:1:end,1)-0.15, T,'b-','MarkerSize',MarkerSize,'LineWidth',LineWidth);
        
        %%no Mod
        if(y_plus==1)
            filename2 = sprintf('Results/%s_%s_noMod_yplus_const/',cas,mod);
            dataRans = ReadRansX(filename2,4);
            T = 1.0*(dataRans(1:1:end,5)-1)+1;
            h9=plot(dataRans(1:1:end,1)-0.15, T,'r-','MarkerSize',MarkerSize,'LineWidth',LineWidth);
        else
        
            filename2 = sprintf('Results/%s_%s_noMod_const/',cas,mod);
            dataRans = ReadRansX(filename2,4);
            T = 1.0*(dataRans(1:1:end,5)-1)+1;
            h9=plot(dataRans(1:1:end,1)-0.15, T,'r-','MarkerSize',MarkerSize,'LineWidth',LineWidth);
        end
        
        %if(j==1)
            filename2 = sprintf('Results/Cluster/%s_%s_noMod_Prt/',cas,mod);
            dataRans = ReadRansX(filename2,4);
            T = 1.0*(dataRans(1:1:end,5)-1)+1;
            h8=plot(dataRans(1:1:end,1)-0.15, T,'r:','MarkerSize',MarkerSize,'LineWidth',LineWidth);

            if(y_plus==1)

% NOOOOT SOLVED YET IN THE CLUSTER
                filename2 = sprintf('Results/%s_%s_noMod_yplus_DWX/',cas,mod);
                dataRans = ReadRansX(filename2,4);
                T = 1.0*(dataRans(1:1:end,5)-1)+1;
                h6=plot(dataRans(1:1:end,1)-0.15, T,'rx','MarkerSize',MarkerSize,'LineWidth',LineWidth);
            else

                filename2 = sprintf('Results/Cluster/%s_%s_noMod_DWX/',cas,mod);
                dataRans = ReadRansX(filename2,4);
                T = 1.0*(dataRans(1:1:end,5)-1)+1;
                h6=plot(dataRans(1:1:end,1)-0.15, T,'rx','MarkerSize',MarkerSize,'LineWidth',LineWidth);
                if(ncases>1)
                    filename2 = sprintf('Results/Cluster/%s_%s_noMod_DWXZh/',cas,mod); %filename2 = sprintf('Results/Cluster/%s_%s_noMod_NK/',cas,mod);
                    dataRans = ReadRansX(filename2,4);
                    T = 1.0*(dataRans(1:1:end,5)-1)+1;
                    h7=plot(dataRans(1:1:end,1)-0.15, T,'r-.','MarkerSize',MarkerSize,'LineWidth',LineWidth);
                end
            
            end
            filename2 = sprintf('%s/',mod);
             dataRans = ReadRansX(filename2,4);
             T = 1.0*(dataRans(1:1:end,5)-1)+1;
             %h10=plot(dataRans(1:1:end,1)-0.15, T,'g-','MarkerSize',MarkerSize,'LineWidth',LineWidth);
            
            if(ncases==1)
                %legend([h2 h4 h5 h6 h8 h9 h10 h1],'New+DWX','New+Prt','New+cPrt','Conv+DWX','Conv+Prt','Conv+cPrt','Simulating','DNS','Location','northwest');
            else
                %legend([h2 h3 h4 h5 h6 h7 h8 h9 h10 h1],'New+DWX','New+DWX+Zhang','New+Prt','New+cPrt','Conv+DWX','Conv+DWX+Zhang','Conv+Prt','Conv+cPrt','Simulating','DNS','Location','northwest');
                %legend([h2 h3 h4 h5 h6 h7 h8 h9 h1],'New+DWX','New+DWX+Zhang','New+Prt','New+cPrt','Conv+DWX','Conv+DWX+Zhang','Conv+Prt','Conv+cPrt','DNS','Location','northwest');

            end


        if((i==3) || (i==4))
            xlabel('Streamwise distance, $x/D$','Interpreter','latex'); 
        end
        if((i==1) || (i==3))
            ylabel('Wall temperature, $T/T_0$','Interpreter','latex');
        end
        if(i<2) 
            %legend([h2 h3 h4 h1 h11],'Conventional','Present Study','Density corrections','DNS-Nemati','DNS-Bae','Location','northwest');
            %legend([h2 h3 h4 h1 h11],'Conventional','Present Study','Present Study with term','DNS-Nemati','DNS-Bae','Location','northwest');
        end
       title(sprintf(titleSpec,turbMod_name{i})); %,'FontWeight','bold');
        titletext = get(gca,'title');

        %set(get(gca,'title'),'Position',titletext.Position + [0.0 0.5 0.0])
        % legend([h3 h2  h1],'Original ystar','"Semi-local" model','DNS','Location','northwest');
        set(gca,'fontsize', 24)

        set(gca, 'XMinorTick', 'on')
        % set(gca, 'XTick', [0:5:30])

        set(gca, 'YMinorTick', 'on')
        % set(gca, 'YTick', [0:5:30])

        ymax = tmax + (tmax-1);
        if(ymax<1.2)
            ymax=1.2;
        end
        xlim([0 60]);
        ylim([1 ymax]);


         %title

    %     titletext.Position = titletext.Position + [0.1 0.0 0.0]
        %if(plot_var==1);

        %else
         %   set(get(gca,'title'),'Position',titletext.Position + [0.0 0.04 0.0])
        %end

    end

    % saveas(gcf,'MK_orig_caseA.png')
    set(gcf,'PaperUnits','points');
    factor1 = 0.85; %0.78;
    factor2 = 0.85; %0.78;
    %set(gcf,'PaperPosition',factor*[0 0 1400 800]);
    set(gcf,'PaperPosition',[0 0 factor1*1400 factor2*800]);
    %filename  = sprintf('%s_%s_bothDNS.eps',cas,data);
    filename  = sprintf('%s_%s.eps',cas, turbMod_name{nturb});
    print(filename,'-depsc')
end