
clc

% 'VARIABLES ="Z","Rad","U","W",                 (1-4)
%             "C","T","k","eps",               (5-8)
%             "v2","omega","nuSA","yplus",     (9-12)
%             "RHO","Pe","mu","mut",           (13-16)
%             "lamcp","cp","alphat","kt",      (17-20)
%             "epst"'                          (21)
plot_var = 19; 
ncore = 4;
imax = 96; kmax=96;

rmax = 0.50005; rmin = 0.0;

%% plotted case:
cas  = 'caseA';   % cAL,      cBL,   cCL, cDL, CEL
mod  = 'MK';      % MK,        VF,    OM,  SA, 
vers = 'm';       % noMod, modNew, Aupoix
Tmod = 'DWX';     % const, Prt, DWX, NK


%% reading data:
if(vers=='m')
    filename  = sprintf('Results/%s_%s_%s/',mod,cas,Tmod);
else
    filename  = sprintf('Results/%s%s_%s_%s/',mod,vers,cas,Tmod);
end

filename  = sprintf('Results/MK/',mod,cas,Tmod);
filename2 = sprintf(filename);
data = readTecplot(filename2,ncore, imax, kmax);


%% ploting data
i=2+1;
figure(i)
contourf(data(:,:,1),data(:,:,2),data(:,:,plot_var),20)
colorbar

%caxis([0. 0.01])        
axis([0 60 rmin rmax])  % 0 x/L and 0 r

%%+z=data(:,:,1);
z=data(:,:,1);              r=data(:,:,2);
T=1.0*(data(:,:,6)-1)+1;    H=data(:,:,5);
Uz=data(:,:,4);             Ur=data(:,:,3);
k=data(:,:,7);              eps=data(:,:,8);
kt=data(:,:,20);            epst=data(:,:,21);
alphat=data(:,:,19);
n= size(z,1);
%return

i=i+1;
figure(i); 
plot(z(1,:),alphat(1,:),'rx'); hold on;
plot(z(n,:),alphat(n,:),'bo'); 
xlabel('Axial direction'); 
ylabel('Eddy diffusivity, alpha t'); 
legend('Centerline value', 'Wall value')

i=i+1;
figure(i); 
plot(z(1,:),T(1,:),'rx'); hold on;
plot(z(n,:),T(n,:),'bo'); hold off;
xlabel('Axial direction'); 
ylabel('Temperature, T'); 
legend('Centerline value', 'Wall value')

i=i+1;
figure(i); 
plot(z(1,:),k(1,:),'rx'); hold on;
plot(z(n,:),k(n,:),'bo'); hold off;
xlabel('Axial direction'); 
ylabel('Turb kin En., k'); 
legend('Centerline value', 'Wall value')

i=i+1;
figure(i); 
plot(z(1,:),kt(1,:),'rx'); hold on;
plot(z(n,:),kt(n,:),'bo'); hold off;
xlabel('Axial direction'); 
ylabel('Fluctuating temp, kt'); 
legend('Centerline value', 'Wall value')

i=i+1;
figure(i); 
plot(z(1,:),eps(1,:),'rx'); hold on;
plot(z(n,:),eps(n,:),'bo'); hold off;
xlabel('Axial direction'); 
ylabel('Dissipation of Turb kin En.. eps'); 
legend('Centerline value', 'Wall value')

i=i+1;
figure(i); 
plot(z(1,:),epst(1,:),'rx'); hold on;
plot(z(n,:),epst(n,:),'bo'); hold off;
xlabel('Axial direction'); 
ylabel('Dissipation of Fluctuating temp, et'); 
legend('Centerline value', 'Wall value')

i=i+1;
figure(i); 
plot(z(1,:),H(1,:),'rx'); hold on;
plot(z(n,:),H(n,:),'bo'); hold off;
xlabel('Axial direction'); 
ylabel('Enthalpy, H'); 
legend('Centerline value', 'Wall value')

i=i+1;
figure(i); 
plot(z(1,:),Uz(1,:),'rx'); hold on;
plot(z(n,:),Uz(n,:),'bo'); hold off;
xlabel('Axial direction'); 
ylabel('Velocity'); 
legend('Centerline value', 'Wall value')

%%
filename  = sprintf('Results/MK_caseA_Prt/',mod,cas,Tmod);
filename2 = sprintf(filename);
data = readTecplot(filename2,ncore, imax, kmax);

z=data(:,:,1);              r=data(:,:,2);
T=1.0*(data(:,:,6)-1)+1;    H=data(:,:,5);
Uz=data(:,:,4);             Ur=data(:,:,3);
k=data(:,:,7);              eps=data(:,:,8);
kt=data(:,:,20);            epst=data(:,:,21);
alphat=data(:,:,19);
n= size(z,1);

i=2+1+1

figure(i); 
plot(z(1,:),alphat(1,:),'r-');
plot(z(n,:),alphat(n,:),'b-'); 
xlabel('Axial direction'); 
ylabel('Eddy diffusivity, alpha t'); 
legend('Centerline value', 'Wall value')