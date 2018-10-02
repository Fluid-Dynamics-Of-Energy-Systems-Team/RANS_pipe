clear all;
px = 4; % number of processors
imax = 96; kmax = imax*3;   % imax*1 = 96, imax*2 = 192, imax*4 = 384, 
nvar = 19; %for tecplot
model = 'SA';  style= 'b-';
kstart = 1; %kmax/12;

%'VARIABLES ="X","Y","U","W",           (1-4)
% "C","T","RHO","mu","mut",             (5-9)
% "lamcp","cp","dqTdx","dqHdx",         (10-13)     
% "dqTdr","dqHdr","diff_x","diff_r"     (14-17)
% ,"Reldiff_x","Reldiff_r"'             (18-19)  
                  
                    
rmax = 0.50005; rmin = 0.0;
amp = 100000;   % amplification factor for such small Qs!
filename  =  sprintf('%s/%s_%s/','Results',model,num2str(kmax));

data = readHeat(filename,px, imax, kmax,nvar); %64, 192);

X     = data(:,:,1);   R = data(:,:,2);
H     = data(:,:,5);   T = data(:,:,6);
lamcp = data(:,:,10); cp = data(:,:,11);

dqTdx  = data(:,:,12);     dqHdx = data(:,:,13);
dqTdr  = data(:,:,14);     dqHdr = data(:,:,15);
%diff_x_code = data(:,:,16);     diff_r_code = data(:,:,17);
%reldiff_x_code = data(:,:,18);  reldiff_r_code = data(:,:,19);

lam = lamcp.*cp;

nr = size(T,1);
nx = size(T,2);

%% derivatives
dTdx = zeros(nr,nx-1);   dTdr = zeros(nr-1,nx);
dHdx = zeros(nr,nx-1);   dHdr = zeros(nr-1,nx);

% derivatives in x
dTdx(:,:)   = (T(1:nr,1:nx-1)-T(1:nr,2:nx))./(X(1:nr,1:nx-1)-X(1:nr,2:nx));
dHdx(:,:)   = (H(1:nr,1:nx-1)-H(1:nr,2:nx))./(X(1:nr,1:nx-1)-X(1:nr,2:nx));
% derivatives in r
dTdr(:,:)   = (T(1:nr-1,1:nx)-T(2:nr,1:nx))./(R(1:nr-1,1:nx)-R(2:nr,1:nx));
dHdr(:,:)   = (H(1:nr-1,1:nx)-H(2:nr,1:nx))./(R(1:nr-1,1:nx)-R(2:nr,1:nx));


%% diffusion factors: lambda for T and lambda/cp for H
lamx   = zeros(nr,nx-1);   lamr = zeros(nr-1,nx);
lamcpx = zeros(nr,nx-1); lamcpr = zeros(nr-1,nx);

%average in x
lamx(:,:)   = 0.5*(lam(1:nr,1:nx-1)  +lam(1:nr,2:nx));
lamcpx(:,:) = 0.5*(lamcp(1:nr,1:nx-1)+lamcp(1:nr,2:nx));
%average in r
lamr(:,:)   = 0.5*(lam(1:nr-1,1:nx)  +lam(2:nr,1:nx));
lamcpr(:,:) = 0.5*(lamcp(1:nr-1,1:nx)+lamcp(2:nr,1:nx));

%% diffusion term
qTx = lamx   .* dTdx *amp;   qTr = lamr   .* dTdr *amp;
qHx = lamcpx .* dHdx *amp;   qHr = lamcpr .* dHdr *amp;

% absolute difference
diffx = qTx-qHx;  diffr = qTr-qTr;
diffx2 = (qTx-qHx)./(qHx+1e-21);
diffr2 = (qTr-qHr)./(qHr+1e-21);

qTOTx= 0;
abs_diffx=0; rel_diffx = 0  ;
for i=1:nr
    for j=kstart:nx-1
        abs_diffx= abs_diffx+ abs(diffx(i,j));
        rel_diffx= rel_diffx+ abs(diffx2(i,j));
        qTOTx= qTOTx + qHx(i,j);
    end
end
abs_diffx
abs_diffx/qTOTx*100
%rel_diffx*100

qTOTr= 0;
abs_diffr=0; rel_diffr=0;
for i=1:nr-1
    for j=kstart:nx
       abs_diffr= abs_diffr+ abs(diffr(i,j));
       rel_diffr= rel_diffr+ abs(diffr2(i,j));
       qTOTr= qTOTr + qHr(i,j);
    end
end
abs_diffr
abs_diffr/qTOTr*100
%rel_diffr*100




%%from the code
diffx = dqTdx-dqHdx;  diffr = dqTdr-dqTdr;

qTOTx= 0;
abs_diffx=0; 
for i=1:nr
    for j=kstart:nx-1
        abs_diffx= abs_diffx+ abs(diffx(i,j));
        qTOTx= qTOTx + qHx(i,j);
    end
end
abs_diffx
abs_diffx/qTOTx*100
%rel_diffx*100

qTOTr= 0;
abs_diffr=0; 
for i=1:nr-1
    for j=kstart:nx
       abs_diffr= abs_diffr+ abs(diffr(i,j));
       qTOTr= qTOTr + qHr(i,j);
    end
end
abs_diffr
abs_diffr/qTOTr*100
%rel_diffr*100

return
%% ploting wall temperature
FSZ = 24;
LW = 3;
MSZ = 10;
TickLength = [0.025 0.05];
MarkerSize = 8;
LineWidth = 2;

dataJ = dlmread('../RANS_pipe/DNS_Hassan/Wall_Temp/Jbulk');   % A 60 long 
skipH = 100;
figure(10);
h1=plot(dataJ(1:skipH:end,1), dataJ(1:skipH:end,px),'ko','MarkerSize',MarkerSize,'LineWidth',LineWidth); hold on;
tmax = max(dataJ(1:end,px));

dataRans = ReadRansX(filename,4);
Tr = 1.0*(dataRans(:,5)-1)+1;
h2=plot(dataRans(1:1:end,1)-0.15, Tr,style,'MarkerSize',MarkerSize,'LineWidth',LineWidth);  hold on;
        size(dataRans(:,5)); 

xlim([0 60]);
ylim([1 1.3]);