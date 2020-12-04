clear all

fmt='real*8';
accuracy='real*8';
Ieee='b';


%---> Parametres

nx=60;
ny=140;
nz=27;
onedegx=nx/60;
onedegy=ny/140;

lice  = 10 ; %notice that this is number of ice grid points plus one (because first on is land)
           % 3 for 2deg of ice, 5 for 4 deg, 7 for 6 deg, 9 for 8 deg
Qice=4.0;


lati=1*((0:ny)-ny/2);

%================================================
%=============>   TOPOGRAPHY    <================
%================================================

topo = zeros(nx,ny);
topo(:,:) = -4000. ;

%south
topo(:,1) = 0.;
%north
topo(:,ny) = 0.;
%America  
topo(1 ,23:ny) = 0.;
%topo(nx,12:ny) = 0.;
%Europe+Africa
%topo(30,12:ny) = 0.;
%Asia
%topo(30:nx,ny-8:ny) = 0.;
%Scotia Ridge
%topo(1 ,2:22) = -3000.; %!cyc: remove the sill 
%topo(nx,2:14) = -3500.;

mask = topo;

for j=1:ny
    for i=1:nx
        if mask(i,j) == 0
            mask(i,j) = 0. ; 
        else
            mask(i,j) = 1. ; 
        end
    end
end


fid=fopen('topog.bin','w','b'); fwrite(fid,topo,'real*4'); fclose(fid);


figure(1)
pcolor(topo'), shading flat; colorbar
%pcolor(mask), shading flat; colorbar

%================================================
%=============>   Wind Stress    <===============
%================================================

wind  = zeros(nx,ny);
%tau0=0.17;
tau0=0.2; %!cyc: set the maximum wind stress the same as NV11 and M15

for j=1:22
	y=j-2;
	wind(:,j)=tau0*(sin(pi*y/20));
end

% for j=1:40
% 	y=j-1.2;
% 	wind(:,j)=1.007*tau0*(sin(3.5*pi*y/(ny-2)));
% end
% for j=41:60
% 	y=j-40.8;
% 	wind(:,j)=-0.65*tau0*(sin(3.5*pi*y/(ny-2)));
% end
% for j=61:80
% 	wind(:,j)=-0.65*tau0;
% end
% for j=80:100
% 	y=j-100.2;
% 	wind(:,j)=0.65*tau0*(sin(3.5*pi*y/(ny-2)));
% end
% for j=101:140
% 	y=j-119;
% 	wind(:,j)=0.58*tau0*(cos(3.5*pi*(y)/(ny-3)));
% end
% 
% for j=46:94
% 	y=j-70.5;
% 	wind(:,j)=wind(:,j)+0.56*tau0*(cos(3*pi*y/(ny+0.5)))^2;
% end
% 
% for i=1:nx
%   wind(i,:)=smooth(wind(i,:),5);
% end

wind(:,1:2)= 0;
wind(:,139:140)= 0;

fid=fopen('wind_x.bin','w','b'); fwrite(fid,wind,'real*4'); fclose(fid);

figure(2)
%pcolor(wind'),shading flat; colorbar
%colorbar
plot(wind(2,:))
hold off
%stop

%%
%=====================================================
%=================>  GM2D map     <===================
%=====================================================
% 
% GM2D = (1/3)+(2/3)*abs(wind/max(wind(10,:)));
% pcolor(GM2D'),shading flat; colorbar
% 
% fid=fopen('GM_bol2dFile.bin','w','b'); fwrite(fid,wind,'real*4'); fclose(fid);

%=====================================================
%=================>      SST      <===================
%=====================================================

SST = zeros(nx,ny);

%theta = ny-3; %width

%for i=1:nx
% for j=lice+1:ny-1
%    y=j-(ny+1)/2;
%	SST(i,j) = mask(i,j)*(29*cos(pi*y/theta));
% end
% SST(i,lice+1:ceil(ny/2))=SST(i,lice+1:ceil(ny/2))...
%                        +(6-SST(i,lice+1))*linspace(1,0,ceil(ny/2)-lice).^1.5;
%end

SST(:,1:22)=0;            %!cyc: set boundary condition to fixed-flux for the entire channel
SST(:,23:end)=10.0; %!cyc: set boundary condition to a constant SST for the basin
SST=repmat(SST,[1,1,nz]);

fid=fopen('rbcs_T.bin','w','b'); fwrite(fid,SST,'real*4'); fclose(fid);


figure(5)
plot(SST(10,:,1));
hold off;

SSTmask=zeros(nx,ny,nz);
SSTmask(:,:,1)=mask;

SSTmask(:,[1:lice,end],1)=0;
SSTmask(:,[lice+1:22,end],1)=0; %!cyc: set boundary condition to fixed-flux for the rest of the channel 

fid=fopen('rbcs_mask.bin','w','b'); fwrite(fid,SSTmask,'real*4'); fclose(fid);

figure(25)
clf
subplot(2,1,1)
imagesc(SSTmask(:,:,1));
subplot(2,1,2)
imagesc(squeeze(SSTmask(5,:,:))');



%====================================================================
%===================> Surface heat flux <============================
%====================================================================


Qc=zeros(nx,ny);
Qc(:,2:lice)=Qice;
Qc(:,lice+1:22)=0; %!cyc: set zero flux boundary condition for the rest of the channel 

fid=fopen('Qo','w','b'); fwrite(fid,Qc,'real*4'); fclose(fid);


figure(12)
imagesc(Qc(:,:));
hold off;


% =============== Taken out all tracer and float stuff ==================
% %====================================================================
% %=======================>   ptracers stuff  <========================
% %====================================================================
% 
% 
% genparam
% z = 1:length(zc);
% z = zc(1,1,:);
% y = 1:length(squeeze(yc(1,:,1)));
% y = yc(1,:,1);
% 
% x = 1:length(xc(:,1,1));
% x = xc(:,1,1);
% nx=105;
% ny=70;
% nz=40;
% 
% % dye
% 
% 
% dye=zeros(nx,ny,nz);
% dye(2:5,ny-12:ny-8,37:38)=1;
% fid=fopen('dye1.bin','w','b'); fwrite(fid,dye,'real*4'); fclose(fid);
% figure(71)
% subplot(211)
% pcolor(x,y,dye(:,:,37)');shading flat; colorbar
% subplot(212)
% D(:,:)=dye(3,:,:);
% pcolor(y,z,D');shading flat; colorbar
% 
% 
% %dye=zeros(nx,ny,nz);
% %dye(10:20,55:60,37:37)=1;
% %fid=fopen('dye.bin','w','b'); fwrite(fid,dye,'real*4'); fclose(fid);
% %figure(7)
% %subplot(211)
% %pcolor(x,y,dye(:,:,37)');shading flat; colorbar
% %subplot(212)
% %D(:,:)=dye(15,:,:);
% %pcolor(y,z,D');shading flat; colorbar
% 
% %================================================================
% %=====================>  Floats initial <========================
% %================================================================
% 
% %fid=fopen('flt_ini_pos.bin','r','b');
% %FLTini = fread(fid,'real*8');
% Nfloatsx = 5;
% Nfloatsy = 5;
% Nfloats = Nfloatsx*Nfloatsy;
% 
% Xi = 1:Nfloats;
% Yj = 1:Nfloats;
% Xi(:)=0;
% Yj(:)=0;
% dxfloats = 60/Nfloatsx;
% dyfloats = 10/Nfloatsy;
% 
% nn=1;
% for ii=1:Nfloatsx
% for jj=1:Nfloatsy
%     Xi(nn) = (ii-1)*dxfloats+dxfloats/2;
%     Yj(nn) = 130+(jj-1)*dyfloats+dyfloats/2;
%     nn=nn+1;
% end
% end
% 
% FLTini = zeros( 9, Nfloats+1 );
% 
% %FLTini(1) = Nfloats ;
% %- fill-up 1rst reccord
% FLTini(1,1) = Nfloats;
% FLTini(6,1) = Nfloats;
% 
% for i=2:Nfloats+1
%     
%    n = i-1;
%    
% %c     - npart   A unique float identifier (1,2,3,...)
% 
%    FLTini(1,i)=n ;
% 
% %c     - tstart  start date of integration of float (in s)
% %c               Note: If tstart=-1 floats are integrated right from the 
% %c               beginning
% 
%    FLTini(2,i)= -1 ;
% 
% %c     - xpart   x position of float (in units of XC)
%    
%    FLTini(3,i)= Xi(n) ;
% 
% %c     - ypart   y position of float (in units of YC)
%    
%    FLTini(4,i)= Yj(n) ;
% 
% %c     - kpart   actual vertical level of float
%    
%    FLTini(5,i)= -500.-0.5*n;
% 
% %c     - kfloat  target level of float (should be the same as kpart at
% %c               the beginning)
%    
%    FLTini(6,i)= -500.-0.5*n;
% 
% %c     - iup     flag if the float 
% %c               - should profile   ( >  0 = return cycle (in s) to surface) 
% %c               - remain at depth  ( =  0 )
% %c               - is a 3D float    ( = -1 ).
% %c               - should be advected WITHOUT additional noise ( = -2 ). 
% %c                 (This implies that the float is non-profiling)
% %c               - is a mooring     ( = -3 ), i.e. the float is not advected
%    
%    FLTini(7,i)=-1 ;
% 
% %c     - itop    time of float the surface (in s)
%    
%    FLTini(8,i)=0 ;
% 
% %c     - tend    end  date of integration of float (in s)
% %c               Note: If tend=-1 floats are integrated till the end of
% %c               the integration   
% 
%    FLTini(9,i)=-1 ;
%    
% 
% end
% 
% fid=fopen('flt_ini_pos.bin','w','b'); fwrite(fid,FLTini,'real*4');
% fclose(fid);
%
%
% ====================      end of trancer and float stuff ===============
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                 START FROM U V T & S FILES
% 
% 
%     
% Tz1(:,:,:) = zeros(nx,ny,nz) +25 + 0.1*randn(nx,ny,nz);
% 
% Tz=Tz1;
% for k=1:nz
%  Tz(:,:,k) = Tz1(:,:,k).*mask(:,:);
% end
% 
% figure(111)
% clf
% imagesc(Tz(:,:,1));shading flat; colorbar;
% fid=fopen('T.init','w','b'); fwrite(fid,Tz,'real*4'); fclose(fid);

% Uz(:,:,:) = zeros(nx,ny,nz);
% Uz1 =rdmds(['U.0019500000']);
% 
% for i=1:nx-2
% for j=1:ny
%     ii=i*2;
%     jj=j*2;
%     Uz(i,j,:)=Uz1(ii,jj,:)*mask(i,j);
% end
% Uz(nx-1,j,:)=Uz(nx-2,j,:);
% Uz(nx,j,:)=Uz(nx-2,j,:);
% end
% figure(112)
% pcolor(Uz(:,:,1)');shading flat; colorbar;
% fid=fopen('U.init','w','b'); fwrite(fid,Uz,'real*4'); fclose(fid);
% 
% 
% Vz(:,:,:) = zeros(nx,ny,nz);
% Vz1 =rdmds(['V.0019500000']);
% 
% for i=1:nx-2
% for j=1:ny
%     ii=i*2;
%     jj=j*2;
%     Vz(i,j,:)=Vz1(ii,jj,:)*mask(i,j);
% end
% Vz(nx-1,j,:)=Vz(nx-2,j,:);
% Vz(nx,j,:)=Vz(nx-2,j,:);
% end
% figure(112)
% pcolor(Vz(:,:,1)');shading flat; colorbar;
% fid=fopen('V.init','w','b'); fwrite(fid,Vz,'real*4'); fclose(fid);

% 
% % % !! ADD in data !!
% %  uVelInitFile=    'U.init',
% %  vVelInitFile=    'V.init',
% %  hydrogThetaFile= 'T.init',
% %  hydrogSaltFile = 'S.init',


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

