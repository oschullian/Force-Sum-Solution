clear all
close all
nzall=[101,201,301,401]
readmax=10000000;

for runindex=1
for inz=1:4
fileending=sprintf('_ee_%i.txt',runindex);
pathgro=sprintf('init%i.gro',runindex);
pathbox=sprintf('box_ee_%i.xvg',runindex);
pathforce=sprintf('force_ee_%i.xvg',runindex);
pathvel=sprintf('vel_ee_%i.xvg',runindex);
pathpos=sprintf('pos_ee_%i.xvg',runindex);

fr=fopen(pathgro,'r');
dat=fgetl(fr);
n=str2num(fgetl(fr));

mass=zeros(n,1);

for i=1:n
    i/n
dat=fgetl(fr);
dat=dat(10:end);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
dat=dat(~isspace(dat));
val=dat(1);

if (val=='H')
    mass(i)=1.008;
elseif (val=='C')
    mass(i)=12.011;
elseif (val=='N')
    mass(i)=14.0067;
elseif (val=='O')
    mass(i)=15.9994;
elseif (val=='P')
    mass(i)=30.9738;
else
    val
    error('what element?')
    
end

end

fclose(fr);
% define constants
xfa=1e-9;
vfa=1e-9/1e-12;
ffa=1.66053892103219e-12;
mfa=1.66e-27;

% define the z-vector
fbox=fopen(pathbox,'r');
for i=1:24
   fgetl(fbox) ;
end
ba=transpose(fscanf(fbox,'%f',[7,Inf]));
fclose(fbox);


nz=nzall(inz);
z=linspace(-0.5*max(ba(:,4)),0.5*max(ba(:,4)),nz);
dz=z(2)-z(1);

fz=fopen(horzcat(sprintf('zvec_nz=%i',nz),fileending),'w');

fprintf(fz,'%e ',z);
fclose(fz);

m_unique=unique(mass);

% open files and define the number of frames 
nt=1;

% replicate the mass 
m_part=repmat(mass,1,nt);
m_index=m_part;
for i=1:length(m_unique)
    m_index(m_part==m_unique(i))=i;
end
% define output arrays
dens=zeros(length(m_unique),nz);

%
bap=fopen(pathbox,'r');
for i=1:24
   fgetl(bap);
end
fa=fopen(pathforce,'r');
for i=1:n*3+18
   fgetl(fa);
end
va=fopen(pathvel,'r');
for i=1:n*3+18
   fgetl(va);
end
xa=fopen(pathpos,'r');
for i=1:n*3+18
   fgetl(xa);
end


fkxx=fopen(horzcat(sprintf('ksigmaxx_nz=%i',nz),fileending),'w');
fkxy=fopen(horzcat(sprintf('ksigmaxy_nz=%i',nz),fileending),'w');
fkxz=fopen(horzcat(sprintf('ksigmaxz_nz=%i',nz),fileending),'w');
fkyy=fopen(horzcat(sprintf('ksigmayy_nz=%i',nz),fileending),'w');
fkyz=fopen(horzcat(sprintf('ksigmayz_nz=%i',nz),fileending),'w');
fkzz=fopen(horzcat(sprintf('ksigmazz_nz=%i',nz),fileending),'w');
fcxz=fopen(horzcat(sprintf('csigmaxz_nz=%i',nz),fileending),'w');
fcyz=fopen(horzcat(sprintf('csigmayz_nz=%i',nz),fileending),'w');
fczz=fopen(horzcat(sprintf('csigmazz_nz=%i',nz),fileending),'w');

fd1=fopen(horzcat(sprintf('dens1_nz=%i',nz),fileending),'w');
fd2=fopen(horzcat(sprintf('dens2_nz=%i',nz),fileending),'w');
if (length(m_unique)>2)
    fd3=fopen(horzcat(sprintf('dens3_nz=%i',nz),fileending),'w');
    fd4=fopen(horzcat(sprintf('dens4_nz=%i',nz),fileending),'w');
    fd5=fopen(horzcat(sprintf('dens5_nz=%i',nz),fileending),'w');
end





framecounter=0;
run2=0;
for readin=1:readmax
    readin
    kinsigma=zeros(nz,3,3);
confsigma=zeros(nz,3);

dens=zeros(length(m_unique),nz);


% read in data 
bp=fscanf(bap,'%f',[7,nt]);
xp=fscanf(xa,'%f',[3*n+1,nt]);
vp=fscanf(va,'%f',[3*n+1,nt]);

fp=fscanf(fa,'%f',[3*n+1,nt]);
z_part=xp(4:3:end,:);

if (sum(size(z_part))==0)
        break 
end


framecounter=framecounter+size(z_part,2);
v_x=vp(2:3:end,:)*vfa;
v_y=vp(3:3:end,:)*vfa;
v_z=vp(4:3:end,:)*vfa;

f_x=fp(2:3:end,:)*ffa;
f_y=fp(3:3:end,:)*ffa;
f_z=fp(4:3:end,:)*ffa;

b_center=0*bp(4,:);
b_center=repmat(b_center,n,1);

b_disp=bp(4,:);
b_disp=repmat(b_disp,n,1);

area=bp(2,:).*bp(3,:);
area=repmat(area,n,1)*xfa^2;

vol=area*dz*xfa;


% loop over displacements
for k=-2:2

z_ind_center=round((z_part+b_center+k*b_disp-z(1))/dz)+1;
z_ind_integral=ceil((z_part+b_center+k*b_disp-z(1))/dz)+1;

z_ind_center(z_ind_center>nz)=nan;
z_ind_integral(z_ind_integral>nz)=nan;

z_ind_center(z_ind_center<1)=nan;
z_ind_integral(z_ind_integral<1)=nan;


ind_center_cons=find(~isnan(z_ind_center));
ind_integral_cons=find(~isnan(z_ind_integral));


for i=1:length(ind_center_cons)
    ind_cons=ind_center_cons(i);
    m_ind=m_index(ind_cons);
    z_ind=z_ind_center(ind_cons);
    dens(m_ind,z_ind)=dens(m_ind,z_ind)+mfa*m_unique(m_ind)/vol(ind_cons);
    
    kinsigma(z_ind,1,1)=kinsigma(z_ind,1,1)-mfa*m_unique(m_ind)*v_x(ind_cons)*v_x(ind_cons)/vol(ind_cons);
    kinsigma(z_ind,2,2)=kinsigma(z_ind,2,2)-mfa*m_unique(m_ind)*v_y(ind_cons)*v_y(ind_cons)/vol(ind_cons);
    kinsigma(z_ind,3,3)=kinsigma(z_ind,3,3)-mfa*m_unique(m_ind)*v_z(ind_cons)*v_z(ind_cons)/vol(ind_cons);
    kinsigma(z_ind,2,3)=kinsigma(z_ind,2,3)-mfa*m_unique(m_ind)*v_y(ind_cons)*v_z(ind_cons)/vol(ind_cons);
    kinsigma(z_ind,1,2)=kinsigma(z_ind,1,2)-mfa*m_unique(m_ind)*v_x(ind_cons)*v_y(ind_cons)/vol(ind_cons);
    kinsigma(z_ind,1,3)=kinsigma(z_ind,1,3)-mfa*m_unique(m_ind)*v_x(ind_cons)*v_z(ind_cons)/vol(ind_cons);

end


for i=1:length(ind_integral_cons)
    ind_cons=ind_integral_cons(i);
    z_ind=z_ind_integral(ind_cons);
    
    confsigma(z_ind,1)=confsigma(z_ind,1)+f_x(ind_cons)/area(ind_cons);
    confsigma(z_ind,2)=confsigma(z_ind,2)+f_y(ind_cons)/area(ind_cons);
    confsigma(z_ind,3)=confsigma(z_ind,3)+f_z(ind_cons)/area(ind_cons);
end




end

for i=1:3
    confsigma(:,i)=cumsum(confsigma(:,i));
end



fprintf(fkxx,'%e ',kinsigma(:,1,1));
fprintf(fkxx,'\n');
fprintf(fkyy,'%e ',kinsigma(:,2,2));
fprintf(fkyy,'\n');
fprintf(fkzz,'%e ',kinsigma(:,3,3));
fprintf(fkzz,'\n');
fprintf(fkyz,'%e ',kinsigma(:,2,3));
fprintf(fkyz,'\n');
fprintf(fkxy,'%e ',kinsigma(:,1,2));
fprintf(fkxy,'\n');
fprintf(fkxz,'%e ',kinsigma(:,1,3));
fprintf(fkxz,'\n');


fprintf(fcxz,'%e ',confsigma(:,1));
fprintf(fcxz,'\n');
fprintf(fcyz,'%e ',confsigma(:,2));
fprintf(fcyz,'\n');
fprintf(fczz,'%e ',confsigma(:,3));
fprintf(fczz,'\n');

fprintf(fd1,'%e ',[m_unique(1), dens(1,:)]);
fprintf(fd1,'\n');
fprintf(fd2,'%e ',[m_unique(2), dens(2,:)]);
fprintf(fd2,'\n');
if (length(m_unique)>2)
  fprintf(fd3,'%e ',[m_unique(3), dens(3,:)]);
  fprintf(fd3,'\n');
  fprintf(fd4,'%e ',[m_unique(4), dens(4,:)]);
  fprintf(fd4,'\n'); 
  fprintf(fd5,'%e ',[m_unique(5), dens(5,:)]);
  fprintf(fd5,'\n'); 
end



end

fclose all;
end
end