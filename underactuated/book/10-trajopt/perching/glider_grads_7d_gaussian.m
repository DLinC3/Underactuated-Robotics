% ============================================================
% This function evaulates the gradients at a particular x,u
% ============================================================

function [dfdx,dfdu,dfdtheta] = glider_grads_7d_gaussian(t,x,u,param)

Sw=param.S_w;
Se=param.S_e;
l=param.l_h;
lw=param.l_w;
le=param.l_e;

m=param.m;
g=param.g;
rho=param.rho;
I=param.I;
dt=param.dt;

q3 = x(3);
q4 = x(4);
qdot10 = x(5);
qdot20 = x(6);
qdot3 = x(7);
qdot4 = 0*u;

c3=cos(q3);
c34=cos(q3+q4);
s3=sin(q3);
s34=sin(q3+q4);
s4=sin(q4);
c4=cos(q4);

vx0=param.vx0;
vz0=param.vz0;

qdot1=qdot10+vx0;
qdot2=qdot20+vz0;

xwdot = qdot1 - lw*qdot3*sin(q3);
zwdot=  qdot2 + lw*qdot3*cos(q3);
alpha_w = q3 - atan2(zwdot,xwdot);
Fw = rho*Sw*sin(alpha_w)*(zwdot*zwdot+xwdot*xwdot);

xedot = qdot1 + l*qdot3*sin(q3) + le*(qdot3+qdot4)*sin(q3+q4); 
zedot = qdot2 - l*qdot3*cos(q3) - le*(qdot3+qdot4)*cos(q3+q4);
alpha_e = q3+q4-atan2(zedot,xedot);
Fe = rho*Se*sin(alpha_e)*(zedot*zedot+xedot*xedot);

dxwdot_dq3=-lw*qdot3*cos(q3);

dxwdot_dqdot3=-lw*sin(q3);

dzwdot_dq3=-lw*qdot3*sin(q3);

dzwdot_dqdot3=lw*cos(q3);

U=1/(1+(zwdot/xwdot)*(zwdot/xwdot));

dalpha_w_dq3=1-U*(xwdot*dzwdot_dq3-zwdot*dxwdot_dq3)/(xwdot*xwdot);

dalpha_w_dqdot1=-U*(-zwdot)/(xwdot*xwdot);
dalpha_w_dqdot2=-U*(xwdot)/(xwdot*xwdot);
dalpha_w_dqdot3=-U*(xwdot*dzwdot_dqdot3-zwdot*dxwdot_dqdot3)/(xwdot*xwdot);

dxedot_dq3=l*qdot3*c3+le*(qdot3+qdot4)*c34;
dxedot_dq4=le*(qdot3+qdot4)*c34;

dxedot_dqdot3=l*sin(q3)+le*s34;
dxedot_dqdot4=le*s34;

dzedot_dq3=l*qdot3*sin(q3)+le*(qdot3+qdot4)*s34;
dzedot_dq4=le*(qdot3+qdot4)*s34;

dzedot_dqdot3=-l*c3-le*c34;
dzedot_dqdot4=-le*c34;

U=1/(1+(zedot/xedot)*(zedot/xedot));

dalpha_e_dq3=1-U*(xedot*dzedot_dq3-zedot*dxedot_dq3)/(xedot*xedot);
dalpha_e_dq4=1-U*(xedot*dzedot_dq4-zedot*dxedot_dq4)/(xedot*xedot);
dalpha_e_dqdot1=-U*(-zedot)/(xedot*xedot);
dalpha_e_dqdot2=-U*(xedot)/(xedot*xedot);
dalpha_e_dqdot3=-U*(xedot*dzedot_dqdot3-zedot*dxedot_dqdot3)/(xedot*xedot);
dalpha_e_dqdot4=-U*(xedot*dzedot_dqdot4-zedot*dxedot_dqdot4)/(xedot*xedot);

W=rho*Sw;
U=cos(alpha_w)*zwdot*zwdot;
P=sin(alpha_w)*2*zwdot;
Z=cos(alpha_w)*xwdot*xwdot;
Y=sin(alpha_w)*2*xwdot;

dFw_dq3=W*(U*dalpha_w_dq3+P*dzwdot_dq3+Z*dalpha_w_dq3+Y*dxwdot_dq3);

dFw_dqdot1=W*(U*dalpha_w_dqdot1+Z*dalpha_w_dqdot1+Y);
dFw_dqdot2=W*(U*dalpha_w_dqdot2+P+Z*dalpha_w_dqdot2);
dFw_dqdot3=W*(U*dalpha_w_dqdot3+P*dzwdot_dqdot3+Z*dalpha_w_dqdot3+Y*dxwdot_dqdot3);

W=rho*Se;
U=cos(alpha_e)*zedot*zedot;
P=sin(alpha_e)*2*zedot;
Z=cos(alpha_e)*xedot*xedot;
Y=sin(alpha_e)*2*xedot;

dFe_dq3=W*(U*dalpha_e_dq3+P*dzedot_dq3+Z*dalpha_e_dq3+Y*dxedot_dq3);
dFe_dq4=W*(U*dalpha_e_dq4+P*dzedot_dq4+Z*dalpha_e_dq4+Y*dxedot_dq4);
dFe_dqdot1=W*(U*dalpha_e_dqdot1+Z*dalpha_e_dqdot1+Y);
dFe_dqdot2=W*(U*dalpha_e_dqdot2+P+Z*dalpha_e_dqdot2);
dFe_dqdot3=W*(U*dalpha_e_dqdot3+P*dzedot_dqdot3+Z*dalpha_e_dqdot3+Y*dxedot_dqdot3);
dFe_dqdot4=W*(U*dalpha_e_dqdot4+P*dzedot_dqdot4+Z*dalpha_e_dqdot4+Y*dxedot_dqdot4);


[F H]=rbf_grads(q3,q4,qdot1,qdot2,param);%zeros(1,12);%

df5_dq3=(-(dFw_dq3*s3+Fw*c3+dFe_dq3*s34+Fe*c34)+F(1))/m;
df5_dq4=(-(dFe_dq4*s34+Fe*c34)+F(2))/m;
df5_dqdot1=(-(dFw_dqdot1*s3+dFe_dqdot1*s34)+F(3))/m;
df5_dqdot2=(-(dFw_dqdot2*s3+dFe_dqdot2*s34)+F(4))/m;
df5_dqdot3=-(dFw_dqdot3*s3+dFe_dqdot3*s34)/m;
df5_dqdot4=-(dFe_dqdot4*s34)/m;

df6_dq3=(dFw_dq3*c3-Fw*s3+dFe_dq3*c34-Fe*s34+F(5))/m;
df6_dq4=(dFe_dq4*c34-Fe*s34+F(6))/m;
df6_dqdot1=(dFw_dqdot1*c3+dFe_dqdot1*c34+F(7))/m;
df6_dqdot2=(dFw_dqdot2*c3+dFe_dqdot2*c34+F(8))/m;
df6_dqdot3=(dFw_dqdot3*c3+dFe_dqdot3*c34)/m;
df6_dqdot4=(dFe_dqdot4*c34)/m;

df7_dq3=(dFw_dq3*lw-dFe_dq3*(l*c4+le)+F(9))/I;
df7_dq4=(-dFe_dq4*(l*c4+le)+Fe*(l*s4)+F(10))/I;
df7_dqdot1=(dFw_dqdot1*lw-dFe_dqdot1*(l*c4+le)+F(11))/I;
df7_dqdot2=(dFw_dqdot2*lw-dFe_dqdot2*(l*c4+le)+F(12))/I;
df7_dqdot3=(dFw_dqdot3*lw-dFe_dqdot3*(l*c4+le))/I;
df7_dqdot4=(-dFe_dqdot4*(l*c4+le))/I; 

dfdx(5,1)=0;
dfdx(5,2)=0;    
dfdx(5,3)=df5_dq3;
dfdx(5,4)=df5_dq4;
dfdx(5,5)=df5_dqdot1;
dfdx(5,6)=df5_dqdot2;
dfdx(5,7)=df5_dqdot3;

dfdx(6,1)=0;
dfdx(6,2)=0;
dfdx(6,3)=df6_dq3;
dfdx(6,4)=df6_dq4;
dfdx(6,5)=df6_dqdot1;
dfdx(6,6)=df6_dqdot2;
dfdx(6,7)=df6_dqdot3;

dfdx(7,1)=0;
dfdx(7,2)=0;
dfdx(7,3)=df7_dq3;
dfdx(7,4)=df7_dq4;
dfdx(7,5)=df7_dqdot1;
dfdx(7,6)=df7_dqdot2;
dfdx(7,7)=df7_dqdot3;

dfdx(1,1)=0;
dfdx(1,2)=0;
dfdx(1,3)=0;
dfdx(1,4)=0;
dfdx(1,5)=1;
dfdx(1,6)=0;
dfdx(1,7)=0;


dfdx(2,1)=0;
dfdx(2,2)=0;
dfdx(2,3)=0;
dfdx(2,4)=0;
dfdx(2,5)=0;
dfdx(2,6)=1;
dfdx(2,7)=0;

dfdx(3,1)=0;
dfdx(3,2)=0;
dfdx(3,3)=0;
dfdx(3,4)=0;
dfdx(3,5)=0;
dfdx(3,6)=0;
dfdx(3,7)=1;

dfdx(4,1)=0;
dfdx(4,2)=0;
dfdx(4,3)=0;
dfdx(4,4)=0;
dfdx(4,5)=0;
dfdx(4,6)=0;
dfdx(4,7)=0;

dfdu(1)=0;
dfdu(2)=0;
dfdu(3)=0;
dfdu(4)=1;
dfdu(5)=df5_dqdot4;
dfdu(6)=df6_dqdot4;
dfdu(7)=df7_dqdot4;

dfdu=dfdu';

dfdtheta=H;

%dfdx=eye(length(dfdx))+dfdx*dt;


end

function [F H]=rbf_grads(q3,q4,qdot1,qdot2,param)

thetaCl=param.thetaCl;
thetaCd=param.thetaCd;
thetaCm=param.thetaCm;
muCl=param.muCl;
sigmaCl=param.sigmaCl;
muCd=param.muCd;
sigmaCd=param.sigmaCd;
muCm=param.muCm;
sigmaCm=param.sigmaCm;
S_w=param.S_w;

NB=length(muCl);

m=param.m;
g=param.g;
rho=param.rho;
I=param.I;

alpha = q3 - atan(qdot2/qdot1); % aoa

UU=1/(1+(qdot2*qdot2)/(qdot1*qdot1)); %

dalpha_dqdot1=-UU*(-qdot2/(qdot1*qdot1));%
dalpha_dqdot2=-UU*(1/qdot1);%
dalpha_dq3=1;%

V_norm = sqrt(qdot1^2+qdot2^2);% velocity norm

dVnorm_dqdot1=(0.5/(V_norm))*(2*qdot1);
dVnorm_dqdot2=(0.5/(V_norm))*(2*qdot2);

C=0.5*rho*S_w;

%display('glidergrads');

n_L1=-qdot2/V_norm;
n_L2= qdot1/V_norm;
n_D1=-qdot1/V_norm;
n_D2=-qdot2/V_norm;

dn_L1_dqdot1=(-1)*(-qdot2)/(V_norm^2)*dVnorm_dqdot1;
dn_L1_dqdot2=(-1)/V_norm+(-1)*(-qdot2)/(V_norm^2)*dVnorm_dqdot2;
dn_L2_dqdot1=1/V_norm+(-1)*qdot1/(V_norm^2)*dVnorm_dqdot1;
dn_L2_dqdot2=(-1)*qdot1/(V_norm^2)*dVnorm_dqdot2;

dn_D1_dqdot1=(-1)/V_norm+(-1)*(-qdot1)/(V_norm^2)*dVnorm_dqdot1;
dn_D1_dqdot2=(-1)*(-qdot1)/(V_norm^2)*dVnorm_dqdot2;
dn_D2_dqdot1=(-1)*(-qdot2)/(V_norm^2)*dVnorm_dqdot1;
dn_D2_dqdot2=(-1)/V_norm+(-1)*(-qdot2)/(V_norm^2)*dVnorm_dqdot2;

q=C*V_norm*V_norm;

dq_dqdot1=C*2*V_norm*dVnorm_dqdot1;
dq_dqdot2=C*2*V_norm*dVnorm_dqdot2;

Cm=0;
Cl=0;
Cd=0;

dCm_dq3=0;
dCm_dq4=0;
dCm_dqdot1=0;
dCm_dqdot2=0;

dCl_dq3=0;
dCl_dq4=0;
dCl_dqdot1=0;
dCl_dqdot2=0;

dCd_dq3=0;
dCd_dq4=0;
dCd_dqdot1=0;
dCd_dqdot2=0;

PHI=[];

for k=1:NB
    
ExCl=EXP(alpha,q4,muCl(1,k),muCl(2,k),sigmaCl(1));
ExCd=EXP(alpha,q4,muCd(1,k),muCd(2,k),sigmaCd(1));
ExCm=EXP(alpha,q4,muCm(1,k),muCm(2,k),sigmaCm(1));

dphiCl_dalpha=dGaussiandx(alpha,q4,muCl(1,k),muCl(2,k),sigmaCl(1),ExCl);
dphiCl_dq4=dGaussiandy(alpha,q4,muCl(1,k),muCl(2,k),sigmaCl(1),ExCl);
phiCl=gaussian(alpha,q4,muCl(1,k),muCl(2,k),sigmaCl(1),ExCl);

dphiCd_dalpha=dGaussiandx(alpha,q4,muCd(1,k),muCd(2,k),sigmaCd(1),ExCd);
dphiCd_dq4=dGaussiandy(alpha,q4,muCd(1,k),muCd(2,k),sigmaCd(1),ExCd);
phiCd=gaussian(alpha,q4,muCd(1,k),muCd(2,k),sigmaCd(1),ExCd);

dphiCm_dalpha=dGaussiandx(alpha,q4,muCm(1,k),muCm(2,k),sigmaCm(1),ExCm);
dphiCm_dq4=dGaussiandy(alpha,q4,muCm(1,k),muCm(2,k),sigmaCm(1),ExCm);
phiCm=gaussian(alpha,q4,muCm(1,k),muCm(2,k),sigmaCm(1),ExCm);


Cm=Cm+thetaCm(k)*phiCm;
Cl=Cl+thetaCl(k)*phiCl;
Cd=Cd+thetaCd(k)*phiCd;

AA=dphiCl_dalpha*dalpha_dq3;
BB=dphiCl_dalpha*dalpha_dqdot1;
CC=dphiCl_dalpha*dalpha_dqdot2;

dCl_dq3=dCl_dq3+thetaCl(k)*AA;
dCl_dq4=dCl_dq4+thetaCl(k)*dphiCl_dq4;
dCl_dqdot1=dCl_dqdot1+thetaCl(k)*BB;
dCl_dqdot2=dCl_dqdot2+thetaCl(k)*CC;

AA=dphiCd_dalpha*dalpha_dq3;
BB=dphiCd_dalpha*dalpha_dqdot1;
CC=dphiCd_dalpha*dalpha_dqdot2;

dCd_dq3=dCd_dq3+thetaCd(k)*AA;
dCd_dq4=dCd_dq4+thetaCd(k)*dphiCd_dq4;
dCd_dqdot1=dCd_dqdot1+thetaCd(k)*BB;
dCd_dqdot2=dCd_dqdot2+thetaCd(k)*CC;

AA=dphiCm_dalpha*dalpha_dq3;
BB=dphiCm_dalpha*dalpha_dqdot1;
CC=dphiCm_dalpha*dalpha_dqdot2;

dCm_dq3=dCm_dq3+thetaCm(k)*AA;
dCm_dq4=dCm_dq4+thetaCm(k)*dphiCm_dq4;
dCm_dqdot1=dCm_dqdot1+thetaCm(k)*BB;
dCm_dqdot2=dCm_dqdot2+thetaCm(k)*CC;

PHI=[PHI; phiCl];

end

PHI=[PHI;1];

%PHI

df5dtheta=[q*n_L1*PHI'/m q*n_D1*PHI'/m 0*PHI'/m];
df6dtheta=[q*n_L2*PHI'/m q*n_D2*PHI'/m 0*PHI'/m];
df7dtheta=[0*PHI'/I 0*PHI'/I q*PHI'/I];

H=[zeros(4,3*length(PHI));df5dtheta;df6dtheta;df7dtheta];

Cm=Cm+thetaCm(NB+1);
Cl=Cl+thetaCl(NB+1);
Cd=Cd+thetaCd(NB+1);

dFx_dq3=q*n_L1*dCl_dq3+q*n_D1*dCd_dq3;
dFx_dq4=q*n_L1*dCl_dq4+q*n_D1*dCd_dq4;
dFx_dqdot1=q*Cl*dn_L1_dqdot1+q*dCl_dqdot1*n_L1+dq_dqdot1*Cl*n_L1+q*Cd*dn_D1_dqdot1+q*dCd_dqdot1*n_D1+dq_dqdot1*Cd*n_D1;
dFx_dqdot2=q*Cl*dn_L1_dqdot2+q*dCl_dqdot2*n_L1+dq_dqdot2*Cl*n_L1+q*Cd*dn_D1_dqdot2+q*dCd_dqdot2*n_D1+dq_dqdot2*Cd*n_D1;

dFz_dq3=q*n_L2*dCl_dq3+q*n_D2*dCd_dq3;
dFz_dq4=q*n_L2*dCl_dq4+q*n_D2*dCd_dq4;
dFz_dqdot1=q*Cl*dn_L2_dqdot1+q*dCl_dqdot1*n_L2+dq_dqdot1*Cl*n_L2+q*Cd*dn_D2_dqdot1+q*dCd_dqdot1*n_D2+dq_dqdot1*Cd*n_D2;
dFz_dqdot2=q*Cl*dn_L2_dqdot2+q*dCl_dqdot2*n_L2+dq_dqdot2*Cl*n_L2+q*Cd*dn_D2_dqdot2+q*dCd_dqdot2*n_D2+dq_dqdot2*Cd*n_D2;

dM_dq3=q*dCm_dq3;
dM_dq4=q*dCm_dq4;
dM_dqdot1=dq_dqdot1*Cm+q*dCm_dqdot1;
dM_dqdot2=dq_dqdot2*Cm+q*dCm_dqdot2;

F=[dFx_dq3 dFx_dq4 dFx_dqdot1 dFx_dqdot2 dFz_dq3 dFz_dq4 dFz_dqdot1 dFz_dqdot2 dM_dq3 dM_dq4 dM_dqdot1 dM_dqdot2];


end

function a=dGaussiandx(x,y,mux,muy,sigma,Ex)

p=pi;

a=(2*mux - 2*x)/(4*p*sigma^2)*Ex;

end

function a=dGaussiandy(x,y,mux,muy,sigma,Ex)

p=pi;

a=(2*muy - 2*y)/(4*p*sigma^2)*Ex;

end

function phi=gaussian(x,y,mux,muy,sigma,Ex)
   
p=pi;

    phi=1/(2*p*sigma)*Ex;

end

function Ex=EXP(x,y,mux,muy,sigma)

Ex=exp(-((x-mux)^2+(y-muy)^2)/(2*sigma));

end
