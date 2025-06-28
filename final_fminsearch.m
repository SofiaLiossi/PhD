% fmin search for Smooth Nonlinear	

%lsqcurve fit
function final_fminsearch
initial=[4035584,5672242,22000,8000,15000,5000,15000,5000,4622,1111358,13206];

     function C=kineticsfmin(theta,t)
        
%initial=[4035584,5504723,13848,8000,14500,7500,13500,7200,1180,1070474,13206];

%19/8 exoume energa=15000 kroysmata ara E=30000
%14/8 exoume 20000=Id
opt=odeset('NonNegative',1:11);
%opt = odeset('RelTol',1e-15,'AbsTol',1e-15);
[T,Cv]=ode45(@DifEq1,t,initial,opt);

    function dC=DifEq1(t,c)
      
%parameters
m=1/(82*365);
L=176; %influx rate
em=0.79;%em=0.79; %efficacy of mask
%pd=0; %proportion who social distancing theta3
na=0.752; %mod parameter 
l=0.926; %sensitivy of self test
ev=0.75; %effectiveness of vaccines
v=0.0019; %(apo 14/08 ews 12/10)
psiu=1/180;
sa=1;
su=1/5.8;
%sicu=1/7;
%f1=0.5;
f2=4/100;
ff2=4/1000;
%f2=5.7/100;
%ff2=5.7/1000;
gaa=1/6;
gi=1/6;
gh=1/18;
p=0.254;
dh=1/18; %20.4 merew noshleias
sh=1/4;
a=0.17;

%variables
S=c(1);
U=c(2);
E1=c(3);
E2=c(4);
E=c(3)+c(4);
I1=c(5);
I2=c(6);
I=c(5)+c(6);
A1=c(7);
A2=c(8);
A=c(7)+c(8);
H1=c(9);
R=c(10);
D=c(11);

%theta1=b
%theta2=k pososto anixneyshs
%theta3=f1 pososto sympt
%theta4=pmask 

 N=S+U+E+I+A+H1+R;
  ls=(theta(1).*(1-theta(5)).*(1-em.*theta(4)).*(E+(1-theta(2)).*na.*A+p*I))./((N-((1-p)*I+theta(2).*A+H1+ev*U)));
  lu= (1-ev).*ls;
%disp(ls);

%differential equations
dcdt = zeros (11, 1);
dcdt(1)= L-ls.*S+psiu.*U-v.*S-m.*S; %s
dcdt(2)=-lu.*U-psiu.*U+v.*S-m.*U; %u
dcdt(3)=ls.*S-su.*E1-m.*E1; %e1
dcdt(4)=lu.*U-su.*E2-m.*E2; %e2

dcdt(5)=theta(3).*su.*E1-theta(6).*sh.*I1-(1-theta(6)).*gi.*I1+theta(2).*l.*sa.*A1-m.*I1;%i1
dcdt(6)=theta(3).*su.*E2-(theta(6)/10).*sh.*I2-(1-theta(6)/10).*gi.*I2+theta(2).*l.*sa.*A2-m.*I2;%i2

dcdt(7)=(1-theta(3)).*su.*E1-gaa.*A1-theta(2).*l.*sa.*A1-m.*A1; %A1
dcdt(8)=(1-theta(3)).*su.*E2-gaa.*A2-theta(2).*l.*sa.*A2-m.*A2;

dcdt(9)=theta(6).*sh.*I1+(theta(6)/10).*sh.*I2-a.*dh.*H1-(1-a).*gh.*H1-m.*H1; %h1

%dxdt(11)=(1-f3).*H1-f5.*gh.*Hnonicu-(1-f5).*dh.*Hnonicu-m.*Hnonicu;%hnonicu
%dxdt(11)=dxdt(9)-dxdt(10);%hnonicu

dcdt(10)=(1-theta(6)).*gi.*I1+(1-(theta(6))).*gi.*I2+gaa.*A+(1-a).*gh.*H1-m.*R; %recovered

dcdt(11)=a.*dh.*H1;% %deads

    dC=dcdt;
    end
C=Cv(:,11);
    end

data = xlsread('FullEodyData.xlsx');
t=data((514:594),1);
c=data((514:594),6);%nekroi cumul

theta0=[0.7 0.4 0.5 0.5 0.6 0.04]; %b,κ, f1, pmask,pd 0.8 0.1 0.5 0.5

options=optimset( 'MaxFunEvals',2e3);

[theta,fval,exitflag] = fminsearch(@(theta)norm(c - kineticsfmin(theta,t)),theta0,options);


save('theta')
tv=linspace(min(t),max(t));
Cfit = kineticsfmin(theta, tv);
figure(1)
plot(t, c,'r*','MarkerSize',10) %real data
hold on
plot(tv, Cfit,'g--','LineWidth',3);
legend('real data Deaths','fitted curve')
hold off
xlabel('Days (from day zero of pandemic)')
ylabel('Cumulative Deceased Individuals')
title('Delta Variant (8/15/21 - 10/15/21)')
grid on
fprintf(1,'\tRate Constants:\n')
    for theta1 = 1:length(theta)
    fprintf(1, '\t\ttheta(%d) = %8.5f\n', theta1, theta(theta1))
    end
   
%plot(t,Rsd,'r*')
end



