%lsqcurve fit
function new_model_multistart 

function X=kineticsmulti(theta,t)  

%ιnitial=[4035584,5672242,25000,25000,15000,5000,50000,30000,4622,1111358,13206];
%initial=[4035584,5672242,30000,15000,17000,10000,15000,15000,4622,1111358,13206];
%initial=[4035584,5672242,35000,20000,15000,5000,15000,5000,4622,1111358,13206];
initial=[4035584,5672242,50000,40000,30000,30000,30000,30000,4300,1111358,24];


opt=odeset('NonNegative',1:11,'RelTol',1e-8,'AbsTol',1e-10); %mono sthn ode45 trexei
[T,Xv]=ode15s(@final,t,initial,opt);
    function dxdt = final(t,X)
%parameters
m=1/(82*365);
L=176; %influx rate
em=0.79;%em=0.79; %efficacy of mask
na=0.752; %mod parameter 
l=0.926; %sensitivy of self test
ev=0.75; %effectiveness of vaccines
v=0.0019; %(apo 14/08 ews 12/10)
psiu=1/180;
sa=1;
su=1/5.8;
f2=4/100;
ff2=4/1000;
%f3=10.8/100;
%f4=23.8/100;
%f5=91.7/100;
gaa=1/6;
gi=1/6;
%gicu=1/10;
gh=1/18;
p=0.254;
dh=1/18; %20.4 merew noshleias
%dicu=1/10;
sh=1/4;
a=0.17;
k=0;

%variables
S=X(1);
U=X(2);
E1=X(3);
E2=X(4);
E=X(3)+X(4);
I1=X(5);
I2=X(6);
I=X(5)+X(6);
A1=X(7);
A2=X(8);
A=X(7)+X(8);
H1=X(9);
R=X(10);
D=X(11);

%theta1=b
%theta2=f1 pososto sympt
%theta3=pmask 
%theta4=pd
%theta5=k pososto anixneyshs

N=S+U+E+I+A+H1+R;
  ls=(theta(1).*(1-theta(4)).*(1-em.*theta(3)).*(E+(1-theta(2)).*na.*A+p*I))./((N-((1-p)*I+theta(2).*A+H1+ev.*U)));
  lu= (1-ev).*ls;
%disp(ls);

%differential equations
dxdt = zeros (11, 1);
dxdt(1)= L-ls.*S+psiu.*U-v.*S-m.*S; %s
dxdt(2)=-lu.*U-psiu.*U+v.*S-m.*U; %u
dxdt(3)=ls.*S-su.*E1-m.*E1; %e1
dxdt(4)=lu.*U-su.*E2-m.*E2; %e2

dxdt(5)=theta(3).*su.*E1-theta(6).*sh.*I1-(1-theta(6)).*gi.*I1+theta(2).*l.*sa.*A1-m.*I1;%i1
dxdt(6)=theta(3).*su.*E2-(theta(6)/10).*sh.*I2-(1-theta(6)/10).*gi.*I2+theta(2).*l.*sa.*A2-m.*I2;%i2

dxdt(7)=(1-theta(3)).*su.*E1-gaa.*A1-theta(2).*l.*sa.*A1-m.*A1; %A1
dxdt(8)=(1-theta(3)).*su.*E2-gaa.*A2-theta(2).*l.*sa.*A2-m.*A2;

dxdt(9)=theta(6).*sh.*I1+(theta(6)/10).*sh.*I2-a.*dh.*H1-(1-a).*gh.*H1-m.*H1; %h1

%dxdt(11)=(1-f3).*H1-f5.*gh.*Hnonicu-(1-f5).*dh.*Hnonicu-m.*Hnonicu;%hnonicu
%dxdt(11)=dxdt(9)-dxdt(10);%hnonicu

dxdt(10)=(1-theta(6)).*gi.*I1+(1-(theta(6))).*gi.*I2+gaa.*A+(1-a).*gh.*H1-m.*R; %recovered

dxdt(11)=a.*dh.*H1;% %deads
    end
X=Xv(:,11)
end

%  data for FITTING
data = xlsread('FullEodyData.xlsx');
t=data((650:749),1);
data2=data((650:749),9);%nekroi
X=data2;

theta0=[0.5,0.05,0.5,0.4,0.4,0.04];

%multistart optimisation%
problem = createOptimProblem ('lsqcurvefit',...
                              'objective', @kineticsmulti, ...
                              'x0',theta0, ...
                                'lb', [0,0,0.2,0,0,0], 'ub',[1.5,1,0.8,1,1,0.1],...
                                'xdata',t,'ydata',X)
                            %'options',optimset('OutputFcn',@PlotIterates));

ms = MultiStart('PlotFcns',@gsplotbestf);

%tic,[theta,fval,exitflag,output,solutions]=run(ms,problem,50);toc
[theta,fval,exitflag,output,solutions]=run(ms,problem,200);
save('theta')

%plot(t,X,'r');
%theta=lsqcurvefit(problem);
end
