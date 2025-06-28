%ga method final 

%lsqcurve fit
function final_ga_method
    function C=kinetics_ga(theta,t)
        
%function value=residual norm 
%initial=[4835584,4672242,90000,10000,90000,10000,90000,10000,5400,100000,7906];%12000
%initial=[2918000,7400000,84000,30000,84000,30000,84000,30000,2800,100000,7906];%12000
initial=[3318000,6400000,330000,170000,160000,85000,160000,85000,2800,50000,9500];%12000

%19/8 exoume energa=15000 kroysmata ara E=30000
%14/8 exoume 20000=Id
%20/9 12600

opt=odeset('NonNegative',1:11);
%opt = odeset('RelTol',1e-9,'AbsTol',1e-10);
[T,Cv]=ode23(@DifEq1,t,initial,opt);

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
%f2=4/100;
%ff2=4/1000;
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
%t=data((650:749),1); %494-604
%c=data((650:749),9);%nekroi cumul

t=data((670:730),1); %494-604 δελτα
c=data((670:730),9);%nekroi cumul δελτα


ftns = @(theta) norm(c-kinetics_ga(theta,t));
PopSz = 500;
Parms = 6;
%opts = optimoptions('ga', 'InitialPopulationMatrix',[randi([20 80], 50, 5)
%randi([95 110], 50,
%6)*1E-2]);'InitialPopulationMatrix',randi(1E+4,PopSz,Parms)*1E-3,or randi[0
%2],PopSz,Parms
opts = optimoptions('ga', 'PopulationSize',PopSz, 'MaxGenerations',500, 'PlotFcn','gaplotbestf');
t0 = clock;
fprintf('\nStart Time: %4d-%02d-%02d %02d:%02d:%07.4f\n', t0)
[theta,fval,exitflag,output] = ga(ftns, Parms, [],[],[],[],[0,0,0.3,0.3,0.3, 0],[1.5,1,0.8,1,1,0.1],[],[],opts)
%[theta,fval,exitflag,output] = ga(ftns, Parms, [],[],[],[],[],[],[],[],opts)

t1 = clock;
fprintf('\nStop Time: %4d-%02d-%02d %02d:%02d:%07.4f\n', t1)
GA_Time = etime(t1,t0)

fprintf(1,'\tRate Constants:\n')
for k1 = 1:length(theta)
    fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1, theta(k1))
end
save('theta');
tv = linspace(min(t), max(t));
Cfit = kinetics_ga(theta, tv);

figure(1)
hd = plot(t, c, 'pr');
%legend(hd,'')
hold on
hlp = plot(tv, Cfit,'b');
hold off
grid
xlabel('Time (days)')
ylabel('Deceased Individuals')
legend('real data','fitted curve')
end