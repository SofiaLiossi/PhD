%lsqcurve fit delta variant
function final_countdown_omicron

    function X=kinetics(theta,t)
   
%initial=[4035584,4672242,35000,12000,32000,8000,32000,8000,4200,1111358,295];
initial=[3318000,6400000,186000,94000,100000,50000,100000,50000,2100,27000,7247];%12000
%active cases =100000
%hospitalized=300 atoma *6 = 1800 problepsi=1600-2800
%recovered= 56000*2=100000
%3months 416117 cases + 7000000 vaccination
               %s u          e1 e2      i1 i2      a1    a2    h1    r     d

%opt=odeset('NonNegative',1:11,'RelTol',1e-9,'AbsTol',1e-10); %mono sthn ode45 trexei
[T,Xv]=ode23(@final,t,initial);
    function dxdt = final(t,X)       
%parameters
m=1/(82*365);
L=176; %influx rate
em=0.79;%em=0.79; %efficacy of mask
na=0.752; %mod parameter 
l=0.78; %sensitivy of self test in omicron
ev=0.7; %effectiveness of vaccines in omicron
v=0.0019; %(apo 14/08 ews 12/10)
psiu=1/180;%140-180
sa=1;
su=1/5.8;

%f3=10.8/100;
%f4=23.8/100;
%f5=91.7/100;
gaa=1/6;
gi=1/6;
%gicu=1/10;
gh=1/18;
p=0.254;
dh=1/21; %20.4 merew noshleias 21
%dicu=1/10;
sh=1/4;
a=0.17;


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
  ls=(theta(1).*(1-theta(5)).*(1-em.*theta(4)).*(E+(1-theta(2)).*na.*A+p.*I))./((N-((1-p).*I+theta(2).*l.*A+H1+ev.*U)));
  lu= (1-ev).*ls;
%disp(ls);

%differential equations
dxdt = zeros (11, 1);
dxdt(1)= L-ls.*S+psiu.*U-v.*S-m.*S; %s
dxdt(2)=-lu.*U-psiu.*U+v.*S-m.*U; %u
dxdt(3)=ls.*S-su.*E1-m.*E1; %e1
dxdt(4)=lu.*U-su.*E2-m.*E2; %e2

dxdt(5)=theta(3).*su.*E1-theta(6).*sh.*I1-(1-theta(6)).*gi.*I1+theta(2).*l.*sa.*A1-m.*I1;%i1
dxdt(6)=theta(3).*su.*E2-(theta(6)./10).*sh.*I2-(1-(theta(6)./10)).*gi.*I2+theta(2).*l.*sa.*A2-m.*I2;%i2

dxdt(7)=(1-theta(3)).*su.*E1-gaa.*A1-theta(2).*l.*sa.*A1-m.*A1; %A1
dxdt(8)=(1-theta(3)).*su.*E2-gaa.*A2-theta(2).*l.*sa.*A2-m.*A2;

dxdt(9)=theta(6).*sh.*I1+(theta(6)./10).*sh.*I2-a.*dh.*H1-(1-a).*gh.*H1-m.*H1; %h1

dxdt(10)=(1-theta(6)).*gi.*I1+(1-(theta(6)./10)).*gi.*I2+gaa.*A+(1-a).*gh.*H1-m.*R; %recovered

dxdt(11)=a.*dh.*H1;% %deads


    end
X=Xv(:,11);
%X=Xv(:,[11,13]);
%X = cumsum(Xv(:,11));
    end
data = xlsread('FullEodyData.xlsx');
t=data((640:730),1);%551-690
data2=data((640:730),9);%nekroi
X=data2;

theta0=[0.8,0.3,0.8,0.7,0.5,0.01];
lb=[0,0,0,0,0,0];
ub=[1.9,1,1,1,1,1];

options=optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt', 'MaxFunctionEvaluations',200000, 'MaxIterations', 200000, 'StepTolerance',1e-12,'TolFun',1e-12);
[theta,Rsdnrm,residual,ExFlg,OptmInfo,Lmda,Jmat]=lsqcurvefit(@kinetics,theta0,t,X,lb,ub,options); %zeros(size(theta0))

save('theta')
%tv=linspace(min(t),max(t));
tv=(min(t):1:max(t));
Cfit = kinetics(theta, tv);
tiledlayout(2,1)
nexttile
plot(t, X,'c*','MarkerSize',6) %real data
hold on
plot(tv, Cfit,'m','LineWidth',2);

%plot conf interval
%[ypred,delta] = nlpredci(@kinetics,tv,theta,residual,'Jacobian',Jmat,'PredOpt','observation');
%lower=ypred-delta;
%upper=ypred+delta;
%plot(tv,[lower,upper],'g--','linewidth',2)
%save('[lower,upper]')
legend('real data','fitted curve','95% conf. bounds','FontSize',12)

hold off
%xlabel('Days (from day zero of pandemic)')
ylabel('Cumulative Deceased Individuals','FontSize',12)
title('Omicron Variant (12/19/21 - 03/24/22)','FontSize',12)
grid on
fprintf(1,'\tRate Constants:\n')
    for theta1 = 1:length(theta)
    fprintf(1, '\t\ttheta(%d) = %8.5f\n', theta1, theta(theta1))
    end
 
%confidence interval
%format long g
%fprintf(1,'\tconf interval:\n')
%ci = nlparci(theta,residual,'jacobian',Jmat);
%disp(ci)

 lastwarn(''); %Clear warning memory
 alpha = 0.05; %significance level
 df = length(residual) - numel(theta); %degrees of freedom or numel(theta)
 crit = tinv(1-alpha/2,df);       %critical value
 covm =pinv(Jmat'*Jmat) * var(residual); %covariance matrix
 %covm =(Jmat'.*Jmat) \ var(residual); %covariance matrix

 [~, warnId] = lastwarn; %detect if inv threw 'nearly singular' warning.
 covmIdx = sub2ind(size(covm),1:size(covm,1),1:size(covm,2));  %indices of the diag of covm
 %CI = nan(numel(theta),2);
 
 if ~strcmp(warnId, 'MATLAB:nearlySingularMatrix')
     CI(:,1) = theta - crit * sqrt(covm(covmIdx));
     CI(:,2) = theta + crit * sqrt(covm(covmIdx));
     sd(:,1)=covm(covmIdx)
     save('sd')
     save('CI')
     disp(CI)
 end  
nexttile
 sfalma=Cfit./data2;
 figure(1)
plot(t,sfalma,'m','LineWidth',2)
xlabel('Days (from day zero of pandemic)','FontSize',12)
ylabel('fitted curve/real data','FontSize',12)
grid on
hold off
%save('comv')
 %save('CI')
 %save('Cfit')
 figure(2)
stem(t,residual,'r*')
grid on
end

