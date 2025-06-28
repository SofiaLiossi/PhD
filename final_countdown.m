%lsqcurve fit delta variant
function final_countdown

    function X=kinetics(theta,t)
   
initial=[4835584,4672242,80000,16000,50000,10000,50000,10000,750,23000,539];
 %initial=[4835584,4672242,20000,20000,12000,12000,10000,10000,800,53000,539];
   
%hospitalized= 600-900
                 %s u       e1 e2      i1 i2      a1    a2    h1    r     d

%opt=odeset('NonNegative',1:11,'RelTol',1e-9,'AbsTol',1e-10); %mono sthn ode45 trexei
[T,Xv]=ode23(@final,t,initial);
    function dxdt = final(t,X)       
%parameters
m=1/(82*365);
L=176; %influx rate
em=0.79;%em=0.79; %efficacy of mask
na=0.752; %mod parameter 
%na=1;
l=0.926; %sensitivy of self test
ev=0.75; %effectiveness of vaccines
v=0.0019; %(apo 14/08 ews 12/10)
psiu=1/180; %140-180
sa=1;
su=1/5.8;
gaa=1/6;
gi=1/6;
%gicu=1/10;
gh=1/18;
p=0.254;
dh=1/21; %20.4 merew noshleias
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
%theta2=k
%theta3=f1
%theta4=pm
%theta5=pd pososto anixneyshs
%theta6=f2

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
t=data((501:580),1);%1/7 - 15/10
data2=data((501:580),9);%nekroi
X=data2;

%theta1 = 0.3 + (1-0.3) .* rand(1);
%theta2 = 0.05 + (0.3-0.05) .* rand(1);
%theta3 = 0.2 + (0.8-0.2) .* rand(1);
%theta4 = 0.2 + (0.8-0.2) .* rand(1);
%theta5 = 0.2 + (0.8-0.2) .* rand(1);
%theta6 = 0.02 + (0.1-0.02) .* rand(1);

%theta0=[theta1,theta2,theta3,theta4,theta5,theta6];
theta0=[0.7,0.1,0.68,0.52,0.68,0.015];
%theta0=[0.67,0.28,0.7,0.55,0.61,0.02];

%theta0=[0.6,0.3,0.6,0.52,0.62,0.04];
%theta0=[0.49,0.17,0.6,0.75,0.72,0.02];

lb=[0,0,0,0,0,0];
ub=[1.8,1,1,1,1,1];
%lb=[];
%ub=[];
options=optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt', 'MaxFunctionEvaluations',200000, 'MaxIterations', 200000, 'StepTolerance',1e-12,'TolFun',1e-12);
[theta,Rsdnrm,residual,ExFlg,OptmInfo,Lmda,Jmat]=lsqcurvefit(@kinetics,theta0,t,X,lb,ub,options); %zeros(size(theta0))

%confidence interval
%format long g
%fprintf(1,'\tconf interval:\n')
%ci = nlparci(theta,residual,'jacobian',Jmat);
%disp(ci)
 set(gca, 'FontName', 'Arial')

 lastwarn(''); %Clear warning memory
 alpha = 0.05; %significance level
 df = length(residual) - numel(theta); %degrees of freedom or numel(theta)
 crit = tinv(1-alpha/2,df);       %critical value
 covm =pinv(Jmat'*Jmat) * var(residual); %covariance matrix
 %covm =(Jmat'.*Jmat) \ var(residual); %covariance matrix

 [~, warnId] = lastwarn; %detect if inv threw 'nearly singular' warning.
 covmIdx = sub2ind(size(covm),1:size(covm,1),1:size(covm,2));  %indices of the diag of covm
 save('covmIdx');
 save('covm(covmIdx)')
 %CI = nan(numel(theta),2);
 if ~strcmp(warnId, 'MATLAB:nearlySingularMatrix')
     CI(:,1) = theta - crit * sqrt(covm(covmIdx));
     CI(:,2) = theta + crit * sqrt(covm(covmIdx));
     disp(CI)
      sd(:,1)=covm(covmIdx)
     save('sd')
 end 
 

save('theta')
%tv=linspace(min(t),max(t));
tv=(min(t):max(t));
Cfit = kinetics(theta, tv);
save('Cfit')
sfalma=Cfit./data2;

tiledlayout(2,1)
nexttile
plot(t, X,'c*','MarkerSize',6) %real data
hold on
%param=linspace(0.2,0.6,20);
%arrayfun(@(t)plot(tv,kinetics(t,tv),'-'),param)
plot(tv, Cfit,'m','LineWidth',2); %show the best fit line


%[ypred,delta] = nlpredci(@kinetics,tv,theta,residual,'Jacobian',Jmat,'PredOpt','observation');
%lower=ypred-delta;
%upper=ypred+delta;
%save('[lower,upper]')
%plot(tv,[lower,upper],'g--','linewidth',1)

legend('real data','fitted curve','95% conf. bounds','FontSize',12)
hold off

%xlabel('Days (from day zero of pandemic)')
ylabel('Cumulative Deceased Individuals','FontSize',12)
title('Delta Variant (08/01/21 - 10/20/21)','FontSize',12)
grid on
fprintf(1,'\tRate Constants:\n')
    for theta1 = 1:length(theta)
    fprintf(1, '\t\ttheta(%d) = %8.5f\n', theta1, theta(theta1))
    end
nexttile
plot(t,sfalma,'m','linewidth',2)
grid on
xlabel('Days (from day zero of pandemic)','FontSize',12)
ylabel('fitted curve/real data','FontSize',12)
 %save('comv')
 %save('CI')
 %save('Cfit')
 figure(2)
%scatter(Cfit,residual,'r*')
scatter(t,residual,'r*')

grid on

end

