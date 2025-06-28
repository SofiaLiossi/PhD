
% grafima new cases - deads - variants
t=FullEodyData{(1:803),1};
deads=FullEodyData{(1:803),3};
cases= FullEodyData{(1:803),2};

yyaxis left
plot(t,cases,'LineWidth',1,'Color',[0.4660, 0.6740, 0.1880])
ylabel('new cases','FontSize',12)
hold on
x1 = datetime(2020,1,1);
x2 = datetime(2021,1,17);
x3=datetime(2021,6,20);
x4=datetime(2021,12,12);
x5=datetime(2022,5,31);
x = [x1 x2 x3 x4;x1 x2 x3 x4; x2 x3 x4 x5;x2 x3 x4 x5];
y=[0 0 0 0;6*10^4 6*10^4 6*10^4 6*10^4;6*10^4 6*10^4 6*10^4 6*10^4;0 0 0 0];
c=[4 3 2 1];
fill(x,y,c,'FaceAlpha',0.1)

yyaxis right
plot(t,deads,'LineWidth',1,'Color',[0.4940, 0.1840, 0.5560])
ylabel('deaths')
xlabel('days')
title('COVID-19 Evolution in Greece')
annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', "Others")
annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', "Alpha")
annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', "Delta")
annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', "Omicron")

hold off

