function sensitivity 
syms  su m f1 f2 gi sh ff2 k l sa gaa a dh gh A B p na b pd em pm ev

A=b*(1-pd).*(1-em.*pm);
B=(1-ev)*A;

R=(A.*m.^3 + A.*B.*m.^3 + A.*gaa.*gi.^2 + A.*gaa.*m.^2 + 2.*A.*gi.*m.^2 + A.*gi.^2.*m + A.*B.*gaa.*gi.^2 + A.*B.*gaa.*m.^2 + 2.*A.*B.*gi.*m.^2 + A.*B.*gi.^2.*m - A.*f2.*gaa.*gi.^2 - A.*ff2.*gaa.*gi.^2 - A.*f2.*gi.*m.^2 - A.*f2.*gi.^2.*m - A.*ff2.*gi.*m.^2 - A.*ff2.*gi.^2.*m + A.*f2.*m.^2.*sh + A.*ff2.*m.^2.*sh + A.*gi.^2.*na.*su + A.*m.^2.*na.*su + 2.*A.*gaa.*gi.*m + 2.*A.*B.*gaa.*gi.*m - A.*f2.*gaa.*gi.*m - A.*ff2.*gaa.*gi.*m + A.*f2.*gaa.*gi.*sh + A.*ff2.*gaa.*gi.*sh + A.*f2.*gaa.*m.*sh + A.*f2.*gi.*m.*sh + A.*ff2.*gaa.*m.*sh + A*ff2*gi*m*sh + 2*A*gi*m*na*su - A*B*f2*gaa*gi^2 - A*B*ff2*gaa*gi^2 - A*B*f2*gi*m^2 - A*B*f2*gi^2*m - A*B*ff2*gi*m^2 - A*B*ff2*gi^2*m + A*B*f2*m^2*sh + A*B*ff2*m^2*sh + A*B*gi^2*na*su + A*B*m^2*na*su + A*f2*ff2*gaa*gi^2 + A*f2*ff2*gi^2*m + A*f2*ff2*gaa*sh^2 + A*f2*ff2*m*sh^2 - A*f1*gi^2*na*su - A*f2*gi^2*na*su - A*ff2*gi^2*na*su + A*gi^2*k*l*sa - A*gi^2*k*na*su - A*f1*m^2*na*su + A*f1*m^2*p*su + A*k*l*m^2*sa - A*k*m^2*na*su + A*f1*f2*gi^2*na*su + A*f1*ff2*gi^2*na*su + A*f2*ff2*gi^2*na*su - A*f2*gi^2*k*l*sa - A*ff2*gi^2*k*l*sa + A*f1*gi^2*k*na*su + A*f2*gi^2*k*na*su + A*ff2*gi^2*k*na*su + A*f1*k*m^2*na*su + A*f2*ff2*na*sh^2*su + A*B*f2^2*gi^2*na*su - A*B*f2*gaa*gi*m - A*B*ff2*gaa*gi*m + A*B*f2*gaa*gi*sh + A*B*ff2*gaa*gi*sh + A*B*f2*gaa*m*sh + A*B*f2*gi*m*sh + A*B*ff2*gaa*m*sh + A*B*ff2*gi*m*sh + 2*A*B*gi*m*na*su - 2*A*f2*ff2*gaa*gi*sh - 2*A*f2*ff2*gi*m*sh + A*f1*gaa*gi*p*su - 2*A*f1*gi*m*na*su - A*f2*gi*m*na*su - A*ff2*gi*m*na*su + A*f1*gaa*m*p*su + A*f1*gi*m*p*su + 2*A*gi*k*l*m*sa - 2*A*gi*k*m*na*su + A*f2*gi*na*sh*su + A*ff2*gi*na*sh*su + A*f2*m*na*sh*su + A*ff2*m*na*sh*su + A*B*f2*ff2*gaa*gi^2 + A*B*f2*ff2*gi^2*m + A*B*f2*ff2*gaa*sh^2 + A*B*f2*ff2*m*sh^2 - 2*A*B*f2*gi^2*na*su - A*B*ff2*gi^2*na*su + A*B*gi^2*k*l*sa - A*B*gi^2*k*na*su - A*B*f2*m^2*na*su + A*B*f1*m^2*p*su + A*B*k*l*m^2*sa - A*B*k*m^2*na*su + A*B*f1*gaa*m*p*su + A*B*f1*gi*m*p*su + 2*A*B*gi*k*l*m*sa - 2*A*B*gi*k*m*na*su + A*B*f2*gi*na*sh*su + A*B*ff2*gi*na*sh*su + A*B*f2*m*na*sh*su + A*B*ff2*m*na*sh*su - A*f1*ff2*gaa*gi*p*su + A*f1*f2*gi*m*na*su + A*f1*ff2*gi*m*na*su - A*f1*ff2*gi*m*p*su - A*f2*gi*k*l*m*sa - A*ff2*gi*k*l*m*sa + 2*A*f1*gi*k*m*na*su + A*f2*gi*k*m*na*su + A*ff2*gi*k*m*na*su - A*f1*f2*gi*na*sh*su - A*f1*ff2*gi*na*sh*su - 2*A*f2*ff2*gi*na*sh*su + A*f1*ff2*gaa*p*sh*su + A*f2*gi*k*l*sa*sh + A*ff2*gi*k*l*sa*sh - A*f2*gi*k*na*sh*su - A*ff2*gi*k*na*sh*su - A*f1*f2*m*na*sh*su - A*f1*ff2*m*na*sh*su + A*f1*ff2*m*p*sh*su + A*f2*k*l*m*sa*sh + A*ff2*k*l*m*sa*sh - A*f2*k*m*na*sh*su - A*ff2*k*m*na*sh*su + A*gi*k*l*p*sa*su + A*k*l*m*p*sa*su + 2*A*B*f2*ff2*gi^2*na*su - A*B*f2*gi^2*k*l*sa - A*B*ff2*gi^2*k*l*sa + 2*A*B*f2*gi^2*k*na*su + A*B*ff2*gi^2*k*na*su + A*B*f2^2*gi*m*na*su + A*B*f2*k*m^2*na*su + A*B*f2*ff2*na*sh^2*su - A*B*f2^2*gi*na*sh*su - A*B*f2^2*m*na*sh*su - A*f1*f2*ff2*gi^2*na*su + A*f2*ff2*gi^2*k*l*sa - A*f1*f2*gi^2*k*na*su - A*f1*ff2*gi^2*k*na*su - A*f2*ff2*gi^2*k*na*su - A*f1*f2*ff2*na*sh^2*su + A*f2*ff2*k*l*sa*sh^2 - A*f2*ff2*k*na*sh^2*su - A*B*f2^2*ff2*gi^2*na*su - A*B*f2^2*gi^2*k*na*su - A*B*f2^2*ff2*na*sh^2*su - 2*A*B*f2*ff2*gaa*gi*sh - 2*A*B*f2*ff2*gi*m*sh + A*B*f1*gaa*gi*p*su - 3*A*B*f2*gi*m*na*su - A*B*ff2*gi*m*na*su + A*B*f2^2*ff2*gi^2*k*na*su + A*B*f2^2*ff2*k*na*sh^2*su - A*B*f1*f2*gaa*gi*p*su + A*B*f2*ff2*gi*m*na*su - A*B*f1*f2*gi*m*p*su - A*B*f2*gi*k*l*m*sa - A*B*ff2*gi*k*l*m*sa + 3*A*B*f2*gi*k*m*na*su + A*B*ff2*gi*k*m*na*su - 3*A*B*f2*ff2*gi*na*sh*su + A*B*f1*f2*gaa*p*sh*su + A*B*f2*gi*k*l*sa*sh + A*B*ff2*gi*k*l*sa*sh - A*B*f2*gi*k*na*sh*su - A*B*ff2*gi*k*na*sh*su - A*B*f2*ff2*m*na*sh*su + A*B*f1*f2*m*p*sh*su + A*B*f2*k*l*m*sa*sh + A*B*ff2*k*l*m*sa*sh - A*B*f2*k*m*na*sh*su - A*B*ff2*k*m*na*sh*su + A*B*gi*k*l*p*sa*su + A*B*k*l*m*p*sa*su - A*f1*f2*gi*k*m*na*su - A*f1*ff2*gi*k*m*na*su + 2*A*f1*f2*ff2*gi*na*sh*su - 2*A*f2*ff2*gi*k*l*sa*sh + A*f1*f2*gi*k*na*sh*su + A*f1*ff2*gi*k*na*sh*su + 2*A*f2*ff2*gi*k*na*sh*su + A*f1*f2*k*m*na*sh*su + A*f1*ff2*k*m*na*sh*su - A*ff2*gi*k*l*p*sa*su + A*ff2*k*l*p*sa*sh*su + A*B*f2*ff2*gi^2*k*l*sa - 2*A*B*f2*ff2*gi^2*k*na*su - A*B*f2^2*gi*k*m*na*su + 2*A*B*f2^2*ff2*gi*na*sh*su + A*B*f2*ff2*k*l*sa*sh^2 - A*B*f2*ff2*k*na*sh^2*su + A*B*f2^2*gi*k*na*sh*su + A*B*f2^2*k*m*na*sh*su + A*f1*f2*ff2*gi^2*k*na*su + A*f1*f2*ff2*k*na*sh^2*su - 2*A*B*f2^2*ff2*gi*k*na*sh*su + A*B*f2^2*gi*k*l*p*sa*su - A*B*f2^2*k*l*p*sa*sh*su - A*B*f2*ff2*gi*k*m*na*su - 2*A*B*f2*ff2*gi*k*l*sa*sh + 3*A*B*f2*ff2*gi*k*na*sh*su + A*B*f2*ff2*k*m*na*sh*su + A*B*f1*gi*k*l*p*sa*su - 2*A*B*f2*gi*k*l*p*sa*su + A*B*f1*k*l*m*p*sa*su - A*B*f2*k*l*m*p*sa*su + A*B*f2*k*l*p*sa*sh*su - 2*A*f1*f2*ff2*gi*k*na*sh*su - A*B*f1*f2*gi*k*l*p*sa*su + A*B*f1*f2*k*l*p*sa*sh*su)/(gaa*m^3 + 2*gi*m^3 + m^3*su + m^4 + gi^2*m^2 + gaa*gi^2*su + f2*m^3*sh + ff2*m^3*sh + gaa*m^2*su + 2*gi*m^2*su + gi^2*m*su - f2*gi^2*m^2 - ff2*gi^2*m^2 - f2*gi*m^3 - ff2*gi*m^3 + 2*gaa*gi*m^2 + gaa*gi^2*m + 2*gaa*gi*m*su - f2*gaa*gi*m^2 - f2*gaa*gi^2*m - ff2*gaa*gi*m^2 - ff2*gaa*gi^2*m - f2*gaa*gi^2*su - ff2*gaa*gi^2*su + f2*gaa*m^2*sh + f2*gi*m^2*sh - f2*gi*m^2*su - f2*gi^2*m*su + ff2*gaa*m^2*sh + ff2*gi*m^2*sh - ff2*gi*m^2*su - ff2*gi^2*m*su + k*l*m^3*sa + f2*m^2*sh*su + ff2*m^2*sh*su + f2*ff2*gi^2*m^2 + f2*ff2*m^2*sh^2 + f2*gaa*gi*m*sh - f2*gaa*gi*m*su + ff2*gaa*gi*m*sh - ff2*gaa*gi*m*su + f2*gaa*gi*sh*su + ff2*gaa*gi*sh*su + f2*gaa*m*sh*su + f2*gi*m*sh*su + ff2*gaa*m*sh*su + ff2*gi*m*sh*su + f2*ff2*gaa*gi^2*m + f2*ff2*gaa*gi^2*su + f2*ff2*gaa*m*sh^2 - 2*f2*ff2*gi*m^2*sh + f2*ff2*gi^2*m*su + f2*ff2*gaa*sh^2*su + 2*gi*k*l*m^2*sa + gi^2*k*l*m*sa + f2*ff2*m*sh^2*su + gi^2*k*l*sa*su + k*l*m^2*sa*su - f2*gi*k*l*m^2*sa - f2*gi^2*k*l*m*sa - ff2*gi*k*l*m^2*sa - ff2*gi^2*k*l*m*sa - f2*gi^2*k*l*sa*su - ff2*gi^2*k*l*sa*su + f2*k*l*m^2*sa*sh + ff2*k*l*m^2*sa*sh - 2*f2*ff2*gaa*gi*m*sh - 2*f2*ff2*gaa*gi*sh*su - 2*f2*ff2*gi*m*sh*su + 2*gi*k*l*m*sa*su + f2*gi*k*l*m*sa*sh - f2*gi*k*l*m*sa*su + ff2*gi*k*l*m*sa*sh - ff2*gi*k*l*m*sa*su + f2*gi*k*l*sa*sh*su + ff2*gi*k*l*sa*sh*su + f2*k*l*m*sa*sh*su + ff2*k*l*m*sa*sh*su + f2*ff2*gi^2*k*l*m*sa + f2*ff2*gi^2*k*l*sa*su + f2*ff2*k*l*m*sa*sh^2 + f2*ff2*k*l*sa*sh^2*su - 2*f2*ff2*gi*k*l*m*sa*sh - 2*f2*ff2*gi*k*l*sa*sh*su);

sensx=diff(R,ev)*(ev/R)

%omicron
%b=0.85889;
%k=0.69675;
%f1=0.19211;
%pm=0.63441;
%pd=0.64568;
%f2=0.02177;
%delta
b=0.49735;
k=0.1726;
f1=0.36467;
pm=0.45146;
pd=0.62221;
f2=0.0375;

m=1/(82*365);
L=176; %influx rate
em=0.79; 
na=0.752; %mod parameter 
l=0.926; %sensitivy of self test in delta
%l=0.78; %sensitivity in omicron
ev=0.75; %effectiveness of vaccines in delta
%ev= 0.7; %efectivness in omicron
v=0.0019; %(apo 14/08 ews 12/10)
psiu=1/180;
sa=1;
su=1/5.8;
gaa=1/6;
gi=1/6;
gh=1/18;
p=0.254;
dh=1/18; %20.4 merew noshleias
sh=1/4;
a=0.17;
ff2=f2/10;
subs(vpa(sensx))
%subs(vpa(R))
end   
