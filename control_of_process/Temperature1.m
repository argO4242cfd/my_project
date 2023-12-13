function [dTdt]=Temperature1(t, T)
%costanti
a=6; b=0; c=2; d=3; e=7; f=6;
MAX=max([a,b,c,d,e,f]);
%parametri serbatoi
V1=4*MAX/(a+b+c+d+e+f); %m3
V2=5*MAX/(a+b+c+d+e+f); %m3
T0=25; 
Qn=V1/MAX; %m3/min

rho=1000; %kg/m3

T1n=T0+20+2*MAX;
T2n=T1n+20+MAX;
lambda=2085; %KJ/kg
cp=4.186; %KJ/KG*K
  %°C temperatura acqua in saturazione a sei bar

%variabili T

% potenza fornita vapore

W1=Qn*rho*cp*(T1n-T0)/lambda ; %kg/min portata di vapore
W2=Qn*rho*cp*(T2n-T1n)/lambda ; %kg/min portata di vapore

Pw1=W1*lambda; %Kw
Pw2=W2*lambda; 

Tv1=T(1);
Tv2=T(2);
if t>=100
          Qn=2*Qn;
end

[dTdt]=[(rho*cp*Qn*(T0-Tv1)+Pw1)/(rho*cp*V1); (rho*cp*Qn*(Tv1-Tv2)+Pw2)/(rho*cp*V2)];
end

   




