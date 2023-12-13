clc
clear all
close all
clc
%%ESERCIZIO1
a=6; b=0; c=2; d=3; e=7; f=6;
MAX=max([a,b,c,d,e,f]);
%parametri serbatoi
V1=4*MAX/(a+b+c+d+e+f); %m3
V2=5*MAX/(a+b+c+d+e+f); %m3
T0=25; 
Qn=V1/MAX; %m3/min
T1n=T0+20+2*MAX;
T2n=T1n+20+MAX;
P=2; %pressione costante nei serbatioi BAR
rho=1000; %kg/m3
Pvap_satur=6; %bar
lambda=2085; %KJ/kg
cp=4.186; %KJ/KG*K
Tvap=159;  %°C temperatura acqua in saturazione a sei bar
%calcolare la portatata di vapore(W1 e W2) per riscaldare i due sebatoi 
W1=Qn*rho*cp*(T1n-T0)/lambda ; %kg/min portata di vapore
W2=Qn*rho*cp*(T2n-T1n)/lambda ; %kg/min portata di vapore

Pw1=W1*lambda; %Kw
Pw2=W2*lambda; %Kw

Ti=[T1n T2n];
for i=1:7
    if i==1
        tspan=0:300;
        [tout, Tout]= ode45('Temperature1', tspan, Ti);
        figure(i)
        plot(tout, Tout(:,1),'r',tout, Tout(:,2), 'b')
        xlabel('Tempo [min]')
        ylabel('Temperatura[°C]')
        grid on
        title('Andamento temperature serbatoi')
        legend('Serbatoio 1','Serbatoio 2')
 
    elseif i==2
        tspan=0:300;
        [tout, Tout]= ode45('Temperature2', tspan, Ti);
        figure(i)
        plot(tout, Tout(:,1),'r',tout, Tout(:,2), 'b')
        xlabel('Tempo [min]')
        ylabel('Temperatura[°C]')
        grid on
        title('Andamento temperature serbatoi')
        legend('Serbatoio 1','Serbatoio 2')
    
    elseif i==3
        tspan=0:300;
        [tout, Tout]= ode45('Temperature3', tspan, Ti);
        figure(i)
        plot(tout, Tout(:,1),'r',tout, Tout(:,2), 'b')
        xlabel('Tempo [min]')
        ylabel('Temperatura[°C]')
        grid on
        title('Andamento temperature serbatoi')
        legend('Serbatoio 1','Serbatoio 2')
    
    elseif i==4
        tspan=0:300;
        [tout, Tout]= ode45('Temperature4', tspan, Ti);
        figure(i)
        plot(tout, Tout(:,1),'r',tout, Tout(:,2), 'b')
        xlabel('Tempo [min]')
        ylabel('Temperatura[°C]')
        grid on
        title('Andamento temperature serbatoi')
        legend('Serbatoio 1','Serbatoio 2')

       elseif i==5
            tspan=0:300;
            [tout, Tout]= ode45('Temperature5', tspan, Ti);
            figure(i)
            plot(tout, Tout(:,1),'r',tout, Tout(:,2), 'b')
            xlabel('Tempo [min]')
            ylabel('Temperatura[°C]')
            grid on
            title('Andamento temperature serbatoi')
            legend('Serbatoio 1','Serbatoio 2')
      elseif i==6
          tspan=0:300;
            [tout, Tout]= ode45('Temperature6', tspan, Ti);
            figure(i)
            plot(tout, Tout(:,1),'r',tout, Tout(:,2), 'b')
            xlabel('Tempo [min]')
            ylabel('Temperatura[°C]')
            grid on
            title('Andamento temperature serbatoi')
            legend('Serbatoio 1','Serbatoio 2')
    else 
        tspan=0:300;
        [tout, Tout]= ode45('Temperature7', tspan, Ti);
        figure(i)
        plot(tout, Tout(:,1),'r',tout, Tout(:,2), 'b')
        xlabel('Tempo [min]')
        ylabel('Temperatura[°C]')
        grid on
        title('Andamento temperature serbatoi')
        legend('Serbatoio 1','Serbatoio 2')
    end
end
