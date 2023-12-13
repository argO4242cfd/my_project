function [dTdt]=Temperature3(t,T)
        a=6; b=0; c=2; d=3; e=7; f=6;
        
        cp=4.186; % [Kj/Kg*K]
       
        lambda=2085; % [Kj/Kg]
        rho=1000; % [Kg/m^3]
        V1=4*max([a,b,c,d,e,f])/(a+b+c+d+e+f); % [m^3]
        V2=5*max([a,b,c,d,e,f])/(a+b+c+d+e+f); % [m^3]
        Qn=V1/max([a,b,c,d,e,f]); % Portata nominale di acqua [m^3/min]
        T0=25; % [°C]
        T1n=T0+20+2*max([a,b,c,d,e,f]); % Temperatura nominale [°C]
        T2n=T1n+20+max([a,b,c,d,e,f]); % Temperatura nominale [°C]
        W1n=Qn*rho*cp*(T1n-T0)/lambda; % Portata nominale di vapore 1 [Kg/min]
        W2n=Qn*rho*cp*(T2n-T1n)/lambda; % Portata nominale di vapore 2 [Kg/min]

        W1=W1n;
        W2=W2n;
        T1=T(1);
        T2=T(2);
        A=0.1*Qn;
        Periodo=(V1+V2)/Qn;
        omega=2*3.14/Periodo/4;
        Qn0=Qn;

        % Ingresso sinusoidale per t=10
          if t>=10
           Qn=Qn0+A*sin(omega*(t-10));
          end
        
        [dTdt]=[(rho*Qn*cp*(T0-T1)+W1*lambda)/(rho*cp*V1); (rho*Qn*cp*(T1-T2)+W2*lambda)/(rho*cp*V2)];
end