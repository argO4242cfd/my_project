clc
clear all
close all
vettoreT=load('Misure_Esercitazione_Mat_602376.mat');
T0=vettoreT.DatiRTD_T0;
T100=vettoreT.DatiRTD_T100;
T=[T0 , T100];
%R0
R0=vettoreT.DatiRTD_R0; %valore di resistenza associata al ghiaccio
uA_R0=std(R0)/sqrt(length(R0)); %incertezza categoria A
a_R0=0.1; %accuracy
uB_R0=a_R0/sqrt(3); %incertezza di categoria B
uC_R0=sqrt(uA_R0^2+uB_R0^2);%incertezza composta
U_R0=1.95*uC_R0; %incertezza estesa con fattore di copertura k=1.95
%R100
R100=vettoreT.DatiRTD_R100; %vettore delle resistenza a 100 gradi
uA_R100=std(R100)/sqrt(length(R100));
a_R100=0.1;
uB_R100=a_R100/sqrt(3);
uC_R100=sqrt(uA_R100^2+uB_R100^2);
U_R100=1.95*uC_R100;
%alpha
alfa=(mean(R100)-mean(R0))/(100);%coeff angolare retta di taratura
Ualfa=(1.95/100)*(uC_R0^2+uC_R100^2)^0.5 %errore propagato dalle due incertezze composte
%curva di taratura
R=[R0,R100];
figure(1)
hold on
plot(T,R,'b')
grid on
xlim([0,100])
xlabel('T(°C)')
ylabel('R(Ω)')
T1=[0,100];
R1=mean(R0)*(1+0.0039083*T1-5.77500*10^-7*T1.^2);
plot(T1,R1,'r')
legend('curva di taratura ottenuta','modello dataheet')
box on
hold off
T0media=mean(T0);
T100media=mean(T100);
R0media=mean(R0);
R100media=mean(R100);
%%metodo montecarlo
%distribuzioni di probabilità uniformi da cui sono estratti i dati di R0 e
%R100
pd1=makedist('uniform','lower',mean(R0)-uC_R0,'upper',mean(R0)+uC_R0);
pd2=makedist('uniform','lower',mean(R100)-uC_R100,'upper',mean(R100)+uC_R100);
R0mc=random(pd1,1,1e6);  %dati estratti
R100mc=random(pd2,1,1e6);
alfa_mc=(R100mc-R0mc)/(100); 
figure(2)
histogram(alfa_mc) 
xlabel('coeff angolare')
ylabel('Numero di estrazioni')
grid on 
figure(3)
histogram(alfa_mc,'normalization','cdf') %istogramma secondo la funzione di distribuzione cumulata
hold on
grid on
xlabel('coeff angolare')
ylabel('cdf')
xlim([min(alfa_mc) max(alfa_mc)])
dvst=std(alfa_mc); %dev standard del coefficiente angolare
alfamedio=mean(alfa_mc);
Umc=(0.3843-0.3825)/2 %incertezza estesa 
%0.3834 e 0.3825 presi vedendo l'istogramma della cdf e in prossimità di
%cdf=0.025 e cdf=0.975
k=Umc/dvst %fattore di copertura 
