%seconda parte
clc
close all
clear all
vettoreT=load('Misure_Esercitazione_Mat_602376.mat');
ti=vettoreT.time_s-min(vettoreT.time_s); %vettore tempo che parte da zero
T_RTDi=vettoreT.T_vs_time_RTD;
T_NTCi=vettoreT.T_vs_time_NTC;
T_TCi=vettoreT.T_vs_time_TC;
tstart_RTD=14.19; %tempo morto
tstart_NTC=13.8;
tstart_TC=12.57;
tend=81; %tempo in cui assumo essere a regime
t=vettoreT.time_s-min(vettoreT.time_s); 
%il vettore tempo è stato tagliato prendendo i valori compresi tra tstart e
%tend
t_RTD=ti(ti>tstart_RTD & ti<tend);
t_NTC=ti(ti>tstart_NTC & ti<tend);
t_TC=ti(ti>tstart_TC & ti<tend);
T_RTD=vettoreT.T_vs_time_RTD(ti>tstart_RTD & ti<tend);
T_NTC=vettoreT.T_vs_time_NTC(ti>tstart_NTC & ti<tend);
T_TC=vettoreT.T_vs_time_TC(ti>tstart_TC & ti<tend);
%ho tagliato i vettori delle temperature in base a quelle che corrispondono
%ai tempi compresi tra tstart e tend cosi che non considero la parte in cui
%sono fuori dal transitorio 
figure(1) %grafici temperature sul tempo totale
hold on
plot(ti,T_RTDi,'sb')
plot(ti,T_NTCi,'or')
plot(ti,T_TCi,'vk')
xlabel('t(s)')
ylabel('T(°C)')
grid on
legend('T.RTD','T.NTC','T.TC')
hold off
figure(2) %grafici sul tempo tagliato
hold on
subplot(1,3,1)
plot(t_RTD-min(t_RTD),T_RTD,'sb')
xlabel('t(s)')
ylabel('T(°C)')
title('T.RTD')
grid on
grid minor
hold off
hold on
subplot(1,3,2)
plot(t_NTC-min(t_NTC),T_NTC,'or')
xlabel('t(s)')
ylabel('T(°C)')
title('T.NTC')
grid on
grid minor
hold off
hold on
subplot(1,3,3)
plot(t_TC-min(t_TC),T_TC,'vk')
grid on
grid minor
xlabel('t(s)')
ylabel('T(°C)')
title('T.TC')
hold off
%confronto andamento dati vs andamento del modello adottato
%RTD
figure(3)
hold on
x1=t_RTD-min(t_RTD); %vettore del tempo che inizia da 0
y1=T_RTD;
plot(x1,y1)
Tamb=num2str(y1(1));
f=fit(x1',y1',strcat('(',num2str(y1(1)),'-Tf)*exp(-a*x)+Tf'));
plot(f)
xlabel('t(s)')
ylabel('T(°C)')
grid on
legend('andamento dati','andamento teorico')
title('Termoresistenza')
Ttau_RTD=y1(1)-0.9*(y1(1)-y1(end)); %temperatura a cui faccio il 90%del salto
tau_RTD=max(x1(y1>Ttau_RTD)); %tempo corrispondente a quella temperatura 
hold off
% %NTC
figure(4)
hold on
x2=t_NTC-min(t_NTC);
y2=T_NTC;
plot(x2,y2)
Tamb=num2str(y2(1))
f2=fit(x2',y2',strcat('(',num2str(y2(1)),'-Tf)*exp(-a*x)+Tf'));
plot(f2)
Ttau_NTC=y2(1)-0.9*(y2(1)-y2(end));
tau_NTC=max(x2(y2>Ttau_NTC)); 
xlabel('t(s)')
ylabel('T(°C)')
grid on
legend('andamento dati','andamento teorico')
title('Termistore')
hold off
%TC
figure(5)
hold on
x3=t_TC-min(t_TC);
y3=T_TC;
plot(x3,y3)
Tamb=num2str(y3(1));
f3=fit(x3',y3',strcat('(',num2str(y3(1)),'-Tf)*exp(-a*x)+Tf'));
plot(f3)
xlabel('t(s)')
ylabel('T(°C)')
grid on
legend('andamento dati','andamento teorico')
title('Termocoppia')
Ttau_TC=y3(1)-0.9*(y3(1)-y3(end));
tau_TC=max(x3(y3>Ttau_TC));
tau_tot=[tau_RTD,tau_NTC,tau_TC]
hold off
%temperature medie finali
TfinRTD=mean(T_RTD(t_RTD>20));
TfinNTC=mean(T_NTC(t_NTC>20));
TfinTC=mean(T_TC(t_TC>20));
Tfin=[TfinRTD TfinNTC TfinTC]