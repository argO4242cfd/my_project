clc
clear all
close all
vettoreT=load('Misure_Esercitazione_Mat_602376.mat');
%accuracy
a_RTD=0.5; a_NTC=2; a_TC=2; 
a_tot=[ a_RTD , a_NTC , a_TC];
T_RTD=vettoreT.DatiRTD_Tfixed_100deg;
T_NTC=vettoreT.DatiNTC_Tfixed_100deg;
T_TC=vettoreT.DatiTC_Tfixed_100deg;
T_tot=[T_RTD ; T_NTC ; T_TC];
%valori medi e incertezze di tipo A,composta,estesa
Tm_RTD=mean(T_RTD); Tm_NTC=mean(T_NTC); Tm_TC=mean(T_TC);
dvst_RTD=std(T_RTD); dvst_NTC=std(T_NTC); dvst_TC=std(T_TC);
dvst=[dvst_RTD , dvst_NTC , dvst_TC];
for i=1:3 %u=(u_RTD u_NTC u_TC)
    uA(i)=dvst(i)/sqrt(length(T_RTD));
    uB(i)=a_tot(i)/sqrt(3);%incertezza calcolta considerando una
                           ...distribuzione uniforme  
    uC(i)=sqrt(uA(i)^2+uB(i)^2);%incertezza composta
    k=3; %fattore di copertura
    U(i)=k*uC(i); %incertezza estesa
end
%ISTOGRAMMI 
%RTD
subplot(1,3,1);  
str_RTD=histogram(T_RTD,7,'normalization','pdf');
hold on
x = [min(T_RTD):0.01:max(T_RTD)]; % vettore delle ascisse rispetto a cui si calcola la gaussiana
f_gauss = exp(-(x-Tm_RTD).^2./(2*dvst_RTD^2))./(dvst_RTD*sqrt(2*pi));%gaussiana
plot(x,f_gauss,'-r')
xlabel('T[°C]')
ylabel('pdf')
title('T_.RTD');
hold off
%NTC
subplot(1,3,2); 
str_NTC=histogram(T_NTC,8,'normalization','pdf');
hold on
x = [min(T_NTC):0.01:max(T_NTC)];
f_gauss = exp(-(x-Tm_NTC).^2./(2*dvst_NTC^2))./(dvst_NTC*sqrt(2*pi));
plot(x,f_gauss,'-r')
xlabel('T[°C]')
ylabel('pdf')
title('T_.NTC')
hold off
%TC
subplot(1,3,3); 
str_TC=histogram(T_TC,7,'normalization','pdf');
hold on
x = [min(T_TC):0.01:max(T_TC)]; 
f_gauss = exp(-(x-Tm_TC).^2./(2*dvst_TC^2))./(dvst_TC*sqrt(2*pi));
plot(x,f_gauss,'-r')
xlabel('T[°C]')
ylabel('pdf')
title('T_.TC');
sgtitle('ISTOGRAMMI');
hold off
%compatibilità
figure(2)
hold on
errorbar(Tm_RTD,1,[],[],U(1),U(1),'ok')
errorbar(Tm_NTC,2,[],[],U(2),U(2),'ok')
errorbar(Tm_TC,3,[],[],U(3),U(3),'ok')
ylim([0 4])
xlabel('T(°C)')
set(gca,'YTick',[1 2 3 ],'YTickLabel',{'RTD' 'NTC' 'TC' });
grid on 
hold off
