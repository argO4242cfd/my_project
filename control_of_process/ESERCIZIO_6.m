clc
clear all
close all
a=6; b=0; c=2; d=3; e=7; f=6;
M=7; m=(b+d+f)/3; q=(a+c+e)/3;
%funzioni di trasferimento
%prima fdt
alpha1=0.5*M;
tau1=M+m+q;
tau2=q;
K=m+q;
eta=m/q;
P11=tf(K*[alpha1,1],[tau2^2, 2*tau2*eta, 1]);
P12=tf(1,[tau1, 1]);
P=P11*P12;
%seconda fdt
 theta=0.25*M;
 P21=tf(K,[tau1^2,2*tau1,1]);
 P22=tf(1,[tau2,1],'IOdelay',theta);
 Pd=P21*P22;
 %CONTROLLORE PID 1
KCC1=K;
thetaCC1=4;
tauCC1=14;
P1a=tf(KCC1,[tauCC1,1],'IOdelay',thetaCC1);
KcPID=tauCC1/(KCC1*thetaCC1)*(4/3+thetaCC1/(4*tauCC1));
tauIPID=thetaCC1*(32+6*thetaCC1/tauCC1)/(13+8*thetaCC1/tauCC1);
tauDPID=thetaCC1*(4/(11+2*thetaCC1/tauCC1));
C=tf([KcPID*tauDPID*tauIPID,KcPID*tauIPID,KcPID],[tauIPID,0]);
G1=C*P;
T1=G1/(1+G1);
S1=1/(1+G1);


%SENZA FEEDFORWARD
Cff=minreal(-Pd/P);
Sy=Pd/(1+P*C);
Su=-(C*Pd)/(1+P*C);
figure(1)
step(Sy)
hold on 
step(Su)
grid on
title('Risposta in anello chiuso senza controllo in avanti')
legend('uscita y(t)','ingresso u(t)')
% CON FEEDFORWARD

Syff=(Pd+P*Cff)/(1+P*C);
Suff = (Cff-C*Pd)/(1+P*C) ;

Su1 = -C*(Pd+Cff*P)/(1+P*C) ;
Su2 = Cff;

figure(2)
step(Syff)
hold on 
step(Sy)
grid on
title('Confronto di y(t) con e senza feedforward')
legend('Con feedforward','Senza feedforward')

figure(3)
step(Suff)
hold on
step(Syff)
grid on
legend('Ingresso u(t)','Uscita y(t)')
title('Confronto tra ingresso e uscita con feedforward')



figure(4)
step(Suff)
hold on
step(Su1)
hold on 
step(Su2)
grid on
title('Confronto degli ingressi')
legend('u(t)','u1(t)','u2(t)')

figure(5)
subplot(1,2,1)
step(Su1)
grid on
title('Risposta controllore feedback')
subplot(1,2,2)
step(Su2)
grid on
title('Risposta controllore feedforward')
sgtitle('Confronto degli ingressi')







