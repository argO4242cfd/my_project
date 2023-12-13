clc
clear all
close all
a=6; b=0; c=2; d=3; e=7; f=6;
M=max([a,b,c,d,e,f]); m=(b+d+f)/3; q=(a+c+e)/3;
%funzioni di trasferimento
    %prima fdt
    alpha1=0.5*M;
    tau1=M+m+q;
    tau2=q;
    K=m+q;
    eta=m/q;
    P11=tf(K*[alpha1,1],[tau2^2, 2*tau2*eta, 1]);
    P12=tf(1,[tau1, 1]);
    P1=P11*P12;
%seconda fdt
    alpha2=0.25*M;
    P21=tf(K*[-alpha2,1],[tau1^2, 2*tau1 , 1]);
    P22=tf(1,[tau2, 1]);
    P2=P21*P22;
%calcolo dei poli,zeri,guadagno
Pm1=zpk(P1);
Pm2=zpk(P2);
[Zeri1,Poli1,K1] = zpkdata(Pm1);
guadagno1=dcgain(Pm1);
[Zeri2,Poli2,K2] = zpkdata(Pm2);
guadagno2=dcgain(Pm2);
%Studiare numericamente la risposta
...in funzione del tempo per un ingresso a gradino.
figure(1)
step(P1,P2); 
grid on
legend('P1','P2')
title('Risposta al gradino')




%%
%TERZO ESERCIZIO
%%prima fdt
Kc1=1;
C1=tf(Kc1,1);
G1=P1*C1;
figure(2)
rlocus(G1)
title('Luogo delle radici P1')
%seconda fdt
Kc2=1;
C2=tf(Kc2,1);
G2=P2*C2;
figure(3)
rlocus(G2)
title('Luogo delle radici P2')
%ripetere per un controllore proporzionale integrale 


%controllore PI 1
tauI1=20;
Kc11=1;
CI1=tf([tauI1*Kc11,Kc11],[tauI1,0]);
GI1=CI1*P1;
figure(4)
subplot(1,2,1)
rlocus(GI1)
title('Luogo delle radici PI 1 con tau=20')
tauI1=40;
Kc11=1;
CI1=tf([tauI1*Kc11,Kc11],[tauI1,0]);
GI1=CI1*P1;
subplot(1,2,2)
rlocus(GI1)
title('Luogo delle radici PI 1 con tau=40')
%controllore PI 2
tauI2=0.5;
Kc22=1;
CI2=tf([tauI2*Kc22,Kc22],[tauI2,0]);
GI2=CI2*P2;
figure(5)
rlocus(GI2)
title('Luogo delle radici PI 2')

%controllore PDI
%secondo PDI
KI1=40;
KD1=0.5;
Kc111=1;
CID1=tf([KI1*KD1*Kc111,Kc111*KI1,Kc111],[KI1,0]);
GID1=CID1*P1;
figure(6)
subplot(1,2,1)
rlocus(GID1)
title('Luogo delle radici PID 1 con tauI=40 e tauD=0.5')
KI1=40;
KD1=1;
Kc111=1;
CID1=tf([KI1*KD1*Kc111,Kc111*KI1,Kc111],[KI1,0]);
GID1=CID1*P1;
subplot(1,2,2)
rlocus(GID1)
title('Luogo delle radici PID 1 tauI=40 e tauD=1')

%2 PDI
KI2=0.5;
KD2=0.3;
Kc222=1;
CID2=tf([KD2*Kc222*KI2,Kc222*KI2,Kc222],[KI2,0]);
GID2=CID2*P2;
figure(7)
rlocus(GID2)
title('Luogo delle radici PID 2')



%%
%ESERCIZIO 4
%prima FDT
 figure(14)
opt=bodeoptions;      % Opzioni di Bode di default
opt.MagUnits='abs'; % Valore assoluto di Bode. Fa diventare lineare la scala
opt.MagScale='log'; % La scala di Bode diventa assoluta (comando precedente) e logaritmica
opt.PhaseMatching='on';
subplot(1,2,1)
bode(P1,opt)
grid on
subplot(1,2,2)
 nyquist(P1)
 sgtitle('Funzione di trasfermimento P1')
 [Mg1,MP1,w_mg1,w_mf1]=margin(P1);
% Definizione P2
 theta=0.25*M;
 F21=tf(K,[tau1^2,2*tau1,1]);
 F22=tf(1,[tau2,1]);
 [n_pade, d_pade] = pade(theta, 3);
 F23= tf(n_pade, d_pade);
 F2=F21*F22*F23;
figure(15)
opt=bodeoptions;      % Opzioni di Bode di default
opt.MagUnits='abs'; % Valore assoluto di Bode. Fa diventare lineare la scala
opt.MagScale='log'; % La scala di Bode diventa assoluta (comando precedente) e logaritmica
opt.PhaseMatching='on';
 subplot(1,2,1)
 bode(F2,opt)
 grid on
 subplot(1,2,2)
 nyquist(F2)
sgtitle('Funzione di trasfermimento P2')
 [Mg2,MP2,w_mg2,w_mf2]=margin(F2);
%CON REGOLATORE PI
%PRIMO
tau_4_I1=40;
Kc_4_1=1;
CI_4_1=tf([tau_4_I1*Kc_4_1,Kc_4_1],[tau_4_I1,0]);
GI_4_1=CI_4_1*P1;
figure(16)
opt=bodeoptions;      % Opzioni di Bode di default
opt.MagUnits='abs'; % Valore assoluto di Bode. Fa diventare lineare la scala
opt.MagScale='log'; % La scala di Bode diventa assoluta (comando precedente) e logaritmica
opt.PhaseMatching='on';
 subplot(1,2,1)
 bode(GI_4_1,opt)
 grid on
 subplot(1,2,2)
 nyquist(GI_4_1)
sgtitle('Funzione di trasfermimento PI1')
 [Mg1PI,MP1PI,w_mg1PI,w_mf1PI]=margin(GI_4_1);
 %SECONDO
tau_5_I1=0.1;
Kc_5_1=1;
CI_5_1=tf([tau_5_I1*Kc_5_1,Kc_5_1],[tau_5_I1,0]);
GI_5_1=CI_5_1*F2;
figure(17)
opt=bodeoptions;      % Opzioni di Bode di default
opt.MagUnits='abs'; % Valore assoluto di Bode. Fa diventare lineare la scala
opt.MagScale='log'; % La scala di Bode diventa assoluta (comando precedente) e logaritmica
opt.PhaseMatching='on';
 subplot(1,2,1)
 bode(GI_5_1,opt)
 grid on
 subplot(1,2,2)
 nyquist(GI_5_1)
sgtitle('Funzione di trasfermimento PI2')
 [Mg2PI,MP2PI,w_mg2PI,w_mf2PI]=margin(GI_5_1);


