clc
clear all
close all
a=6; b=0; c=2; d=3; e=7; f=6;
M=7; m=(b+d+f)/3; q=(a+c+e)/3;
alpha=0.5*M;
tau1=M+m+q;
tau2=q;
K=m+q;
eta=m/q;
theta=0.25*M;
%prima fdt
P11=tf(K*[alpha,1],[tau2^2, 2*tau2*eta, 1]);
P12=tf(1,[tau1, 1]);
P1=P11*P12;
%seconda fdt
P21=tf(K,[tau1,1]);
P22=tf(1,[tau2^2,2*tau2,1]);
[n_pade, d_pade] = pade(theta, 3);
P23= tf(n_pade, d_pade);
P2=P21*P22*P23;
%RICORDA: l'approssimazione di pade è stata fatta perchè nel tuning ZN
%serve il luogo delle radici, per cui la funzione rlocus funziona solo con
%pade e non con IOdelay

%primo controllore per P1
     %cohen
     %fdt approssimata
    KCC1=K;
    thetaCC1=4;
    tauCC1=14;
    P1a=tf(KCC1,[tauCC1,1],'IOdelay',thetaCC1);
    figure(3)
    step(P1)
    hold on
    title('Confronto di P1 e P1a ')
    step(P1a)
    legend('P1','P1approssimata')
    grid on
%con questi valori della funzione approssimata si fa il tuning del
%controllore PI
Kcc=tauCC1/(KCC1*thetaCC1)*(0.9+thetaCC1/(12*tauCC1));
tauIc=thetaCC1*(30+3*thetaCC1/tauCC1)/(9+20*thetaCC1/tauCC1);
Cc=Kcc*tf([tauIc,1],[tauIc,0]);
Gcc=P1*Cc;
Tcc=Gcc/(1+Gcc);
Scc=1/(1+Gcc);
%Skogestad
KSK=(tauCC1/KCC1)*(1/(2*thetaCC1));
tauSK=min(tauCC1, 8*thetaCC1);
CSK=KSK*tf([tauSK,1],[tauSK,0]);
GSK=P1*CSK;
TSK=GSK/(1+GSK);
figure(4)
step(Tcc)
hold on
step(TSK)
grid on
title('Confronto delle risposte dei controllori PI per T1')
legend('Risposta Cohen-Coon','Risposta Skogestad')
%secondo controllore per P1

%fdt approssimata
% KCC2=K;
% thetaCC2=10;
% tauCC2=17;
KCC2=K;
thetaCC2=theta+3/2*tau2;
tauCC2=tau1+tau2/2;
P2a=tf(KCC2,[tauCC2,1],'IOdelay',thetaCC2);
figure(5)
step(P2)
hold on
title('Confronto di P2 e P2a ')
step(P2a)
legend('P2','P2approssimata')
grid on
 figure(7)
 rlocus(P2)
        KCas=0.765;
        omega=0.196;
Pu=2*pi/omega;
KZN=KCas/2.2;
tauZN=Pu/1.2;
CZN=KZN*tf([tauZN,1],[tauZN,0]);
GZN=P2*CZN;
TZN=GZN/(1+GZN);
%skogestad

KSK2=(tauCC2/KCC2)*(1/(2*thetaCC2));
tauSK2=min(tauCC2, 8*thetaCC2);
CSK2=KSK2*tf([tauSK2,1],[tauSK2,0]);
GSK2=P2*CSK2;
TSK2=GSK2/(1+GSK2);
figure(8)
step(TZN)
hold on
step(TSK2)
title('Confronto delle risposte dei controllori PI per T2')
legend('Risposta Ziegler Nicholson','Risposta Skogestad')
grid on
SSKPI2=1/(1+GSK2);
%controllore PID
%prima fdt 
%partendo dai dati delle approssimazioni fatte in precedenza

    KcPID=tauCC1/(KCC1*thetaCC1)*(4/3+thetaCC1/(4*tauCC1));
    tauIPID=thetaCC1*(32+6*thetaCC1/tauCC1)/(13+8*thetaCC1/tauCC1);
    tauDPID=thetaCC1*(4/(11+2*thetaCC1/tauCC1));
    CccPID=tf([KcPID*tauDPID*tauIPID,KcPID*tauIPID,KcPID],[tauIPID,0]);
    GPID=CccPID*P1;
    TPID=GPID/(1+GPID);
    SPID=1/(1+GPID);
%SECONDA FDT con choen del PID , si parte sempre dalle approssimazioni
%fatte prima 
% KcPID2=tauCC2/(KCC2*thetaCC2)*(4/3+thetaCC2/(4*tauCC2));
% tauIPID2=thetaCC2*(32+6*thetaCC2/tauCC2)/(13+8*thetaCC2/tauCC2);
% tauDPID2=thetaCC2*(4/(11+2*thetaCC2/tauCC2));
% CccPID2=tf([KcPID2*tauDPID2*tauIPID2,KcPID2*tauIPID2,KcPID2],[tauIPID2,0]);
% GPID2=CccPID2*P2a;
% TPID2=GPID2/(1+GPID2);
% SPID2=1/(1+GPID2);
%OPPOURE APPROSSIMAZIONE CON SKOGESTAD
tau1ap=tau1;
tau2ap=tau2+0.5*tau2;
thetaap=theta+0.5*tau2;
P2_SOPTD=tf(K,[tau1ap*tau2ap,tau1ap+tau2ap,1],'IOdelay',thetaap);
figure(10)
step(P2)
hold on
step(P2_SOPTD)
grid on
legend('P2','P2_SOPTD')
title('Confronto di P2 con approssimazione del secondo ordine')
KSK_SOPTD=tau1ap/K*(1/(2*thetaap));
tauISK_SOPTD=min(tau1ap,8*thetaap);
tauDSK_SOPTD=tau2ap;
CSKPID2=tf([KSK_SOPTD*tauDSK_SOPTD*tauISK_SOPTD,KSK_SOPTD*tauISK_SOPTD,KSK_SOPTD],[tauISK_SOPTD,0]);
GSKPID2=CSKPID2*P2;
SSKPID2=1/(1+GSKPID2);
%%%%prima rispista
figure(11)
step(SPID)
hold on 
step(Scc)
grid on 
title('Confronto delle risposte dei controllori PI e PID per S1')
legend('PID-Choen Coon','PI-Choen Coon')
%seconda risposta
figure(12)
step(SSKPID2)
hold on
step(SSKPI2)
title('Confronto delle risposte dei controllori PI e PID per S2')
legend('PID-Skogestad','PI-Skogestad')
grid on


%TERZO QUESITO
%Confrontare i regolatori PID e PI del punto precedente nel caso di abbattimento di disturbo a gradino entrante
%sull%ingresso al processo
    %Si deve definire una nuova G 
  
%Prima fdt

G1PIi=P1/(1+P1*Cc);
G1PIDi=P1/(1+P1*CccPID);
% S1PIi=1/(1+G1PIi);
% S1PIDi=1/(1+G1PIDi);
figure(13)
step(G1PIi);
hold on
step(G1PIDi)
grid on
title('Abbattimento del disturbo entrante in ingresso per P1')
legend('PI-Cohen Coon','PID-Cohen Coon')

%seconda fdt
G2PIi=P1/(1+P2*CSK2);
G2PIDi=P1/(1+P2*CSKPID2);
%S2PIi=1/(1+G2PIi);
%S2PIDi=1/(1+G2PIDi);
figure(14)
step(G2PIi);
hold on
step(G2PIDi)
grid on
title('Abbattimento del disturbo entrante in ingresso per P2')
legend('PI-Skogestad','PID-Skogestad')




