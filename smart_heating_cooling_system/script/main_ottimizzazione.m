clearvars
close all
clc
tic

load Dati_energetici_grosseto.mat
load T_est.mat % [°C];
load Rad_glob.mat % [W/m^2]
load Rad_dir.mat % [W/m^2]
load Rad_diff.mat% [W/m^2]
clear V_wind; clear Clima; clear Months; clear Days; clear Hours; clear Um_rel
%------------------

%Porto i dati di ora in ora in dati di minuto in minuto, ipotizzo valore
%costante per tutto il minuto

T_est_min=zeros(60*length(T_est),1);
Rad_glob_min=zeros(60*length(T_est),1);
Rad_dir_min=zeros(60*length(T_est),1);
Rad_diff_min=zeros(60*length(T_est),1);

for ii=1:length(T_est)
    T_est_min(1+(ii-1)*60:60*ii)=repmat(T_est(ii),[60,1]);
    Rad_glob_min(1+(ii-1)*60:60*ii)=repmat(Rad_glob(ii),[60,1]);
    Rad_dir_min(1+(ii-1)*60:60*ii)=repmat(Rad_dir(ii),[60,1]);
    Rad_diff_min(1+(ii-1)*60:60*ii)=repmat(Rad_diff(ii),[60,1]);
    
end

%-------------------------------------

Cp=4.186; % Kj/Kg*k
T_acquedotto=15; %°C
Carico_elettrico=richieste_elettriche;
Carico_termico=portata_ACS.*Cp.*(T_mandata_ACS-T_acquedotto);

%------------------- DEFINIZIONE SCENARI-------------------- 

n_P_ibridi=[32,36,40,44,48,50,52,56,62,68,80,90,100];% numero di pannelli ibridi
Taglia_CHP=[13,16,23,34,37,40,44];% KW TAGLIA CHP
Tmin=30;% °C T minima TES
V_tes=[0.7,1,1.5,2];%,2];%,0.7,1,1.5]; % m^3 volume TES
T_max=55;%°C T massima per CHP
cont=0;
ro=1000;%Kg/m^3
m_2=[0.07,0.09,0.13,0.18,0.19,1,1];%,1,2.9];% Kg/s
eta_th_nom=[0.698,0.757,0.783,0.735,0.739,0.72,0.695];
eta_th_min=[0.35,0.35,0.35,0.35,0.35,0.35,0.35];
eta_el_nom=[0.318,0.316,0.312,0.32,0.335,0.3,0.315];

K_boll=[3,3.51,4.38,5.37];

%----------INIZIALIZZAZIONE VARIABILI A ZERO---------------

Q_term=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
Tw_out_ST=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
eta_th=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
Q_el_PV=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
eta_eff_PV=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));

T_TS=ones(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes)).*25; 
Tw_out_TES=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
Q_losses=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));

Th_en_eff=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
FC_th=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
T_CHP=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
El_en_eff=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
Pr_en_CHP=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
eta_th_eff=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
Th_en_res=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
El_en_res=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
El_en_surplus_CHP=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
El_en_surplus_PVT=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));

Carico_termico_CHP=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
Carico_elettrico_CHP=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
T_ritorno_CHP=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
m_eff_TES=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));

costo_gas=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
costo_rete_acquisto=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
costo_rete_vendita=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
CO2=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
CO2_CHP=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
CO2_rete=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));

costo_CHP=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
costo_TES=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
costo_PVT=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
costo_capitale_tot=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
costo_operativo_tot=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
costo_tot_year=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
year=20;
CO2_savings=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
CH4_CHP=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
Confronto_costi=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
costo_operativo_CHP=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
NPC=zeros(1,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));

int=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
p=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
m_CHP=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));

%-------------------------------------------------

for i=1:length(n_P_ibridi)
m_1(i)=0.02*(n_P_ibridi(i)/2);%Kg/s
end

mw_in_ST=zeros(length(Carico_termico),length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));
A=zeros(4,length(n_P_ibridi)*length(Taglia_CHP)*length(T_max)*length(V_tes));

V_tes_rif=0.3;% m^3
costo_rif_TES=600;% €
costo_unitario_PVT=800;% €
P_caldaia=200;% KW
costo_capitale_caldaia=2222+56*P_caldaia;% €
costo_operativo_caldaia=0.02*costo_capitale_caldaia;% €
costo_rif_CHP=2500*Taglia_CHP(1);% 2500 €/KW installato, prezzo massimo

Sup_pannello=1.6; %[m^2]
a2=0.007e-03;%KW/m^2*K^2
eta0=0.56;
eta_PV=0.17;% A T=25°C
eta_1=0.1475; 
eta_2=0.125;
T_est_1=30;
T_est_minima=35;
IAM=1.46;
pendenza_efficienza=0.0045; %eff/°C
a1=2.06e-03;%KW/m^2*K

T_r=15;%°C

PC=50/3.6;% PC metano in KWh/Kg, 50 MJ/Kg
price_gas=0.3;% da bolletta €/KWh 
price_el_vendita=0.1; % da GSE, sono €/KWh 
price_el_acquisto=0.2;% PUN da GME

costo_CH4_caldaia=(sum(Carico_termico)/(0.9*60))*price_gas;
costo_en_el_caldaia=(sum(Carico_elettrico)/(60))*price_el_acquisto;
CO2_caldaia=((sum(Carico_termico)/(0.9*60))/PC)*2.75+((sum(Carico_elettrico))/60)*0.3;


%--------------------- ANALISI DI OTTIMIZZAZIONE--------------------

for i=1:length(n_P_ibridi) % for dimensionamento pannelli ibridi
    for j=1:length(Taglia_CHP) % for dimensionamento CHP
        for k=1:length(T_max) % for T di setpoint
            for y=1:length(V_tes) % for V accumulo
            cont=cont+1; % contatore degli scenari
            mw_in_ST(:,cont)=m_1(i); % portata nei pannelli
            m_CHP(:,cont)=m_2(j);% portata nel CHP
            A(:,cont)=[n_P_ibridi(i),Taglia_CHP(j),T_max(k),V_tes(y)]'; % matrice con le condizioni di ogni scenario  
       
%----------------SIMULAZIONE NELL'ANNO MINUTO PER MINUTO------------

            for minute=1:length(Carico_termico)
               
                %----------------SOLARE IBRIDO-------------
                
                [Q_term(minute,cont),Tw_out_ST(minute,cont),eta_th(minute,cont),Q_el_PV(minute,cont),eta_eff_PV(minute,cont),p(minute,cont)] = ...
                    PVT_2(Sup_pannello,n_P_ibridi(i),T_TS(minute,cont),mw_in_ST(minute,cont),Rad_glob_min(minute)/1000,T_est_min(minute),eta0,IAM,a1,a2, eta_1, eta_2, T_est_1,T_est_minima);
                
                if p(minute,cont)==0
                    mw_in_ST(minute,cont)=0;
                end
                
                %-------------------CHP----------------------------
               
                Q_CHP_nom=Taglia_CHP(j)*eta_th_nom(j);
                E_CHP_nom=Taglia_CHP(j)*eta_el_nom(j);
                El_en_surplus_PVT(minute,cont)=(max([0,-richieste_elettriche(minute)+Q_el_PV(minute,cont)]));
                
                %-----------LOGICA GESTIONE CHP--------------------

                if minute==1 && T_TS(minute,cont)<T_max(k)
                    int(minute+1,cont)=1;
                else
                    if T_TS(minute,cont)<=Tmin || (T_TS(minute,cont)<T_max(k) && int(minute,cont)==1) || T_mandata_ACS(minute)>T_TS(minute-1,cont)
                        int(minute+1,cont)=1;
                        Carico_termico_CHP(minute,cont)=min([Q_CHP_nom, 1000*4.2*V_tes(y)*(T_max(k)-T_TS(minute,cont))/60]); 
                        Carico_elettrico_CHP(minute,cont)=(max([0,richieste_elettriche(minute)-Q_el_PV(minute,cont)])); 
                        T_ritorno_CHP(minute,cont)=T_TS(minute,cont)+5;
                    elseif T_TS(minute,cont)>T_max(k) || (T_TS(minute,cont)<T_max(k) && int(minute,cont)==0)
                        int(minute+1,cont)=0;
                        Carico_termico_CHP(minute,cont)=0; 
                        Carico_elettrico_CHP(minute,cont)=0;
                        m_CHP(minute,cont)=0;
                        T_ritorno_CHP(minute,cont)=T_TS(minute,cont);
                    else
                        error('errore')
                    end
                end
                %--------------------------------------------

                [Th_en_eff(minute,cont), FC_th(minute,cont), T_CHP(minute,cont), El_en_eff(minute,cont), Pr_en_CHP(minute,cont), eta_th_eff(minute,cont), ...
                    Th_en_res(minute,cont),El_en_res(minute,cont),El_en_surplus_CHP(minute,cont)] = ...
                    CHP_FTL_2(Carico_termico_CHP(minute,cont),Carico_elettrico_CHP(minute,cont),T_ritorno_CHP(minute,cont), m_CHP(minute,cont),eta_th_nom(j), Q_CHP_nom, E_CHP_nom, Q_CHP_nom*0.2, eta_th_min(j),Carico_elettrico(minute),Q_el_PV(minute,cont));
                
                
                %------------ACCUMULO TERMICO-----------------
                
                [T_TS(minute+1,cont),Tw_out_TES(minute+1,cont),Q_losses(minute,cont),m_eff_TES(minute,cont)] = accumulo_termico_2(mw_in_ST(minute,cont),Tw_out_ST(minute,cont),m_CHP(minute,cont),...
                T_CHP(minute,cont),portata_ACS(minute),T_acquedotto,K_boll(y),T_est_min(minute),T_TS(minute,cont),V_tes(y),T_mandata_ACS(minute));
                
                
                %------------ COSTI, CO2, CH4 -----------------------------
                

            end

            %----------------- FINE SIMULAZIONE ANNUALE PER UN SINGOLO SCENARIO ------------------

            costo_gas(1,cont)=sum((Pr_en_CHP(:,cont)/60)*price_gas); 
            costo_rete_acquisto(1,cont)=sum((El_en_res(:,cont)/60)*price_el_acquisto);  
            costo_rete_vendita(1,cont)=sum(((El_en_surplus_CHP(:,cont)+El_en_surplus_PVT(:,cont))/60)*price_el_vendita);
            CO2(1,cont)=sum(((Pr_en_CHP(:,cont)/60)/PC)*2.75)+sum((El_en_res(:,cont)/60)*0.3);
            CO2_CHP(1,cont)=sum(((Pr_en_CHP(:,cont)/60)/PC)*2.75);
            CO2_rete(1,cont)=sum((El_en_res(:,cont)/60)*0.3);
            costo_CHP(cont)=7789*(E_CHP_nom)^0.6;
            costo_operativo_CHP(cont)=0.015*(sum(Carico_elettrico_CHP(:,cont))/60);
            costo_TES(cont)=costo_rif_TES*(V_tes(y)/V_tes_rif)^0.7;
            costo_PVT(cont)=costo_unitario_PVT*n_P_ibridi(i);
            costo_capitale_tot(cont)=costo_CHP(cont)+costo_TES(cont)+costo_PVT(cont);
            costo_operativo_tot(cont)=costo_operativo_CHP(cont)+costo_rete_acquisto(1,cont)+costo_gas(1,cont)-costo_rete_vendita(1,cont);
            costo_tot_year(cont)=costo_capitale_tot(cont)+costo_operativo_tot(cont)*year;
            CO2_savings(cont)=CO2_caldaia-CO2(1,cont); % risparmio CO2 caso integrato rispetto a solo caldaia
            CH4_CHP(cont)=sum(((Pr_en_CHP(:,cont)/60)/PC));
            Confronto_costi(cont)=costo_tot_year(cont)/(costo_capitale_caldaia+year*(costo_CH4_caldaia+costo_en_el_caldaia+costo_operativo_caldaia));
            c=0;
            for s=1:20 % year=20
              c=c+(1/(1+0.05))^s;% 0.05 tasso di interesse
            end
            NPC(cont)=costo_capitale_tot(cont)+(costo_operativo_tot(cont))*c;
            
           end
       end
   end
end

%------------------- FINE DELL'ANALISI ANNUALE DI TUTTI GLI SCENARI ---------

NPC_caldaia=costo_capitale_caldaia+(costo_CH4_caldaia+costo_en_el_caldaia+costo_operativo_caldaia)*c;

[CO2_ordinata,z]=sort(CO2,'descend');
%NPV_diff=NPV-NPV_caldaia;
NPC_ordinato=NPC(z);

%-------------------- FRONTE DI PARETO --------------------------

figure 
plot(CO2_ordinata*year/1000,NPC_ordinato,'bd')
title(strcat('Risultati ottimizzazione esaustiva, numero di scenari=',num2str(length(A)),''))
xlabel('Ton CO2 prodotti')
ylabel('€ NPC sistema integrato')
legend('dati scenari analizzati')
grid on
box on
hold on

 %----------------- IDENTIFICAZIONE PUNTI SUL FRONTE SINISTRO E DESTRO -------------------

[CO2_min,i]=min(CO2*year/1000);
[NPC_min,j]=min(NPC);
r=0;
index=[];
CO2_supp=CO2*year/1000;
NPC_supp=NPC;
while i~=j
r=r+1;
for s=1:length(CO2)

     if NPC_supp(s)>NPC_supp(i) && CO2_supp(s)>CO2_supp(i) 
       NPC_supp(s)=2e05;
       CO2_supp(s)=10000;
     end

 end
 NPC_supp(i)=2e05;
 CO2_supp(i)=10000;
 index(r)=i;
 [CO2_min,i]=min(CO2_supp);

end
[CO2_min,i]=min(CO2_supp);
index(r+1)=i;
A_ott=[];
for i=1:length(index)
A_ott(:,i)=A(:,index(i));
end
CO2_pareto=[];
NPC_pareto=[];
for k=1:length(index)
    CO2_pareto(k)=CO2(index(k))*year/1000;
    NPC_pareto(k)=NPC(index(k));

 txt = ['n^° PVT = ',num2str(A(1,index(k))),', taglia CHP=',num2str(A(2,index(k))),'KW, T_m_a_x=',num2str(A(3,index(k))),'°C, V_t_e_s=',num2str(A(4,index(k))),'m^3'];
        %plot(rand(10,1),'DisplayName',txt)
        plot(CO2_pareto(k),NPC_pareto(k),'r*','DisplayName',txt)
end
legend show

[CO2_max,x]=max(CO2*year/1000);
[NPC_max,y]=max(NPC);
CO2_supp2=CO2*year/1000;
NPC_supp2=NPC;
r=0;
index2=[];
while y~=x
r=r+1;
for s=1:length(CO2)

     if NPC_supp2(s)<NPC_supp2(y) && CO2_supp2(s)<CO2_supp2(y) 
       NPC_supp2(s)=0;
       CO2_supp2(s)=0;
     end

end
NPC_supp2(y)=0;
 CO2_supp2(y)=0;


 index2(r)=y;
 [NPC_max,y]=max(NPC_supp2);

end
[NPC_max,y]=max(NPC_supp2);
index2(r+1)=y;
A_pegg=[];
for i=1:length(index2)
A_pegg(:,i)=A(:,index2(i));
end

CO2_pegg=[];
NPC_pegg=[];
for k=1:length(index2)
    CO2_pegg(k)=CO2(index2(k))*year/1000;
    NPC_pegg(k)=NPC(index2(k));

 txt = ['n^° PVT = ',num2str(A(1,index2(k))),', taglia CHP=',num2str(A(2,index2(k))),'KW, T_m_a_x=',num2str(A(3,index2(k))),'°C, V_t_e_s=',num2str(A(4,index2(k))),'m^3'];
        %plot(rand(10,1),'DisplayName',txt)
        plot(CO2_pegg(k),NPC_pegg(k),'g*','DisplayName',txt)
end

legend show
hold off

%------------------------- PROPOSTA 3 SCENARI MIGLIORI ------------------------------------

for i=1:length(index)
CO2_pareto_ad(i)=(CO2_pareto(i)-min(CO2_pareto))/(max(CO2_pareto)-min(CO2_pareto));
NPC_pareto_ad(i)=(NPC_pareto(i)-min(NPC_pareto))/(max(NPC_pareto)-min(NPC_pareto));
end
figure
plot([min(CO2_pareto_ad),max(CO2_pareto_ad)],[min(NPC_pareto_ad),min(NPC_pareto_ad)],'b-.')
hold on
plot([min(CO2_pareto_ad),min(CO2_pareto_ad)],[min(NPC_pareto_ad),max(NPC_pareto_ad)],'b-.')
plot(CO2_pareto_ad,NPC_pareto_ad,'rd')

for i=1:length(CO2_pareto_ad)
d(i)=sqrt((CO2_pareto_ad(i)-min(CO2_pareto_ad))^2+(NPC_pareto_ad(i)-min(NPC_pareto_ad))^2);
end
[d_min,indice]=mink(d,3); % 'indice' da gli indici delle distanze minori, relativi ai punti migliori del pareto(CO2_pareto)

plot([min(CO2_pareto_ad),CO2_pareto_ad(indice(1))],[min(NPC_pareto_ad),NPC_pareto_ad(indice(1))],'b-.')
hold on
plot([min(CO2_pareto_ad),CO2_pareto_ad(indice(2))],[min(NPC_pareto_ad),NPC_pareto_ad(indice(2))],'r-.')
plot([min(CO2_pareto_ad),CO2_pareto_ad(indice(3))],[min(NPC_pareto_ad),NPC_pareto_ad(indice(3))],'g-.')
plot(CO2_pareto_ad(indice(1)),NPC_pareto_ad(indice(1)),'bo')
plot(CO2_pareto_ad(indice(2)),NPC_pareto_ad(indice(2)),'ro')
plot(CO2_pareto_ad(indice(3)),NPC_pareto_ad(indice(3)),'go')
axis equal
grid on
box on
title('Proposta di scelta scenari ottimi, minore distanza dall''origine')
xlim([0 1])
ylim([0 1])
hold off

%-------------------- ORE DI DISSERVIZIO ACS --------------------------------

h_disserv=zeros(1,length(index));
for s=1:length(index)
    for m=1:length(Carico_elettrico)
    if T_TS(m,index(s))<T_mandata_ACS(m)
        h_disserv(1,s)=h_disserv(1,s)+1;
    end
    end
end
figure
h_disserv=round(h_disserv/60);
X = categorical({'A','B','C','D','E','F','G','H','I','L','M','N','O','P','Q','R'});
X = reordercats(X,{'A','B','C','D','E','F','G','H','I','L','M','N','O','P','Q','R'});
for i=1:length(index)
b(i)=bar(X(i),h_disserv(i),'b');

hold on
end
grid on
title('Ore di disservizio ACS in un anno di funzionamento')
xlabel('Scenari ottimi')
ylabel('Ore')
print('Ore di disservizio','-djpeg')
hold off

%--------------------- ORE DI FUNZIONAMENTO CHP -----------------------------

h_chp=zeros(1,length(index));
for s=1:length(index) 
          for y=1:length(Carico_elettrico)
          if   m_CHP(y,index(s))>0
              h_chp(1,s)=h_chp(1,s)+1;
          end
          end
end
h_chp=round(h_chp/60);
figure
for i=1:length(index)
b(i)=bar(X(i),h_chp(i),'b');
xtips1 = b(i).XEndPoints;
ytips1 = b(i).YEndPoints;
labels1 = string(b(i).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
hold on
end
grid on
title('Ore di funzionamento CHP in un anno di funzionamento')
xlabel('Scenari ottimi')
ylabel('Ore')
hold off

%----------------------------- OPEX/CAPEX -----------------------------

figure 
hold on
b(1)=bar(1,round(costo_capitale_tot(index(1))/1000,1),0.5,'b');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
b(1)=bar(1+0.2,round(costo_operativo_tot(index(1))*c/1000,1),0.5,'r'); 
xtips2 = b(1).XEndPoints;
ytips2 = b(1).YEndPoints;
labels2 = string(b(1).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
for i=1:length(index)
b(i)=bar(i,round(costo_capitale_tot(index(i))/1000,1),0.5,'b');
xtips1 = b(i).XEndPoints;
ytips1 = b(i).YEndPoints;
labels1 = string(b(i).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
b(i)=bar(i+0.2,round(costo_operativo_tot(index(i))*c/1000,1),0.5,'r');
xtips2 = b(i).XEndPoints;
ytips2 = b(i).YEndPoints;
labels2 = string(b(i).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
end
hold off
grid on
title('Confronto costi capitali e operativi totali attualizzati per le soluzioni ottime')
xlabel('Scenari ottimi')
legend('Opex','Capex')
ylabel('K€')
print('capexopex','-djpeg')

%----------------------- EMISSIONI CO2 RETE e CHP -----------------------

figure
hold on
b(1)=bar(1,round(CO2_rete(index(1)),0),0.5,'r');
xtips2 = b(1).XEndPoints;
ytips2 = b(1).YEndPoints;
labels2 = string(b(1).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
b(1)=bar(1+0.2,round(CO2_CHP(index(1)),0),0.5,'b');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
for i=2:length(index)
b(i)=bar(i,round(CO2_rete(index(i)),0),0.5,'r');
xtips2 = b(i).XEndPoints;
ytips2 = b(i).YEndPoints;
labels2 = string(b(i).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
b(i)=bar(i+0.2,round(CO2_CHP(index(i)),0),0.5,'b');
xtips1 = b(i).XEndPoints;
ytips1 = b(i).YEndPoints;
labels1 = string(b(i).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
end
hold off
grid on
title('Confronto emissioni annuali di CO2 dovute al CHP e all''acquisto di energia elettrica')
xlabel('Scenari ottimi')
legend('CO2 rete','CO2 CHP')
ylabel('Kg CO2 emessa')
print('CO2 confronto','-djpeg')

%--------------------- CO2 EMESSA IN 20 ANNI SCENARI OTTIMI --------------------------

figure
subplot(2,1,1)
for i=1:length(index)
b(i)=bar(X(i),round(CO2(index(i))*year/1000,1),'b');
xtips1 = b(i).XEndPoints;
ytips1 = b(i).YEndPoints;
labels1 = string(b(i).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
hold on
end
hold off
grid on
title('CO2 emessa scenari ottimi in 20 anni di funzionamento impianto')
xlabel('Scenari ottimi')
ylabel('ton CO2 emessa')

%--------------------- NPC in 20 ANNI SCENARI OTTIMI ----------------------------------

subplot(2,1,2)
for i=1:length(index)
b(i)=bar(X(i),round(NPC(index(i))/1000,1),'r');
xtips1 = b(i).XEndPoints;
ytips1 = b(i).YEndPoints;
labels1 = string(b(i).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
hold on
end
hold off
grid on
title('NPC scenari ottimi in 20 anni di funzionamento impianto')
xlabel('Scenari ottimi')
ylabel('K€')
print('particolare punti del pareto','-djpeg')

%---------------------- CONFRONTO NPC CALDAIA/SISTEMA INTEGRATO -------------------------------

CH4_caldaia=((sum(Carico_termico)/(0.9*60))/PC);

figure 
for i=1:length(index)
    b(i)=bar(X(i),round(NPC(index(i))/NPC_caldaia,2),'b');
xtips1 = b(i).XEndPoints;
ytips1 = b(i).YEndPoints;
labels1 = string(b(i).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
hold on
end
grid on
title('rapporto tra il Net present Cost caldaia e sistema integrato')
xlabel('Scenari ottimi')
ylabel('NPC_i_n_t_e_g_r_a_t_o/NPC_c_a_l_d_a_i_a')
hold off
print('NPC confronto','-djpeg')
toc

%------------------- CONFRONTO CAPEX CALDAIA/SISTEMA INTEGRATO -----------------------------

figure 
T=[];
for i=1:length(index)
T=[T;costo_CHP(index(i))/1000,costo_TES(index(i))/1000,costo_PVT(index(i))/1000];
end
bar(T,'stacked')
hold on
X1 = categorical({'Caldaia'});
X1 = reordercats(X1,{'Caldaia'});
b(1)=bar(18,round(costo_capitale_caldaia/1000,0),'r');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
hold off
grid on
title('Confronto tra i costi di investimento impianto tradizionale vs integrato')
xlabel('Scenari ottimi')
ylabel('K€')
legend('Costo capitale CHP','Costo capitale accumulo','Costo capitale PVT','Costo capitale caldaia')

%-------------------- CONFRONTO COSTI OPERATIVI CALDAIA/SISTEMA INTEGRATO -----------------

figure 
T2=[];
for i=1:length(index)
T2=[T2;year*costo_operativo_CHP(index(i))/1000,year*costo_gas(index(i))/1000,year*costo_rete_acquisto(index(i))/1000,-year*costo_rete_vendita(index(i))/1000];
end
bar(T2,'stacked')
hold on
b(1)=bar(18,round(year*(costo_CH4_caldaia+costo_en_el_caldaia+costo_operativo_caldaia)/1000,0),'r');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
grid on
title('Confronto costi operativi non attualizzati complessivi sui 20 anni sistema tradizionale vs sistema integrato')
xlabel('Scenari ottimi')
ylabel('K€')
legend('Costo operativo CHP','Costo acquisto gas','Costo acquisto energia elettrica','Guadagno energia elettrica venduta','Costi operativi totali caldaia')
hold off

%--------------------- BILANCI ENERGIA ELETTRICA CASI OTTIMI ----------------------

figure
t=[];
for i=1:length(index)
t=[t;-sum(El_en_surplus_PVT(:,index(i)))/60/1000,-sum(El_en_surplus_CHP(:,index(i)))/60/1000,sum(El_en_eff(:,index(i)))/60/1000,sum(Q_el_PV(:,index(i)))/60/1000,sum(El_en_res(:,index(i)))/60/1000];
end
bar(t,'stacked')%X,
legend('Surplus PVT','Surplus CHP','En prodotta CHP','En prodotta PVT','En comprata dalla rete')
grid on
title('Bilancio energia elettrica casi ottimali')
xlabel('Scenari ottimi')
ylabel('MWh')

%------------------- CH4 RISPARMIATA SISTEMA INTEGRATO/CALDAIA -----------------------

figure
for i=1:length(index)
CH4_savings_perc(i)=(CH4_caldaia-CH4_CHP(index(i)))/CH4_caldaia*100;
end
bar(X,CH4_savings_perc,'b')
%legend('Surplus PVT','Surplus CHP','En prodotta CHP','En prodotta PVT','En comprata dalla rete')
grid on
title('% CH_4 risparmiata in un anno con sistema integrato rispetto a sistema tradizionale')
ylabel('%')
xlabel('Scenari ottimi')

%------------------ CO2 RISPARMIATA SISTEMA INTEGRATO/CALDAIA--------------

CO2_caldaia_tot=CO2_caldaia+sum(Carico_elettrico)/60*0.3;

figure
for i=1:length(index)
CO2_integrato_tot(i)=CO2_CHP(index(i))+CO2_rete(index(i));
CO2_risparmiata_perc(i)=(CO2_caldaia_tot-CO2_integrato_tot(i))/CO2_caldaia_tot*100;
bar(X(i),CO2_risparmiata_perc(i))
hold on
end
title('Percentuale di CO2 risparmiata con sistema integrato rispetto sistema tradizionale in 1 anno')
xlabel('scenari ottimi')
ylabel('%')
grid on

%-----------------BILANCIO TERMICO CASI OTTIMI --------------------------

Uscita_CHP=zeros(1,length(index));Ingresso_CHP=zeros(1,length(index));Carico_effettivo_ACS=zeros(1,length(index));
Uscita_Solare=zeros(1,length(index));Ingresso_Solare=zeros(1,length(index));Perdite_TS=zeros(1,length(index));
Bilancio_energia_interna=zeros(1,length(index));Bilancio=zeros(1,length(index));

for i=1:length(index)
Uscita_CHP(i)= sum(m_CHP(:,index(i)).*4.2.*T_CHP(:,index(i)))/60; 
Ingresso_CHP(i)= sum(m_CHP(:,index(i)).*4.2.*(T_TS(1:end-1,index(i))+5))/60; 
Uscita_Solare(i) = sum(mw_in_ST(:,index(i)).*4.2.*Tw_out_ST(:,index(i)))/60; 
Ingresso_Solare(i)= sum(mw_in_ST(:,index(i)).*4.2.*T_TS(1:end-1,index(i)))/60; 
Carico_effettivo_ACS(i)= sum(m_eff_TES(:,index(i)).*4.2.*(T_TS(1:end-1,index(i))-T_acquedotto))/60; 
Perdite_TS(i)= sum(K_boll(end).*(T_TS(1:end-1,index(i))-T_est_min(:))/1000/60);
Bilancio_energia_interna(i)=2000*4.2*(T_TS(end,index(i))-T_TS(1,index(i)))/3600;
Bilancio(i)=(Uscita_CHP(i)-Ingresso_CHP(i))+(Uscita_Solare(i)-Ingresso_Solare(i))-Perdite_TS(i)-Carico_effettivo_ACS(i)+Bilancio_energia_interna(i);
end

figure
t2=[];
for i=1:length(index)
t2=[t2;(Uscita_CHP(i)-Ingresso_CHP(i))/1000,(Uscita_Solare(i)-Ingresso_Solare(i))/1000,(-Perdite_TS(i))/1000,Bilancio_energia_interna(i)/1000];
end
bar(X,t2,'stacked')
legend('Contributo CHP','Contributo PVT','Perdite termiche','Energia interna accumulo')
grid on
title('Bilancio energia termica casi ottimali')
ylabel('MWh')
xlabel('Scenari ottimi')

%----------------- sostituzione CHP con caldaia della stessa taglia, per tutti i casi ottimi -----------------------

eta_caldaia=0.9;
Q_caldaia_int=[];
for i=1:length(index)
Q_caldaia_int=[Q_caldaia_int,A(2,index(i))];
end
Fattore_conv_en_prim=1.95;
for i=1:length(index)
En_prim_caldaia_integrato(i)=((Uscita_CHP(i)-Ingresso_CHP(i))/eta_caldaia) + (sum(El_en_eff(:,index(i))))/60*Fattore_conv_en_prim;% KWh di energia primaria 
CO2_caldaia_integrato(i)=((Uscita_CHP(i)-Ingresso_CHP(i))/eta_caldaia/PC)*2.75+(sum(El_en_eff(:,index(i))))/60*0.3;
Costi_tot_caldaia_integrato(i)=(1+0.02)*(2222+56*Q_caldaia_int(i))+((Uscita_CHP(i)-Ingresso_CHP(i))/eta_caldaia)*price_gas*year+(sum(El_en_eff(:,index(i))))/60*price_el_acquisto*year;%€
Costi_tot_CHP_integrato(i)=costo_CHP(index(i))+costo_operativo_CHP(index(i))*year+costo_gas(index(i))*year-sum(((El_en_surplus_CHP(:,cont))/60)*price_el_vendita)*year;
end

%----------------- CONFRONTO SULLA ENERGIA PRIMARIA CONSUMATA -----------------------

figure 
hold on
b(1)=bar(1,round(En_prim_caldaia_integrato(1)/1000,1),0.5,'b');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
b(1)=bar(1+0.2,round(sum(Pr_en_CHP(:,index(1)))/60/1000,1),0.5,'r');
xtips2 = b(1).XEndPoints;
ytips2 = b(1).YEndPoints;
labels2 = string(b(1).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
for i=2:length(index)
b(i)=bar(i,round(En_prim_caldaia_integrato(i)/1000,1),0.5,'b');
xtips1 = b(i).XEndPoints;
ytips1 = b(i).YEndPoints;
labels1 = string(b(i).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
b(i)=bar(i+0.2,round(sum(Pr_en_CHP(:,index(i)))/60/1000,1),0.5,'r');
xtips2 = b(i).XEndPoints;
ytips2 = b(i).YEndPoints;
labels2 = string(b(i).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
end
hold off
grid on
title('Confronto impianto integrato con caldaia in sostituzione CHP Energia primaria usata')
xlabel('Scenari ottimi')
legend('Energia primaria caldaia','Energia primaria CHP')
ylabel('MWh')

%----------------------- CONFRONTO SULLA CO2 CONSUMATA -------------------------

figure 
hold on
b(1)=bar(1,round(CO2_caldaia_integrato(1)/1000,2),0.5,'b');
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
b(1)=bar(1+0.2,round(CO2_CHP(index(1))/1000,2),0.5,'r'); 
xtips2 = b(1).XEndPoints;
ytips2 = b(1).YEndPoints;
labels2 = string(b(1).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
for i=2:length(index)
b(i)=bar(i,round(CO2_caldaia_integrato(i)/1000,2),0.5,'b');
xtips1 = b(i).XEndPoints;
ytips1 = b(i).YEndPoints;
labels1 = string(b(i).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
b(i)=bar(i+0.2,round(CO2_CHP(index(i))/1000,2),0.5,'r');
xtips2 = b(i).XEndPoints;
ytips2 = b(i).YEndPoints;
labels2 = string(b(i).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
end
hold off
grid on
title('Confronto impianto integrato con caldaia in sostituzione CHP, CO2 emessa annuale')
xlabel('Scenari ottimi')
legend('CO_2 emessa caldaia','CO_2 emessa CHP')
ylabel('Ton CO_2')

%----------------------- CONFRONTO COSTI TOTALI SUI 20 ANNI---------------

t3=[];t4=[];
for i=1:length(index)
t3=[t3;(2222+56*Q_caldaia_int(i))/1000,(0.02*(2222+56*Q_caldaia_int(i))+((Uscita_CHP(i)-Ingresso_CHP(i))/eta_caldaia)*price_gas*year+(sum(El_en_eff(:,index(i))))/60*price_el_acquisto*year)/1000];
t4=[t4;costo_CHP(index(i))/1000,(costo_operativo_CHP(index(i))*year+costo_gas(index(i))*year-sum(((El_en_surplus_CHP(:,cont))/60)*price_el_vendita)*year)/1000];
end
figure 
x3=[1:1:length(index)];
x4=[1.2:1:length(index)+0.2];
bar(x4,t4,'stacked')
hold on
bar(x3,t3,'stacked')
% b(i)=bar(i+0.2,round(Costi_tot_CHP_integrato(i)/1000,1),0.5,'r');
% xtips2 = b(i).XEndPoints;
% ytips2 = b(i).YEndPoints;
% labels2 = string(b(i).YData);
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% b(i)=bar(i,round(Costi_tot_caldaia_integrato(i)/1000,1),0.5,'b');
% xtips1 = b(i).XEndPoints;
% ytips1 = b(i).YEndPoints;
% labels1 = string(b(i).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% end
hold off
grid on
title('Confronto impianto integrato con caldaia in sostituzione CHP, Costi totali in 20 anni di funzionamento')
xlabel('Scenari ottimi')
legend('Costi investimento CHP','Costi operativi CHP','Costi investimento caldaia','Costi operativi caldaia')
ylabel('K€')

%---------------------- STUDIO CASO OTTIMO, A=284 ------------------------

%------------------- PARTICOLARE GIORNO ENERGIA ELETTRICA COMPRATA---------

T_est_giornaliera=reshape(T_est,[24,365]);T_est_giornaliera_media=sum(T_est_giornaliera)/24;
[T_min,g]=min(T_est_giornaliera_media);
[T_max_g,g]=max(T_est_giornaliera_media);
Energia_comprata_rete=El_en_res(:,index(7));
En_comprata_rete_giorno=reshape(Energia_comprata_rete,[1440,365]);
En_comprata_rete_363=En_comprata_rete_giorno(:,363);
figure
plot(En_comprata_rete_363,'b')
xlabel('Minuti')
ylabel('KWh')
grid on
title('Energia elettrica comprata dalla rete nel giorno con radiazione solare massima e minima')
En_comprata_rete_174=En_comprata_rete_giorno(:,174);
hold on
plot(En_comprata_rete_174,'r')
legend('29 Dicembre','22 Giugno')
hold off

%------------- ANDAMENTO TEMPERATURA ACCUMULO IN UN GIORNO TIPO-----------

figure
En_comprata_rete_giorno2=reshape(Energia_comprata_rete,[1440,365])';
[minuti,giorni]=meshgrid(0:1439,1:365);
surf(minuti,giorni,En_comprata_rete_giorno2,'EdgeColor','none')
xlabel('Minuti')
ylabel('Giorni')
zlabel('KWh')
title('Energia elettrica comprata dalla rete durante l''anno')
colorbar

T_TS_i=T_TS(1:end-1,index(7));
T_TS_giorno_tipo=reshape(T_TS_i,[1440 365]);
T_mandata_acs_giorno=T_mandata_ACS(1:end);
T_ACS_giorni=reshape(T_mandata_acs_giorno,[1440 365]);
int_giorni=int(1:end-1,index(7));
int_matrice=reshape(int_giorni,[1440,365]);
figure
plot(T_TS_giorno_tipo(:,23),'b')
hold on
plot(T_ACS_giorni(:,23),'r')
plot(int_matrice(:,23)*80,'k')
plot([0 1440],[55 55],'-.g')
title('Confronto temperature in un giorno tipo(23 Gennaio)')
xlabel('Minuti')
ylabel('T[°C]')
ylim([0 80])
grid on
legend('Temperatura accumulo','Temperatura richiesta ACS','interruttore CHP on','limite temperatura CHP')
xlim([0 1440])

%------------------ TORTA BILANCIO TERMICO ----------------------

figure
Q_PVT_caso_ott=Uscita_Solare(7)-Ingresso_Solare(7);
Q_CHP_caso_ott=Uscita_CHP(7)-Ingresso_CHP(7);
ax= gca(); 
pieData = [Q_PVT_caso_ott, Q_CHP_caso_ott];
h=pie(ax,pieData);
newColors = [...
    0,       1,       0.49609;   %spring green
    0.8500 0.3250 0.0980];  % arancione chiaro
ax.Colormap =newColors; 
legend('Energia termica dal PVT','Energia termica dal CHP')

E_CHP_caso_ott=(sum(El_en_eff(:,index(7)))-sum(El_en_surplus_CHP(:,index(7))))/sum(Carico_elettrico);
E_PVT_caso_ott=(sum(Q_el_PV(:,index(7)))-sum(El_en_surplus_PVT(:,index(7))))/sum(Carico_elettrico);
E_res_caso_ott=sum(El_en_res(:,index(7)))/sum(Carico_elettrico);

%---------------------- TORTA BILANCIO ELETTRICO--------------------

figure
ax2 = gca(); 
pieData2 = [E_CHP_caso_ott E_PVT_caso_ott E_res_caso_ott]; 
h2 = pie(ax2, pieData2); 
newColors2 = [...
    0.8500 0.3250 0.0980;   % arancione chiaro
    0,       1,       0.49609;   %spring green
    1,1,0];  % giallo
ax2.Colormap =newColors2; 
legend('Energia elettrica dal CHP','Energia elettrica dal PVT','Energia elettrica dalla rete')

