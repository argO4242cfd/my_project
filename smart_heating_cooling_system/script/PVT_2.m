function [Th_en_u,Tw_out,eta_th,El_en_PV,eta_eff_PV,k] = PVT_2(A_sol_PVT,n_PVT,Tw_in,mw_in,I_sol,T_ext,eta0,IAM,a1,a2, eta_1, eta_2, T_ext_1, T_ext_2)
% PVT (solare cogenerativo)
% Sono previste due modalità di funzionamento. La prima modalità si
% verifica quando la portata circolante all'interno del pannello è maggiore
% di 0, in tal caso il PVT è in grado di fornire sia energia termica che
% elettrica. Il modello per la stima della producibilità di energia termica
% è analogo a quello del solare termico, mentre per la producibilità da
% fotovoltaico si tiene conto che la temperatura della cella sia uguale a
% quella dell'acqua circolante, e che questa influenzi le prestazioni del
% modulo fotovoltaico. 
% Nel caso in cui non ci sia portata circolante, si ha sola produzione di
% energia elettrica con il modello tradizionale del modulo fotovoltaico.

if I_sol >0 
    if mw_in>0
        
        [Th_en_u,Tw_out,eta_th,k] = solare_termico_2(A_sol_PVT,n_PVT,Tw_in,mw_in,I_sol,T_ext,eta0,IAM,a1,a2);
        
        Tavg_PVT=0.5*(Tw_out+Tw_in);
        
        eta_eff_PV = (eta_1*(1-0.0048*(Tavg_PVT-25))); %[-]
        
        El_en_PV = eta_eff_PV*A_sol_PVT*n_PVT*I_sol; %[kW]
        
    else
        
        [El_en_PV,eta_eff_PV] = fotovoltaico(I_sol, A_sol_PVT, n_PVT, eta_1, eta_2, T_ext, T_ext_1, T_ext_2);
        
        Th_en_u=0;
        
        Tw_out=0;
        
        eta_th=eps;        
        
    end
else
    Th_en_u=0;
    Tw_out=15;
    eta_th=eps;
    El_en_PV=0;
    eta_eff_PV=eps;
    k=0;
end

