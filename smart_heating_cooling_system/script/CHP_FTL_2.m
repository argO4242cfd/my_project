function [Th_en_eff, FC_th, T_m_eff, El_en_eff, Pr_en_CHP, eta_th_eff, Th_en_res,El_en_res, El_en_surplus] = CHP_FTL_2(Th_en_nd, El_en_nd, ...
    T_r, m_acqua,eta_th_nom, Q_CHP_nom, E_CHP_nom, Q_CHP_min, eta_th_min,C_el,Q_el_PV)                                                 
% CHP_FTL 
% Function che simula un microcogeneratore di taglia Q_CHP_nom (potenza
% termica) e E_CHP_nom (potenza elettrica). Viene calcolata l'effettiva
% energia termica che può essere fornita, considerando il minimo di
% modulazione (Q_CHP_min) e il rendimento con cui viene fornita
% (eta_th_eff) con legge lineare. L'energia elettrica fornita viene
% calcolata considerando che questa sia proporzionale al fattore di carico
% termico. Si calcola anche l'eventuale surplus o residuo rispetto alle
% richieste elettriche e il residuo sull'energia termica richiesta. 
% Sulla base del rendimento termico calcolato si calcola l'energia primaria
% in ingresso al microcogeneratore.

Th_en_eff=Th_en_nd; %[kW]
FC_th=Th_en_eff/Q_CHP_nom; %[-]
El_en_eff=FC_th*E_CHP_nom; %[kW]
eta_th_eff=eta_th_min+(Th_en_eff-Q_CHP_min)/(Q_CHP_nom-Q_CHP_min)*(eta_th_nom-eta_th_min); %[-]

T_m_eff=T_r+Th_en_eff/(m_acqua*4.2); %[°C]

Th_en_res=max([0,(Th_en_nd-Th_en_eff)]); %[kW] deve essere sempre 0
El_en_res=max([0,(El_en_nd-El_en_eff)]); %[kW] residuo da prendere dalla rete
El_en_surplus=max([0,(El_en_eff-El_en_nd)]); %[kW] surplus da vendere alla rete

if m_acqua==0
    T_m_eff=T_r;
    El_en_res=max([0,(C_el-El_en_eff-Q_el_PV)]); %[kW] residuo da prendere dalla rete
end

Pr_en_CHP=Th_en_eff/eta_th_eff; %[kW]

end


