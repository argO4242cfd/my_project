function [En_el_PV,eta_effettivo_PV] = fotovoltaico(I_sol,S_PV_unitario, n_PV, eta_1, eta_2, T_ext, T_ext_1, T_ext_2)
% Fotovoltaico
% Il rendimento dipende solo dalla temperatura esterna. Devono essere forniti due valori
% di rendimento del fotovoltaico a due diversi valori di temperatura esterna. 
% In base al rendimento, 
% si calcola l'energia elettrica prodotta
% in base alla radiazione solare globale incidente e la superficie
% totale.

eta_effettivo_PV = eta_2+(eta_1-eta_2)*(T_ext-T_ext_2)/(T_ext_1-T_ext_2); %[-]
En_el_PV = eta_effettivo_PV*S_PV_unitario*n_PV*I_sol; %[kW]

end

