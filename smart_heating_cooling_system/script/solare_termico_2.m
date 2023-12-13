function [Q_u,Tw_out,eta,k] = solare_termico_2(A_sol_ST,n_ST,Tw_in,mw_in,I_sol,T_ext,eta0,IAM,a1,a2)
% Solare Termico
% Fornendo in ingresso i valori di eta_0, IAM, a1 e a2 (da dati di catalogo), si calcola l'energia termica fornita
% da solare termico in base al numero di pannelli (n_ST) e all'area di ogni pannello (A_sol_ST).
% Deve essere fornita la temperatura di ingresso al pannello (Tw_in) e la portata (mw_in) e si ricava la temperatura di uscita
% e l'efficienza del sistema.
% if minute==1
%     Tw_in(1,cont)=15;
%     A=a2;
%     B=2*mw_in*4.2+a1*A_sol_ST*n_ST-2*a2*A_sol_ST*n_ST*T_ext;
%     C=-2*mw_in*4.2*Tw_in(1,cont)-eta0*IAM*I_sol*A_sol_ST*n_ST-a1*A_sol_ST*n_ST*T_ext+a2*A_sol_ST*n_ST*T_ext^2;
%     Tw_avg=(-B+sqrt(B^2-4*A*C))/(2*A); %[°C]
%     Tw_out=2*Tw_avg-Tw_in(1,cont); %[°C]
%     if Tw_out>80
%         k=0;
%         Q_u =mw_in*4.2*(Tw_out-Tw_in(1,cont)); %[kW]
%         eta=Q_u/(I_sol*A_sol_ST*n_ST); %[-]
%         Tw_out=15;
%     else
%         Q_u =mw_in*4.2*(Tw_out-Tw_in(1,cont)); %[kW]
%         eta=Q_u/(I_sol*A_sol_ST*n_ST); %[-]
%     end
%     
%     if Q_u<0
%         k=0;
%     else if Q_u>=0
%             k=1;
%         end
%     end
% else

    A=a2*A_sol_ST*n_ST;
    B=2*mw_in*4.2+a1*A_sol_ST*n_ST-2*a2*A_sol_ST*n_ST*T_ext;
    C=-2*mw_in*4.2*Tw_in-eta0*IAM*I_sol*A_sol_ST*n_ST-a1*A_sol_ST*n_ST*T_ext+a2*A_sol_ST*n_ST*T_ext^2;
    Tw_avg=(-B+sqrt(B^2-4*A*C))/(2*A); %[°C]
    Tw_out=2*Tw_avg-Tw_in; %[°C]
    if Tw_out>80
        Tw_out=80;
        %k=0;
        Q_u =mw_in*4.2*(Tw_out-Tw_in); %[kW]
        eta=Q_u/(I_sol*A_sol_ST*n_ST); %[-]
%         Tw_out=15;
    else
        Q_u =mw_in*4.2*(Tw_out-Tw_in); %[kW]
        eta=Q_u/(I_sol*A_sol_ST*n_ST); %[-]
    end
    
    if Q_u<0
        k=0;
        Q_u=0;
        Tw_out=Tw_in;
    else %if Q_u>=0
        k=1;
       
    end

end




