function [T_TS,Tw_out,Q_losses,m_eff_TES] = accumulo_termico_2(m_1,T_in_1,m_2,T_in_2,m_3,T_in_3, K_boll,T_amb,T_TS_0,Vol,T_ACS)
 % Accumulo termico
 % Si suppone perfetto miscelamento (temperatura in tutto il volume
 % dell'acqua= T_TS). L'accumulo, di volume "Vol" in m3, si trova
 % nell'ambiente a temperatura T_amb (esterno, ambiente a temperatura
 % controllata...) e gli scambi termici per trasmissione sono calcolati con
 % il parametro K_boll. Nel calcolo della temperatura di accumulo si
 % considerano anche due ingressi con portate m_1 e m_2 e temperature
 % T_in_1 e T_in_2. 

m_eff_TES=(m_3*(T_ACS-T_in_3))/(T_TS_0-T_in_3);
T_TS=((1000*4.2*Vol/60)*T_TS_0+m_1*4.2*T_in_1+m_2*4.2*T_in_2+(K_boll/1000)*T_amb+m_eff_TES*4.2*T_in_3)/((1000*4.2*Vol/60)+m_eff_TES*4.2+(K_boll/1000)+m_1*4.2+m_2*4.2);

 Tw_out=T_TS; %[°C]
 Q_losses=K_boll*(T_TS-T_amb)/1000; %[kW]
 
end

