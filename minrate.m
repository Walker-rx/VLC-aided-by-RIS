function min_rate = minrate(B,pt,q,R_pd,N,channel_gain)  
%%%%%%% Intelligent Reflecting Surface-Aided Indoor Visible Light Communication Systems
%%%%%%% formula  (7)  
%%%%%%% B:system bandwidth ; pt:optical transmit power ; 
%%%%%%% q:ratio of the transmitted optical power to the electrical power
%%%%%%% R_pd: responsivity of the PD
%%%%%%% N: power spectral density of noise at the PD
%%%%%%% channel_gain: Not multiplied by pt
Numerator = (R_pd*channel_gain*pt/q)^2; 
Denominator = N*B;
k = exp(1)/(2*pi);
min_rate = B*log2(1+ k*Numerator/Denominator );
end