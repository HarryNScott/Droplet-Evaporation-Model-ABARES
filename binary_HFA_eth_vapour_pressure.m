function [act_co_eth, act_co_HFA, P_mix] = binary_HFA_eth_vapour_pressure(X_eth,Temperature,p_vapour_eth,p_vapour_HFA,propellant_name)

% Activity coeffients and mixture pressure calculated using paper from
% Gavtash 2015. see link:
% https://www.researchgate.net/publication/305043499_SATURATED_VAPOUR_PRESSURE_SVP_MEASUREMENT_OF_ETHANOLHFA_BINARY_MIXTURES/download
%
% NOTE: That this will be incorrect for the other mixtures but what else
% can we do.

if strcmp(propellant_name,'HFA134a')==1||strcmp(propellant_name,'HFA152a')==1||strcmp(propellant_name,'HFA1234zeE')==1
a=-57.73;
b=4.325;
c=0.4202;
d=0.01125;
e=-0.005239;
f=-0.000759;

elseif strcmp(propellant_name,'HFA227ea')==1
a=0.613;
b=-6.784;
c=0.0171;
d=-2.666;
e=0.03957;
f=-7.149e-5;

else
    error('incompatible propellant')
end

P_mix = p_vapour_eth + (p_vapour_HFA - p_vapour_eth)*(1-X_eth).*(1 + a*X_eth + b*X_eth.^2 + c*X_eth*Temperature + d*X_eth.^3 + e*X_eth.^2*Temperature + f*X_eth*Temperature^2);

act_co_HFA = (1 + a*X_eth + b*X_eth.^2 + c*X_eth*Temperature + d*X_eth.^3 + e*X_eth.^2*Temperature + f*X_eth*Temperature^2);

act_co_eth = (1-a) + (a-b)*X_eth + (b-d)*X_eth.^2 + d*X_eth.^3 + ...
    -c*Temperature + (c-e)*X_eth*Temperature + e*X_eth.^2*Temperature + f*X_eth*Temperature^2 - f*Temperature^2;



% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Plot for temperature specific value
% % % % 
% % % % X_eth_temp = 0:0.01:1;
% % % % 
% % % % 
% % % % P_mix_temp = p_vapour_eth + (p_vapour_HFA - p_vapour_eth)*(1-X_eth_temp).*(1 + a*X_eth_temp + b*X_eth_temp.^2 + c*X_eth_temp*T + d*X_eth_temp.^3 + e*X_eth_temp.^2*T + f*X_eth_temp*T^2);
% % % % 
% % % % act_co_HFA_temp = (1 + a*X_eth_temp + b*X_eth_temp.^2 + c*X_eth_temp*T + d*X_eth_temp.^3 + e*X_eth_temp.^2*T + f*X_eth_temp*T^2);
% % % % 
% % % % act_co_eth_temp = (1-a) + (a-b)*X_eth_temp + (b-d)*X_eth_temp.^2 + d*X_eth_temp.^3 + ...
% % % %     -c*T + (c-e)*X_eth_temp*T + e*X_eth_temp.^2*T + f*X_eth_temp*T^2 - f*T^2;
% % % % 
% % % % P_mix_AC_temp = (1-X_eth_temp).*act_co_HFA_temp*p_vapour_HFA + X_eth_temp.*act_co_eth_temp*p_vapour_eth;
% % % % 
% % % % 
% % % % figure
% % % % plot(X_eth_temp,P_mix_temp)
% % % % 
% % % % figure
% % % % plot(X_eth_temp,P_mix_AC_temp)
% % % % 
% % % % figure
% % % % plot(X_eth_temp,P_mix_temp-P_mix_AC_temp)
% % % % 
% % % % figure
% % % % plot(X_eth_temp,act_co_eth_temp)
% % % % xlabel('X_{eth}')
% % % % ylabel('\gamma_{eth}')
% % % % 
% % % % figure
% % % % plot(X_eth_temp,act_co_HFA_temp)
% % % % xlabel('X_{eth}')
% % % % ylabel('\gamma_{HFA}')
% % % % 
% % % % pause
end