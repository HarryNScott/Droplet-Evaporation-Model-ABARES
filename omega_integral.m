function omega_D = omega_integral(T_ref,k_B,eps_v,eps_s)
% Written by Harry Scott
% Date: 18-03-17
% Determine the dimensionless collision integral

T_star = k_B*T_ref/sqrt(eps_v*eps_s);

omega_D = 1.06036/((T_star)^0.1561) + 0.193/exp(0.47635*T_star) + 1.03587/exp(1.52996*T_star);
end
    