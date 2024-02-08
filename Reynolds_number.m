function Re_d = Reynolds_number(r_droplet, U_droplet, U_env, rho_v, mu_v)

% Written by Harry Scott
% Date: 18-03-17
% Determine the reynolds number based on the droplet velocity

	Re_d = 2*r_droplet*abs(U_droplet - U_env)*rho_v/mu_v;
