% Written by Harry Scott
% Date: 28-06-22
% This code is designed to simulate the evaporation of droplets from an
% initially stationary position of a monodisperse droplet generator

close all
clear all
clc

%% NOTES %%
%
% 1) B_M is zero if it is negative, ie we only have one way mass transfer
%
% 2) The activity coefficients for 152a and 1234ze are based off teh
% correlation for 134a for mixtures with ethanol
%
% 3) All HFA data and ethanol data is taken off NIST unless otherwsie
% stated
% 
% 4) The mass of each species is limited to zero and teh mass fractions are
% limited between 0 and 1

%% GENERIC USER INPUTS %%

% Data locations
% HOME PC
prop_data = 'H:\FINAL PHD BACKUPS\PhD\PhD\Modelling\Evaporation Model\Propellant Data\';
save_folder = 'H:\Data for AFMC\';
% Which propellant is to be used 'HFA134a','HFA152a','HFA227ea' or 'HFA1234zeE'
user_propellant_name = ["HFA152a"] ;

% Which solver is to be used 'ITC' or 'ETC'
user_solver_method = 'ETC';

% Nozzle diameter
user_nozzle_diameter = 150e-6;

% min mass fraction
min_mass_frac = 1e-9;

%% USER LIQUID INPUTS %%

% Initial droplet temperature [K]
user_T_droplet_i = 273.15-45; 

% Initial droplet radius/diameter
Droplet_Diameters = 1*(0.01162E-9*6/pi)^(1/3); 
user_R_droplet_i = Droplet_Diameters/2;

% Intial droplet velocity (+ve)
user_U_droplet_i = 2;

%%% swept paramters %%%
% Liquid mass diffusivity sensitivity (D = D_data * D_liquid_sweep)
user_D_liquid_sweep = 1e3; %[1e3,1e2,1e1,5e0,1e0,5e-1,1e-1];

% Liquid thermal conductivity sensitivity (K = K_data * K_liquid_sweep)
user_K_liquid_sweep = 1;

% Initial liquid ethanol mass fraction
user_Y_eth_l_i_sweep = 0.15; %[0,0.025,0.05,0.10,0.15,0.20,0.50]; 

%% USER VAPOUR INPUTS %%

% Initial counterflow temperature
user_T_counter_i = 273.15+193; 

% Initial CF relative humidity [0-1] (Note: NOT IMPLEMENTED) 
user_RH_i = 0;

% Initial vapour ethanol mass fraction 
user_Y_eth_v_i =  0;

% Initial vapour mass fraction of propellant
user_Y_HFA_v_i = 0;

% Counterflow Velocity (-ve)
user_U_counter_i = -0.864;

%% CURRENT ITERATION INPUT SWEEP %%

% Sweep possible propellants
for sensitivity_sweep_HFA = 1:length(user_propellant_name)

% Sweep possible K (thermal conductivity) values
for sensitivity_sweep_K = 1:length(user_K_liquid_sweep)
    
% Sweep possible D (mass diffusivity) values
for sensitivity_sweep_D = 1:length(user_D_liquid_sweep)
    
% Sweep possible concentrations (ethanol mass fraction) values
for sensitivity_sweep_Y_eth_l_i = 1:length(user_Y_eth_l_i_sweep)
    
    % Reset all input and output variables for next simulation
    close all
    clearvars -except prop_data user_propellant_name user_solver_method...
        user_nozzle_diameter user_T_droplet_i user_R_droplet_i...
        user_D_liquid_sweep user_K_liquid_sweep user_Y_eth_l_i_sweep...
        user_U_droplet_i user_T_counter_i user_RH_i user_Y_eth_v_i...
        user_Y_HFA_v_i user_U_counter_i min_mass_frac...
        sensitivity_sweep_K sensitivity_sweep_D sensitivity_sweep_Y_eth_l_i...
        sensitivity_sweep_HFA save_folder
    clc
    
%% GET UNIVERSAL CONSTANTS %%

% Define the boltzmann constant
k_B = 1.38064852E-23; % SI units (m^2kgs^-2K^-1)

% Liquid Mass Diffusivity (NEEDS TO BE UPDATED)
D_liquid = 1e-8;

% Gravitationl constant
grav_constant = 9.81;

%% GET REFRIGERANT THERMDYNAMIC DATA %%

% Import refrigerant data for the appropriate propellant
if strcmp(user_propellant_name(sensitivity_sweep_HFA),'HFA134a') == 1

    importfile_general(join([prop_data,'R134aData.mat'],''))
    M_HFA = 102.032;
    sigma_HFA = 4.6893; % NOT SI units (A)
    eps_HFA = 299.363*k_B; % SI units
    
    imp_HFA_Cp_liquid = CplJkgK;
    imp_HFA_Cp_vapour = CpvJkgK;
    imp_HFA_density_liquid = Densitylkgm3;
    imp_HFA_density_vapour = Densityvkgm3;
    imp_HFA_pressure = PressurePa;
    imp_HFA_latent_heat = InternalEnergyvJkg - InternalEnergylJkg;
    imp_HFA_temperature = TemperatureK;
    imp_HFA_conductivity_liquid = ThermCondlWmK;
    imp_HFA_conductivity_vapour = ThermCondvWmK;
    imp_HFA_mu_viscosity_liquid = ViscositylPas;
    imp_HFA_mu_viscosity_vapour = ViscosityvPas;
    imp_HFA_surf_ten = SurfTensionlNm;
    clear CplJkgK CpvJkgK CvlJkgK CvvJkgK Densitylkgm3 Densityvkgm3 EnthalpylJkg...
        EnthalpyvJkg EntropylJkgK EntropyvJkgK InternalEnergylJkg...
        InternalEnergyvJkg JouleThomsonlKPa JouleThomsonvKPa PressurePa...
        SoundSpdlms SoundSpdvms SurfTensionlNm TemperatureK ThermCondlWmK...
        ThermCondvWmK ViscositylPas ViscosityvPas Volumelm3kg Volumevm3kg

elseif strcmp(user_propellant_name(sensitivity_sweep_HFA),'HFA227ea') == 1
    
    importfile_general(join([prop_data,'R227eaData.mat'],''))
    M_HFA = 170.03;
    error('URGENT: NEED TO FIND THE LJ PARAMETRS FOR 227\n')
    sigma_HFA = NaN;
    eps_HFA = NaN ;
    
    imp_HFA_Cp_liquid = CplJkgK;
    imp_HFA_Cp_vapour = CpvJkgK;
    imp_HFA_density_liquid = Densitylkgm3;
    imp_HFA_density_vapour = Densityvkgm3;
    imp_HFA_pressure = PressurePa;
    imp_HFA_latent_heat = InternalEnergyvJkg - InternalEnergylJkg;
    imp_HFA_temperature = TemperatureK;
    imp_HFA_conductivity_liquid = ThermCondlWmK;
    imp_HFA_conductivity_vapour = ThermCondvWmK;
    imp_HFA_mu_viscosity_liquid = ViscositylPas;
    imp_HFA_mu_viscosity_vapour = ViscosityvPas;
    imp_HFA_surf_ten = SurfTensionlNm;
    clear CplJkgK CpvJkgK CvlJkgK CvvJkgK Densitylkgm3 Densityvkgm3 EnthalpylJkg...
        EnthalpyvJkg EntropylJkgK EntropyvJkgK InternalEnergylJkg...
        InternalEnergyvJkg JouleThomsonlKPa JouleThomsonvKPa PressurePa...
        SoundSpdlms SoundSpdvms SurfTensionlNm TemperatureK ThermCondlWmK...
        ThermCondvWmK ViscositylPas ViscosityvPas Volumelm3kg Volumevm3kg

elseif strcmp(user_propellant_name(sensitivity_sweep_HFA),'HFA152a') == 1
    
    importfile_general(join([prop_data,'R152aData.mat'],''))
    M_HFA = 66.05;
    sigma_HFA = 3.845; % NOT SI units (A)
    eps_HFA = 119*k_B; % SI units
    
    imp_HFA_Cp_liquid = CplJkgK;
    imp_HFA_Cp_vapour = CpvJkgK;
    imp_HFA_density_liquid = Densitylkgm3;
    imp_HFA_density_vapour = Densityvkgm3;
    imp_HFA_pressure = PressurePa;
    imp_HFA_latent_heat = InternalEnergyvJkg - InternalEnergylJkg;
    imp_HFA_temperature = TemperatureK;
    imp_HFA_conductivity_liquid = ThermCondlWmK;
    imp_HFA_conductivity_vapour = ThermCondvWmK;
    imp_HFA_mu_viscosity_liquid = ViscositylPas;
    imp_HFA_mu_viscosity_vapour = ViscosityvPas;
    imp_HFA_surf_ten = SurfTensionlNm;
    clear CplJkgK CpvJkgK CvlJkgK CvvJkgK Densitylkgm3 Densityvkgm3 EnthalpylJkg...
        EnthalpyvJkg EntropylJkgK EntropyvJkgK InternalEnergylJkg...
        InternalEnergyvJkg JouleThomsonlKPa JouleThomsonvKPa PressurePa...
        SoundSpdlms SoundSpdvms SurfTensionlNm TemperatureK ThermCondlWmK...
        ThermCondvWmK ViscositylPas ViscosityvPas Volumelm3kg Volumevm3kg

elseif strcmp(user_propellant_name(sensitivity_sweep_HFA),'HFA1234zeE') == 1

    importfile_general(join([prop_data,'R1234zeEData.mat'],'')) %refprop
    M_HFA = 114.04;
    sigma_HFA = 5; % NOT SI units (A) 
    eps_HFA = 340*k_B; % SI units %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5103321/table/T2/
    %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5103321/
 
    imp_HFA_Cp_liquid = CplJkgK;
    imp_HFA_Cp_vapour = CpvJkgK;
    imp_HFA_density_liquid = Densitylkgm3;
    imp_HFA_density_vapour = Densityvkgm3;
    imp_HFA_pressure = PressurePa;
    imp_HFA_latent_heat = InternalEnergyvJkg - InternalEnergylJkg;
    imp_HFA_temperature = TemperatureK;
    imp_HFA_conductivity_liquid = ThermCondlWmK;
    imp_HFA_conductivity_vapour = ThermCondvWmK;
    imp_HFA_mu_viscosity_liquid = ViscositylPas;
    imp_HFA_mu_viscosity_vapour = ViscosityvPas;
    imp_HFA_surf_ten = SurfTensionlNm;
    clear CplJkgK CpvJkgK CvlJkgK CvvJkgK Densitylkgm3 Densityvkgm3 EnthalpylJkg...
        EnthalpyvJkg EntropylJkgK EntropyvJkgK InternalEnergylJkg...
        InternalEnergyvJkg JouleThomsonlKPa JouleThomsonvKPa PressurePa...
        SoundSpdlms SoundSpdvms SurfTensionlNm TemperatureK ThermCondlWmK...
        ThermCondvWmK ViscositylPas ViscosityvPas Volumelm3kg Volumevm3kg

else
    error([user_propellant_name(sensitivity_sweep_HFA),' is not an appropriate propellant name'])
end

%% GET ETHANOL THERMDYNAMIC DATA %% 

importfile_general(join([prop_data,'ethanolData.mat'],'')) %refprop
M_eth = 46.07;
sigma_eth = 3.58; % NOT SI units (A)
eps_eth = 507*k_B; % SI units https://onlinelibrary.wiley.com/doi/pdf/10.1002/9783527676750.app5

imp_eth_Cp_liquid = CplJkgK;
imp_eth_Cp_vapour = CpvJkgK;
imp_eth_density_liquid = Densitylkgm3;
imp_eth_density_vapour = Densityvkgm3;
imp_eth_pressure = PressurePa;
imp_eth_latent_heat = InternalEnergyvJkg - InternalEnergylJkg;
imp_eth_temperature = TemperatureK;
imp_eth_conductivity_liquid = ThermCondlWmK;
imp_eth_conductivity_vapour = ThermCondvWmK;
imp_eth_mu_viscosity_liquid = ViscositylPas;
imp_eth_mu_viscosity_vapour = ViscosityvPas;
imp_eth_surf_ten = SurfTensionlNm;
clear CplJkgK CpvJkgK CvlJkgK CvvJkgK Densitylkgm3 Densityvkgm3 EnthalpylJkg...
        EnthalpyvJkg EntropylJkgK EntropyvJkgK InternalEnergylJkg...
        InternalEnergyvJkg JouleThomsonlKPa JouleThomsonvKPa PressurePa...
        SoundSpdlms SoundSpdvms SurfTensionlNm TemperatureK ThermCondlWmK...
        ThermCondvWmK ViscositylPas ViscosityvPas Volumelm3kg Volumevm3kg
    
%% GET AIR THERMDYNAMIC DATA %%
importfile_general(join([prop_data,'airData.mat'],'')) %refprop

M_air = 28.97; % NOT SI units (g/mol)
sigma_air = 3.711; % NOT SI units (A) 
eps_air = 78.6*k_B; % SI units 

imp_air_Cp_vapour = CpvJkgK;
imp_air_density_vapour = Densityvkgm3;
imp_air_pressure = PressurePa;
imp_air_temperature = TemperatureK;
imp_air_conductivity_vapour = ThermCondvWmK;
imp_air_mu_viscosity_vapour = ViscosityvPas;

clear  CpvJkgK CvvJkgK Densityvkgm3 EnthalpyvJkg EntropyvJkgK...
    InternalEnergyvJkg JouleThomsonvKPa PressurePa SoundSpdvms...
    TemperatureK ThermCondvWmK ViscosityvPas Volumevm3kg

% Define ambient pressure in various units
p_amb = 101325; % SI units (pa)
p_amb_bar = p_amb/100000;% (bar)

%% FINALISE LIQUID INPUT INITIAL CONDITIONS FROM SWEEPING %%

% Initial temperature
current_iter_initial_drop_temperature = user_T_droplet_i; % SI units (K)

% Initial radius
current_iter_initial_radius = user_R_droplet_i; % SI units (m)

% Initial velocities
current_iter_initial_droplet_velocity = user_U_droplet_i; % SI units (m/s)

% Initial liquid mass diffusivity
% This will vary depending on the diffusivity which is a function of temperature
% D_liquid_sweep(sensitivity_sweep_D)*D_liquid 

% Initial liquid thermal conductivity
% This will vary depending on the conductivity which is a function of temperature
% K_liquid_sweep(sensitivity_sweep_K)*lambda_liquid_mean 

% Initial liquid ethanol mass fraction
current_iter_initial_Y_eth_l = user_Y_eth_l_i_sweep(sensitivity_sweep_Y_eth_l_i);

%% FINALISE VAPOUR INPUT INITIAL CONDITIONS FROM SWEEPING %%

% Initial counterflow temperature
current_iter_initial_counter_temperature = user_T_counter_i;

% Initial vapour ethanol mass fractions
current_iter_initial_counter_Y_eth = user_Y_eth_v_i;

% Initial vapour propellant mass fractions
current_iter_initial_counter_Y_HFA = user_Y_HFA_v_i;

% Initial vapour air mass fraction 
current_iter_initial_counter_Y_air = 1-current_iter_initial_counter_Y_HFA-current_iter_initial_counter_Y_eth;

% Initial counterflow velocity
current_iter_initial_counter_velocity = user_U_counter_i; % SI units (m/s)

%% VECTORISE LIQUID INITIAL CONDITIONS %%

% Initialise droplet temperatures. Assume homogeneity.
T_droplet_mean(1) = current_iter_initial_drop_temperature; % SI units (K)
T_droplet_centre(1) = T_droplet_mean(1); % SI units (K)
T_droplet_surface(1) = T_droplet_mean(1); % SI units (K)

% Initialise droplet dynamics
r_droplet(1) = current_iter_initial_radius; % SI units (m)
U_droplet(1) = current_iter_initial_droplet_velocity; % SI units (m/s)

% Initialise evaported mass
m_droplet_evap_eth(1) = 0; % SI units (kg)
m_droplet_evap_HFA(1) = 0; % SI units (kg)
m_droplet_evap(1) = 0; % SI units (kg)

% Initialise mass fractions in the droplet
Y_liquid_eth_mean(1) = current_iter_initial_Y_eth_l;
Y_liquid_HFA_mean(1) = 1 - Y_liquid_eth_mean(1);

% Molar fractions in the droplet
X_liquid_eth_mean(1) = (Y_liquid_eth_mean(1)/M_eth)/((Y_liquid_eth_mean(1)/M_eth)+((1-Y_liquid_eth_mean(1))/M_HFA));
X_liquid_HFA_mean(1) = 1 - X_liquid_eth_mean(1);

% Assume initial homogeneity of species at the surface
Y_liquid_eth_surface(1) = Y_liquid_eth_mean(1);
X_liquid_eth_surface(1) = X_liquid_eth_mean(1);
Y_liquid_HFA_surface(1) = Y_liquid_HFA_mean(1);
X_liquid_HFA_surface(1) = X_liquid_HFA_mean(1);

% Assume initial homogeneity of species at the centre
Y_liquid_eth_centre(1) = Y_liquid_eth_mean(1);
X_liquid_eth_centre(1) = X_liquid_eth_mean(1);
Y_liquid_HFA_centre(1) = Y_liquid_HFA_mean(1);
X_liquid_HFA_centre(1) = X_liquid_HFA_mean(1);

% Intialise droplet locarion
x_droplet_position(1) = 0; 
x_droplet_position_nomralised(1) = 0;

% Initialise droplet volume
V_droplet_i = 4/3*pi*r_droplet(1)^3;

% Initialise droplet density
rho_liquid_eth_mean(1) = property_interpolating_function(T_droplet_mean(1),imp_eth_temperature,imp_eth_density_liquid); % SI units (kg/m^3) 
rho_liquid_HFA_mean(1) = property_interpolating_function(T_droplet_mean(1),imp_HFA_temperature,imp_HFA_density_liquid); % SI units (kg/m^3) 
rho_liquid_mean(1) = (Y_liquid_eth_mean(1)/rho_liquid_eth_mean(1) + Y_liquid_HFA_mean(1)/rho_liquid_HFA_mean(1))^-1; % SI units (kg/m^3) 

% Initial droplet mass
m_droplet_eth(1)=V_droplet_i*rho_liquid_mean(1)*Y_liquid_eth_mean(1); % SI units (kg) 
m_droplet_HFA(1)=V_droplet_i*rho_liquid_mean(1)*Y_liquid_HFA_mean(1); % SI units (kg) 
m_droplet(1)=V_droplet_i*rho_liquid_mean(1); % SI units (kg) 

% Secondary droplet mass calculation
m_droplet_eth_vol_int(1) = m_droplet_eth(1); % SI units (kg) 
m_droplet_HFA_vol_int(1) = m_droplet_HFA(1); % SI units (kg) 
m_droplet_vol_int(1) = m_droplet(1); % SI units (kg) 

% Secondary droplet mass calculation
m_R3_droplet_eth_vol_int(1) = m_droplet_eth(1)/r_droplet(1)^3; % SI units (kg) 
m_R3_droplet_HFA_vol_int(1) = m_droplet_HFA(1)/r_droplet(1)^3; % SI units (kg) 
m_R3_droplet_vol_int(1) = m_droplet(1)/r_droplet(1)^3; % SI units (kg) 

% Secondary radius
radius_eth(1) = r_droplet(1);
radius_HFA(1) = r_droplet(1);

%% VECTORISE VAPOUR INITIAL CONDITIONS %%

% Inital counterflow velocity and temperature
U_counter(1) = current_iter_initial_counter_velocity; % SI units (m/s)
T_counter(1) = current_iter_initial_counter_temperature; % SI units (K)

% Initialise mass fractions in the counterflow
Y_vapour_eth_counter(1) = current_iter_initial_counter_Y_HFA;
Y_vapour_HFA_counter(1) = current_iter_initial_counter_Y_HFA; 
Y_vapour_air_counter(1) = current_iter_initial_counter_Y_air;
    
% Initialise molar fractions in the counterflow
X_tot_temporary = (Y_vapour_air_counter(1)/M_air) + (Y_vapour_eth_counter(1)/M_eth) + (Y_vapour_HFA_counter(1)/M_HFA);
X_vapour_eth_counter(1) = (Y_vapour_eth_counter(1)/M_eth)/X_tot_temporary;
X_vapour_HFA_counter(1) = (Y_vapour_HFA_counter(1)/M_HFA)/X_tot_temporary;
X_vapour_air_counter(1) = (Y_vapour_air_counter(1)/M_air)/X_tot_temporary;

%% TEMPORAL STEP SIZING & INTIALISATION %%

% Define a refined temporal spacing for low diffusivities
if user_D_liquid_sweep(sensitivity_sweep_D) < 1
    dt_coarse = 1e-5; % SI units (s)
    dt_refined = 1e-6; % SI units (s)
    dt_refined2 = 1e-7; % SI units (s)
else
    dt_coarse = 1e-4; % SI units (s)
    dt_refined = 1e-5; % SI units (s)
    dt_refined2 = 1e-6; % SI units (s)
end

% Initialise time
time(1) = 0; % SI units (s)
dt(1) = dt_coarse;

%% SPATIAL STEP SIZING %%

% Define approximate number of radial gridpoints
tot_no_gridpoints = 100;

%% RESET COUNTERS %%

% Define an initial counter
ii = 0;
jj = 1;
warn_count = 0;

%% SET TOLERANCING %% 

% maximum allowable difference in interatively calculating the spalding
% heat transfer number
B_T_tol = 1e-7;

% Mass and radius minimum ratio eg. r_frac = r(i)/r_initial
r_ratio_force_quit = 10;
m_ratio_force_quit = r_ratio_force_quit^3;

% Set minimum radius for droplet radius
min_r_droplet = current_iter_initial_radius/r_ratio_force_quit; % SI units (m) 

% Set minimum droplet mass for each component (Not necessarily used)
min_m_droplet_eth = m_droplet_eth(1)/m_ratio_force_quit; % SI units (kg) 
min_m_droplet_HFA = m_droplet_HFA(1)/m_ratio_force_quit;
min_m_droplet = m_droplet(1)/m_ratio_force_quit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVAPORATION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN EVAPORATION LOOP %%

% Evaporate until the requisite cut off quantity is achieved
while r_droplet(jj) > min_r_droplet % OTHER POSSIBLE CONDITIONS && m_droplet_HFA(jj) >= min_m_droplet_HFA && m_droplet_eth(jj) >= min_m_droplet_eth

    % ii > 0 initialises the code skipping the first step with homogeneity
    if ii > 0
        
        % Mass fraction checks
        if Y_liquid_eth_centre(ii) > 1 || Y_liquid_eth_centre(ii) < 0
            fprintf('Y eth (centre) is %f',Y_liquid_eth_centre(ii),'\')
        end
        
        if Y_liquid_eth_mean(ii) > 1 || Y_liquid_eth_mean(ii) < 0
            fprintf('Y eth (mean) is %f',Y_liquid_eth_mean(ii),'\')
        end
        
        if Y_liquid_eth_surface(ii) > 1 || Y_liquid_eth_surface(ii) < 0
            fprintf('Y eth (surface) is %f',Y_liquid_eth_surface(ii),'\')
        end
        
        if Y_liquid_HFA_centre(ii) > 1 || Y_liquid_HFA_centre(ii) < 0
            fprintf('Y HFA (centre) is %f',Y_liquid_HFA_centre(ii),'\')
        end
        
        if Y_liquid_HFA_mean(ii) > 1 || Y_liquid_HFA_mean(ii) < 0
            fprintf('Y HFA (mean) is %f',Y_liquid_HFA_mean(ii),'\')
        end
        
        if Y_liquid_HFA_surface(ii) > 1 || Y_liquid_HFA_surface(ii) < 0
            fprintf('Y HFA (surface) is %f',Y_liquid_HFA_surface(ii),'\')
        end
        
        if Y_vapour_HFA_surface(ii) > 1 || Y_vapour_HFA_surface(ii) < 0
            fprintf('Y HFA Vapour (surface) is %f',Y_vapour_HFA_surface(ii),'\')
        end
        
        if Y_vapour_eth_surface(ii) > 1 || Y_vapour_eth_surface(ii) < 0
            fprintf('Y ETH Vapour (surface) is %f',Y_vapour_eth_surface(ii),'\')
        end
        
        %% SET TIMESTEP %%
        
        % Refine timestepping based on teh gradient of f_eth
        if ii > 1
            if abs(diff(f_eth(ii-1:ii))./diff(time(ii-1:ii))) > 1
                dt(ii) = dt_refined;
            elseif abs(diff(f_eth(ii-1:ii))./diff(time(ii-1:ii))) > 0.01
                dt(ii) = dt_refined;
            else
                dt(ii) = dt_coarse;
            end
            
            % Slowly increase timestepping
            if dt(ii)> dt(ii-1)
                dt(ii) = dt(ii-1) + dt(ii)/1000;
            end
        end
        
        % Propagate forward in time by dt
        time(ii+1) = time(ii) + dt(ii);
        
        %% FALLING DROPLET DYNAMICS %%  
   
        % Force balance with downward as postive motion and euler timestep
        
        % Get droplet drag coefficient, https://www.scielo.br/pdf/bjce/v35n2/1678-4383-bjce-35-02-0709.pdf
        % Limit drag at low Reynolds number
        if Re_d(ii) <= 0.025
            C_D(ii) = 24/0.025 + 3.72/0.025^(0.5) - (4.83e-3*0.025^0.5)/(1 + 3e-6*0.025^1.5) + 0.49;
        else
            C_D(ii) = 24/Re_d(ii) + 3.72/Re_d(ii)^(0.5) - (4.83e-3*Re_d(ii)^0.5)/(1 + 3e-6*Re_d(ii)^1.5) + 0.49;
        end
        
        % Determine droplet velocity from force balance (down is +ve)
        U_droplet(ii+1) = U_droplet(ii) + grav_constant*dt(ii) - (3/8)*C_D(ii)*(rho_vapour_film(ii)/rho_liquid_mean(ii))*((U_droplet(ii)-U_counter(ii))^2/r_droplet(ii))*dt(ii);

        % Get counterflow velocity
        U_counter(ii+1) = U_counter(1);
        
        % Error Check as counterflow must oppose motion (must be -ve)
        if U_counter(ii+1)>0
            error('Counterflow velocity must be -ve')
        end
        
        % Determine dorplet locaion
        x_droplet_position(ii+1) = x_droplet_position(ii)+U_droplet(ii)*dt(ii);
        dx(ii) = x_droplet_position(ii+1) - x_droplet_position(ii);
        x_droplet_position_nomralised(ii+1) = x_droplet_position(ii+1)/user_nozzle_diameter;
        
        % Counterflow composition constant
        Y_vapour_eth_counter(ii+1) = Y_vapour_eth_counter(1);
        Y_vapour_HFA_counter(ii+1) = Y_vapour_HFA_counter(1);
        Y_vapour_air_counter(ii+1) = Y_vapour_air_counter(1);
        
        % Counterflow temperature fixed as a constant
        T_counter(ii+1) = T_counter(1);
        
        %% EVAPORATION RATE & INTIIAL MASS CALCULATION %%
        
        % -ve evaporation is outwards mass transfer
        % Find the evaporation rates of each component
        m_droplet_evap_eth_rate(ii+1) = f_eth(ii)*(-2)*pi*r_droplet(ii)*rho_vapour_film(ii)*D_vapour_film(ii)*B_M(ii)*Sh(ii);
        m_droplet_evap_HFA_rate(ii+1) = f_HFA(ii)*(-2)*pi*r_droplet(ii)*rho_vapour_film(ii)*D_vapour_film(ii)*B_M(ii)*Sh(ii);
        m_droplet_evap_rate(ii+1) = (-2)*pi*r_droplet(ii)*rho_vapour_film(ii)*D_vapour_film(ii)*B_M(ii)*Sh(ii);
        
        % Calculate the mass evaporation flux
        m_droplet_evap_eth_flux(ii+1) = f_eth(ii)*(-1/2)*rho_vapour_film(ii)*D_vapour_film(ii)*B_M(ii)*Sh(ii)/r_droplet(ii);
        m_droplet_evap_HFA_flux(ii+1) = f_HFA(ii)*(-1/2)*rho_vapour_film(ii)*D_vapour_film(ii)*B_M(ii)*Sh(ii)/r_droplet(ii);
        m_droplet_evap_flux(ii+1) = (-1/2)*rho_vapour_film(ii)*D_vapour_film(ii)*B_M(ii)*Sh(ii)/r_droplet(ii);
        
        % Determine the mass evaporated this step
        m_droplet_step_eth_evap(ii+1) = m_droplet_evap_eth_rate(ii+1)*dt(ii);
        m_droplet_step_HFA_evap(ii+1) = m_droplet_evap_HFA_rate(ii+1)*dt(ii);
        m_droplet_step_evap(ii+1) = m_droplet_evap_rate(ii+1)*dt(ii);
        
        % Determine the mass evaporated up to this point
        m_droplet_tot_eth_evap(ii+1) = sum(m_droplet_step_eth_evap);
        m_droplet_tot_HFA_evap(ii+1) = sum(m_droplet_step_HFA_evap);
        m_droplet_tot_evap(ii+1) = sum(m_droplet_step_evap);
        
        % Determine the remaining mass
        m_droplet_eth(ii+1) = m_droplet_eth(1) + m_droplet_tot_eth_evap(ii+1);
        m_droplet_eth(ii+1) = max([0,m_droplet_eth(ii+1)]);
        m_droplet_HFA(ii+1) = m_droplet_HFA(1) + m_droplet_tot_HFA_evap(ii+1);
        m_droplet_HFA(ii+1) = max([0,m_droplet_HFA(ii+1)]);
        m_droplet(ii+1) = m_droplet(1) + m_droplet_tot_evap(ii+1);
        m_droplet(ii+1) = max([0,m_droplet(ii+1)]);
        
        % Calculate the mean mass and molar fractions in the droplet
        Y_liquid_HFA_mean(ii+1) = max([0,min([1,m_droplet_HFA(ii+1)/(m_droplet_HFA(ii+1) + m_droplet_eth(ii+1))])]);
        Y_liquid_eth_mean(ii+1) = max([0,min([1,m_droplet_eth(ii+1)/(m_droplet_HFA(ii+1) + m_droplet_eth(ii+1))])]);
        X_liquid_HFA_mean(ii+1) = (Y_liquid_HFA_mean(ii+1)/M_HFA)/(Y_liquid_eth_mean(ii+1)/M_eth + Y_liquid_HFA_mean(ii+1)/M_HFA);
        X_liquid_eth_mean(ii+1) = (Y_liquid_eth_mean(ii+1)/M_eth)/(Y_liquid_eth_mean(ii+1)/M_eth + Y_liquid_HFA_mean(ii+1)/M_HFA);
        
        % Calculate new droplet radius based off the mean density.
        rho_liquid_mean(ii+1) = ((Y_liquid_eth_mean(ii+1)/rho_liquid_eth_mean(ii)) + (Y_liquid_HFA_mean(ii+1))/rho_liquid_HFA_mean(ii))^-1;
        r_droplet(ii+1) = ((3/(4*pi))*(m_droplet(ii+1)/rho_liquid_mean(ii+1)))^(1/3);
        
        % Calculate the mean new temperature of the droplet (Q = m*c_P*dT/dt)
        % with updated mass. forward time (Euler). Uniform temprture.
        T_droplet_mean(ii+1) = T_droplet_mean(ii) + (dt(ii)/(m_droplet(ii)*c_p_liquid_mean(ii)))*...
            (2*pi*r_droplet(ii)*lambda_vapour_film(ii)*(T_counter(ii)-T_droplet_surface(ii))*Nu_d(ii) + m_droplet_evap_HFA_rate(ii+1)*L_HFA(ii)+ m_droplet_evap_eth_rate(ii+1)*L_eth(ii));
        
        %% ITC MODEL EVAPORATION %%
        
        % The infinite conductivity model. This model assumes
        % instantaneous heat and mass transfer throughout the droplet.
        
        % Check whether the ITC solver is selected
        if strcmp(user_solver_method,'ITC') == 1
            ITC = 1;
            ETC = 0;
                        
            % ITC stipulates radial homogeneity
            T_droplet_surface(ii+1) = T_droplet_mean(ii+1);
            
            % Define the molar fractions at the liquid surface of the droplet
            % based off radial homogeneity
            Y_liquid_HFA_surface(ii+1) = Y_liquid_HFA_mean(ii+1);
            Y_liquid_eth_surface(ii+1) = Y_liquid_eth_mean(ii+1);
            Y_liquid_HFA_centre(ii+1) = Y_liquid_HFA_mean(ii+1);
            Y_liquid_eth_centre(ii+1) = Y_liquid_eth_mean(ii+1);
            X_liquid_HFA_surface(ii+1) = X_liquid_HFA_mean(ii+1);
            X_liquid_eth_surface(ii+1) = X_liquid_eth_mean(ii+1);

        end
        
        %% ETC MODEL EVAPORATION %%
        
        % The effective conductivity model. This model solves the radial
        % diffusion equation with modified coefficents.
        
        % Check whether the ETC solver is selected
        if strcmp(user_solver_method,'ETC') == 1
            ITC = 0;
            ETC = 1;
            
            %% SETUP THE RADIAL GRID %%
            
            % Define the radial matrix for the new radius
            r_location = linspace(0,r_droplet(ii+1),length(vec));
            r_location = r_location(point);
            
            %% DIMENSIONLESS NUMBERS & THERMODYNAMIC PROPERTIES FOR ETC %%
            
            % Calculate the heat diffusion coefficient
            lambda_liquid_surface(ii) = user_K_liquid_sweep(sensitivity_sweep_K)*lambda_liquid_surface(ii);
            lambda_liquid_mean(ii) = user_K_liquid_sweep(sensitivity_sweep_K)*lambda_liquid_mean(ii);
            alpha_liquid_surface(ii) = lambda_liquid_surface(ii)/(c_p_liquid_surface(ii)*rho_liquid_surface(ii));
            alpha_liquid_mean(ii) = lambda_liquid_mean(ii)/(c_p_liquid_mean(ii)*rho_liquid_mean(ii));
            
            % Calculate mass diffusivity. Lambda already modifed.
            D_liquid_mean(ii) = user_D_liquid_sweep(sensitivity_sweep_D)*D_liquid;

            % Determine surafce velocity with Eueler timestepping.
            % (Reynolds number based off diameter) (Sirignano 2010 pp. 56)
            % This is the droplet reynold snumber in the gas flowe
            U_liquid_surface(ii) = (1/(6*pi)*mu_vapour_film(ii)/mu_liquid_mean(ii)*abs(U_counter(ii)-U_droplet(ii))*12.96/(1+B_M(ii))*(2*r_droplet(ii+1)*rho_liquid_mean(ii)/mu_liquid_mean(ii))^(1/3))^(3/2);
            
            % Determine internal droplet dimensionless numbers
            Re_liquid(ii) = rho_liquid_mean(ii)*U_liquid_surface(ii)*2*r_droplet(ii+1)/mu_liquid_mean(ii);
            Sc_liquid(ii) = nu_liquid_mean(ii)/D_liquid_mean(ii);
            Pr_liquid(ii) = c_p_liquid_mean(ii)*mu_liquid_surface(ii)/lambda_liquid_mean(ii);
            
            % Calculate the effective transfer coefficients (Abramzon and
            % Sirignano 1989)
            ETC_mass_coeff(ii) = 1.86 + 0.86*tanh(2.225*log10(Re_liquid(ii)*Sc_liquid(ii)/30));
            ETC_heat_coeff(ii) = 1.86 + 0.86*tanh(2.225*log10(Re_liquid(ii)*Pr_liquid(ii)/30));
            
            
            %% DE VARIABLES %%
            
            % Normalise radius
            zeta = r_location/r_droplet(ii+1);
            dzeta = diff(zeta);

            % Normalised time
            tor(ii+1) = alpha_liquid_0*time(ii+1)/(R_0)^2;
            dtor(ii+1) = tor(ii+1) - tor(ii);
            
            % Normalise outer radius
            r_s(ii+1) = r_droplet(ii+1)/R_0;
            
            % Find Beta
            beta(ii+1) = 1/2*((r_s(ii+1))^2-(r_s(ii))^2)/dtor(ii+1);
            
            %% SETUP LHS MATRIX %%
            
            % We define a DE IS X_MAT(ii+1)*X_LHS(ii+1) = X_RHS(ii+1)
            
            % Preallocate
            DE_Z_RHS(:,ii+1) = zeros(length(zeta),1);
            DE_Y_eth_RHS(:,ii+1) = zeros(length(zeta),1);
            DE_Y_HFA_RHS(:,ii+1) = zeros(length(zeta),1);
            
            DE_Z_LHS(:,ii+1) = zeros(length(zeta),1);
            DE_Y_eth_LHS(:,ii+1) = zeros(length(zeta),1);
            DE_Y_HFA_LHS(:,ii+1) = zeros(length(zeta),1);
            
            DE_Z_MAT = zeros(length(zeta),length(zeta));
            DE_Y_eth_MAT = zeros(length(zeta),length(zeta));
            DE_Y_HFA_MAT = zeros(length(zeta),length(zeta));
            
            % Centre BC's
            DE_Z_MAT(1,1:2) = [-1,1];
            DE_Y_eth_MAT(1,1:2) = [-1,1];
            DE_Y_HFA_MAT(1,1:2) = [-1,1];
            
            % Diffusion equations
            for nn = 2:length(zeta)-1
                
                % Get spacing of zeta on teh left and right
                h = dzeta(nn);
                H = dzeta(nn-1);
                
                % HEAT MATRIX
              	DE_Z_MAT(nn,nn-1) = -(h^2/(h*H*(h+H)))*(alpha_liquid_0*beta(ii+1)*zeta(nn)/(alpha_liquid_mean(ii)))...
                    + (2*h/(h*H*(h+H)))*ETC_heat_coeff(ii) ...
                    - (2*h^2/(h*H*(h+H)))*(ETC_heat_coeff(ii)/zeta(nn));
                
                DE_Z_MAT(nn,nn) = -alpha_liquid_0*((r_s(ii+1))^2)/(alpha_liquid_mean(ii)*dtor(ii+1))...
                    - (H^2-h^2)/(h*H*(h+H))*(alpha_liquid_0*beta(ii+1)*zeta(nn)/(alpha_liquid_mean(ii)))...
                    - 2*(H^2-h^2)/(h*H*(h+H))*(ETC_heat_coeff(ii)/zeta(nn))...
                    - 2*(h+H)/(h*H*(h+H))*ETC_heat_coeff(ii);
                
                DE_Z_MAT(nn,nn+1) = (H^2/(h*H*(h+H)))*(alpha_liquid_0*beta(ii+1)*zeta(nn)/(alpha_liquid_mean(ii)))...
                    + (2*H/(h*H*(h+H)))*ETC_heat_coeff(ii) ...
                    + (2*H^2/(h*H*(h+H)))*(ETC_heat_coeff(ii)/zeta(nn));
                
                % ETH MATRIX                
                DE_Y_eth_MAT(nn,nn-1) = -(h^2/(h*H*(h+H)))*(alpha_liquid_0*beta(ii+1)*zeta(nn)/(D_liquid_mean(ii)))...
                    + (2*h/(h*H*(h+H)))*ETC_mass_coeff(ii) ...
                    - (2*h^2/(h*H*(h+H)))*(ETC_mass_coeff(ii)/zeta(nn));
                
                DE_Y_eth_MAT(nn,nn) = -alpha_liquid_0*((r_s(ii+1))^2)/(D_liquid_mean(ii)*dtor(ii+1))...
                    - (H^2-h^2)/(h*H*(h+H))*(alpha_liquid_0*beta(ii+1)*zeta(nn)/(D_liquid_mean(ii)))...
                    - 2*(H^2-h^2)/(h*H*(h+H))*(ETC_mass_coeff(ii)/zeta(nn))...
                    - 2*(h+H)/(h*H*(h+H))*ETC_mass_coeff(ii);
                
                DE_Y_eth_MAT(nn,nn+1) = (H^2/(h*H*(h+H)))*(alpha_liquid_0*beta(ii+1)*zeta(nn)/(D_liquid_mean(ii)))...
                    + (2*H/(h*H*(h+H)))*ETC_mass_coeff(ii) ...
                    + (2*H^2/(h*H*(h+H)))*(ETC_mass_coeff(ii)/zeta(nn));
                
                % HFA MATRIX
                DE_Y_HFA_MAT(nn,nn-1) = -(h^2/(h*H*(h+H)))*(alpha_liquid_0*beta(ii+1)*zeta(nn)/(D_liquid_mean(ii)))...
                    + (2*h/(h*H*(h+H)))*ETC_mass_coeff(ii) ...
                    - (2*h^2/(h*H*(h+H)))*(ETC_mass_coeff(ii)/zeta(nn));
                
                DE_Y_HFA_MAT(nn,nn) = -alpha_liquid_0*((r_s(ii+1))^2)/(D_liquid_mean(ii)*dtor(ii+1))...
                    - (H^2-h^2)/(h*H*(h+H))*(alpha_liquid_0*beta(ii+1)*zeta(nn)/(D_liquid_mean(ii)))...
                    - 2*(H^2-h^2)/(h*H*(h+H))*(ETC_mass_coeff(ii)/zeta(nn))...
                    - 2*(h+H)/(h*H*(h+H))*ETC_mass_coeff(ii);
                
                DE_Y_HFA_MAT(nn,nn+1) = (H^2/(h*H*(h+H)))*(alpha_liquid_0*beta(ii+1)*zeta(nn)/(D_liquid_mean(ii)))...
                    + (2*H/(h*H*(h+H)))*ETC_mass_coeff(ii) ...
                    + (2*H^2/(h*H*(h+H)))*(ETC_mass_coeff(ii)/zeta(nn));
                
            end
            
            % Surafce BC
            DE_Z_MAT(end,end-1:end) = [1/dzeta(end), -1/dzeta(end) - (lambda_vapour_film(ii)/(ETC_heat_coeff(ii)*lambda_liquid_mean(ii)))*Nu_d(ii)/2];
            DE_Y_eth_MAT(end,end-1:end) = [1/dzeta(end), -1/dzeta(end) - r_droplet(ii+1)*m_droplet_evap_flux(ii+1)/(ETC_mass_coeff(ii)*D_liquid_mean(ii)*rho_liquid_mean(ii))];
            DE_Y_HFA_MAT(end,end-1:end) = [1/dzeta(end), -1/dzeta(end) - r_droplet(ii+1)*m_droplet_evap_flux(ii+1)/(ETC_mass_coeff(ii)*D_liquid_mean(ii)*rho_liquid_mean(ii))];
            
            
            %% SETUP RHS MATRIX %%
            
            % Centre BC
            DE_Z_RHS(1,ii+1) = 0;
            DE_Y_eth_RHS(1,ii+1) = 0;
            DE_Y_HFA_RHS(1,ii+1) = 0;
            
            % DE
            for nn = 2:length(zeta)-1
                DE_Z_RHS(nn,ii+1) = -alpha_liquid_0*((r_s(ii+1))^2)/(alpha_liquid_mean(ii)*dtor(ii+1))*DE_Z_LHS(nn,ii);
                DE_Y_eth_RHS(nn,ii+1) = -alpha_liquid_0*((r_s(ii+1))^2)/(D_liquid_mean(ii)*dtor(ii+1))*DE_Y_eth_LHS(nn,ii);
                DE_Y_HFA_RHS(nn,ii+1) = -alpha_liquid_0*((r_s(ii+1))^2)/(D_liquid_mean(ii)*dtor(ii+1))*DE_Y_HFA_LHS(nn,ii);
            end
            
            % Surface BC
            DE_Z_RHS(end,ii+1) = - (lambda_vapour_film(ii)/(ETC_heat_coeff(ii)*lambda_liquid_mean(ii)))*Nu_d(ii)*(T_counter(ii) - T_liquid_0)/(2*T_liquid_0)...
                - m_droplet_evap_eth_rate(ii+1)*L_eth(ii)/(4*pi*r_droplet(ii+1)*ETC_heat_coeff(ii)*lambda_liquid_mean(ii)*T_liquid_0)...
                - m_droplet_evap_HFA_rate(ii+1)*L_HFA(ii)/(4*pi*r_droplet(ii+1)*ETC_heat_coeff(ii)*lambda_liquid_mean(ii)*T_liquid_0);
            DE_Y_eth_RHS(end,ii+1) = - r_droplet(ii+1)*m_droplet_evap_flux(ii+1)/(ETC_mass_coeff(ii)*D_liquid_mean(ii)*rho_liquid_mean(ii))*f_eth(ii);
            DE_Y_HFA_RHS(end,ii+1) = - r_droplet(ii+1)*m_droplet_evap_flux(ii+1)/(ETC_mass_coeff(ii)*D_liquid_mean(ii)*rho_liquid_mean(ii))*f_HFA(ii);
            
            %% SOLVE DE %%
            
            DE_Z_LHS(:,ii+1) =  DE_Z_MAT\DE_Z_RHS(:,ii+1);
            DE_Y_eth_LHS(:,ii+1) =  DE_Y_eth_MAT\DE_Y_eth_RHS(:,ii+1);
            DE_Y_HFA_LHS(:,ii+1) =  DE_Y_HFA_MAT\DE_Y_HFA_RHS(:,ii+1);
            
            %% CORRECTIONS & ERROR CALCULATIONS %%
            
            % Estimate the error in the inverse process. basically our
            % inverted matrix shoudl give the same result back
            ETC_error_heat_intitial(:,ii+1) = DE_Z_MAT*DE_Z_LHS(:,ii+1)-DE_Z_RHS(:,ii+1);
            ETC_error_eth_initial(:,ii+1) = DE_Y_eth_MAT*DE_Y_eth_LHS(:,ii+1)-DE_Y_eth_RHS(:,ii+1);
            ETC_error_HFA_initial(:,ii+1) = DE_Y_HFA_MAT*DE_Y_HFA_LHS(:,ii+1)-DE_Y_HFA_RHS(:,ii+1);
            
            % Get percent error
            ETC_percent_error_heat_intitial(:,ii+1) = ETC_error_heat_intitial(:,ii+1)./abs(DE_Z_LHS(:,ii+1)-DE_Z_LHS(:,ii));
            ETC_percent_error_eth_initial(:,ii+1) = ETC_error_eth_initial(:,ii+1)./abs(DE_Y_eth_LHS(:,ii+1)-DE_Y_eth_LHS(:,ii));
            ETC_percent_error_HFA_initial(:,ii+1) = ETC_error_HFA_initial(:,ii+1)./abs(DE_Y_HFA_LHS(:,ii+1)-DE_Y_HFA_LHS(:,ii));
                    
            % Get temperature from normalised value
            T_droplet_full(:,ii+1) = DE_Z_LHS(:,ii+1)*T_liquid_0 + T_liquid_0;
            
            % Correct all of the mass fractions
            for cor = 1:size(DE_Y_HFA_LHS,1)
                
                % Force the massfractions to be between 0 and 1
                DE_Y_eth_LHS(cor,ii+1) = max([min_mass_frac,min([1-min_mass_frac,DE_Y_eth_LHS(cor,ii+1)])]);
                DE_Y_HFA_LHS(cor,ii+1) = max([min_mass_frac,min([1-min_mass_frac,DE_Y_HFA_LHS(cor,ii+1)])]);
                
                % Re normalise values to account for errors so that mass fractions
                % must add to 1
                DE_Y_eth_LHS(cor,ii+1) = DE_Y_eth_LHS(cor,ii+1)/(DE_Y_eth_LHS(cor,ii+1)+DE_Y_HFA_LHS(cor,ii+1));
                DE_Y_HFA_LHS(cor,ii+1) = DE_Y_HFA_LHS(cor,ii+1)/(DE_Y_eth_LHS(cor,ii+1)+DE_Y_HFA_LHS(cor,ii+1));

            end
            
            % Check error after corrections
            if ii>1
                ETC_error_heat_final(:,ii+1) = DE_Z_MAT*DE_Z_LHS(:,ii+1)-DE_Z_RHS(:,ii+1);
                ETC_error_eth_initial(:,ii+1) = DE_Y_eth_MAT*DE_Y_eth_LHS(:,ii+1)-DE_Y_eth_RHS(:,ii+1);
                ETC_error_HFA_initial(:,ii+1) = DE_Y_HFA_MAT*DE_Y_HFA_LHS(:,ii+1)-DE_Y_HFA_RHS(:,ii+1);
            end
            

            %% ETC MODEL OUTPUTS %%
            
            % Temperature of the droplet at the surface and centre
            T_droplet_surface(ii+1) = T_droplet_full(end,ii+1);
            T_droplet_centre(ii+1) = T_droplet_full(1,ii+1);
            
            % Mass and molar fractions at the surface
            Y_liquid_eth_surface(ii+1) =  DE_Y_eth_LHS(end,ii+1);
            Y_liquid_HFA_surface(ii+1) =  DE_Y_HFA_LHS(end,ii+1);
            X_liquid_eth_surface(ii+1) = (Y_liquid_eth_surface(ii+1)/M_eth)/((Y_liquid_HFA_surface(ii+1)/M_HFA) + (Y_liquid_eth_surface(ii+1)/M_eth));
            X_liquid_HFA_surface(ii+1) = (Y_liquid_HFA_surface(ii+1)/M_HFA)/((Y_liquid_HFA_surface(ii+1)/M_HFA) + (Y_liquid_eth_surface(ii+1)/M_eth));
            
            % Mass and molar fractions at the centre
            Y_liquid_eth_centre(ii+1) =  DE_Y_eth_LHS(1,ii+1);
            Y_liquid_HFA_centre(ii+1) =  DE_Y_HFA_LHS(1,ii+1);
            X_liquid_eth_centre(ii+1) = (Y_liquid_eth_centre(ii+1)/M_eth)/((Y_liquid_HFA_centre(ii+1)/M_HFA) + (Y_liquid_eth_centre(ii+1)/M_eth));
            X_liquid_HFA_centre(ii+1) = (Y_liquid_HFA_centre(ii+1)/M_HFA)/((Y_liquid_HFA_centre(ii+1)/M_HFA) + (Y_liquid_eth_centre(ii+1)/M_eth));
            
        end
        %% SOLVER ERROR CHECK %%
        
        % Update jj to check if radius or mass has become suffienty small to exit loop
        jj = ii;
        
        % Check appropriate solver selcted
        if strcmp(user_solver_method,'ETC') == 0  && strcmp(user_solver_method,'ITC') == 0
            ITC = 0;
            ETC = 0;
            error('Invalid solver method')
        end
    end

    %% SURFACE PRESSURES AND MIXTURE FRACTIONS %%
    
    % Get boundary layer temperature. Use the 1/3 rule
    T_BL_v(ii+1) = (T_counter(ii+1) + 2*T_droplet_mean(ii+1))/3;
    
    % Determine counterflow molar fractions
    X_vapour_eth_counter(ii+1) = (Y_vapour_eth_counter(ii+1)/M_eth)/((Y_vapour_air_counter(ii+1)/M_air) + (Y_vapour_HFA_counter(ii+1)/M_HFA) + (Y_vapour_eth_counter(ii+1)/M_eth));
    X_vapour_HFA_counter(ii+1) = (Y_vapour_HFA_counter(ii+1)/M_HFA)/((Y_vapour_air_counter(ii+1)/M_air) + (Y_vapour_HFA_counter(ii+1)/M_HFA) + (Y_vapour_eth_counter(ii+1)/M_eth));
    X_vapour_air_counter(ii+1) = (Y_vapour_air_counter(ii+1)/M_air)/((Y_vapour_air_counter(ii+1)/M_air) + (Y_vapour_HFA_counter(ii+1)/M_HFA) + (Y_vapour_eth_counter(ii+1)/M_eth));
    
    % Find molecular mass of vapour in the total spray
    M_vapour_counter(ii+1) = M_air*X_vapour_air_counter(ii+1) + M_HFA*X_vapour_HFA_counter(ii+1) + M_eth*X_vapour_eth_counter(ii+1); 
    
    % Determine vapour pressures at new temperature.
    p_vapour_eth_surface(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_eth_temperature,imp_eth_pressure);
    p_vapour_HFA_surface(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_HFA_temperature,imp_HFA_pressure);
    
    % Use  empirical relation by Gavtash to get activitty coefficients and
    % partial pressures (2015). NOTE USE CORRELATION FOR ALL MIXTURES 
    [act_co_eth(ii+1), act_co_HFA(ii+1),P_mix(ii+1)] = binary_HFA_eth_vapour_pressure(X_liquid_eth_surface(ii+1),T_droplet_surface(ii+1),p_vapour_eth_surface(ii+1),p_vapour_HFA_surface(ii+1),user_propellant_name(sensitivity_sweep_HFA));
        
    % Determine the partial pressure at the surface using modified
    % Raoults law and daltons law
    P_partial_eth(ii+1) = X_liquid_eth_surface(ii+1)*act_co_eth(ii+1)*p_vapour_eth_surface(ii+1);
    P_partial_HFA(ii+1) = X_liquid_HFA_surface(ii+1)*act_co_HFA(ii+1)*p_vapour_HFA_surface(ii+1);
  	
    % Determine the surface vapour fractions
    X_vapour_eth_surface(ii+1) = P_partial_eth(ii+1)/p_amb;
    X_vapour_HFA_surface(ii+1) = P_partial_HFA(ii+1)/p_amb;
    
    % Determine the molar fractions at the surface if the partial
    % pressures are too high (higher than ambient) make a correction
    if P_partial_HFA(ii+1)+P_partial_eth(ii+1) > p_amb
        warn_count = warn_count + 1;
        warning(join(['Partical pressure too high at surface. Count number',num2str(warn_count),'\n'],''));
        X_vapour_eth_surface(ii+1) = X_vapour_eth_surface(ii+1)/(X_vapour_HFA_surface(ii+1)+X_vapour_eth_surface(ii+1));
        X_vapour_HFA_surface(ii+1) = X_vapour_HFA_surface(ii+1)/(X_vapour_HFA_surface(ii+1)+X_vapour_eth_surface(ii+1));
    end
    
    % Determine the mass fractoin of counterflow at the surface (remaining
    % compomnent)
    X_vapour_counter_surface(ii+1) =  1 - X_vapour_eth_surface(ii+1) - X_vapour_HFA_surface(ii+1);
    
    % Part of the vapour at the surface will be a result of the counterflow
    X_vapour_counter_eth_surface(ii+1) = X_vapour_counter_surface(ii+1)*X_vapour_eth_counter(ii+1);
    X_vapour_counter_HFA_surface(ii+1) = X_vapour_counter_surface(ii+1)*X_vapour_HFA_counter(ii+1);
    X_vapour_counter_air_surface(ii+1) = X_vapour_counter_surface(ii+1)*X_vapour_air_counter(ii+1);
        
    % Sum the molar fractions at the surface with the counterflow
    % composition at the surface
    X_vapour_eth_surface(ii+1) = X_vapour_eth_surface(ii+1) + X_vapour_counter_eth_surface(ii+1);
    X_vapour_HFA_surface(ii+1) = X_vapour_HFA_surface(ii+1) + X_vapour_counter_HFA_surface(ii+1);
    X_vapour_air_surface(ii+1) = X_vapour_counter_air_surface(ii+1);
    
    % Find molecular mass of vapour near the surface
	M_vapour_surface(ii+1) = M_air*X_vapour_air_surface(ii+1) + M_HFA*X_vapour_HFA_surface(ii+1) + M_eth*X_vapour_eth_surface(ii+1); 
    
    % Determine the mass fractions at the surface
    XM_tot_surface = X_vapour_HFA_surface(ii+1)*M_HFA + X_vapour_eth_surface(ii+1)*M_eth +  X_vapour_air_surface(ii+1)*M_air;
    Y_vapour_eth_surface(ii+1) = X_vapour_eth_surface(ii+1)*M_eth/(XM_tot_surface);
    Y_vapour_HFA_surface(ii+1) = X_vapour_HFA_surface(ii+1)*M_HFA/(XM_tot_surface);
    Y_vapour_air_surface(ii+1) = X_vapour_air_surface(ii+1)*M_air/(XM_tot_surface);
    
    % Determine mass fractions in the film using the 1/3 rule
    Y_vapour_eth_film(ii+1) = 2/3*Y_vapour_eth_surface(ii+1) + 1/3*Y_vapour_eth_counter(ii+1);     
    Y_vapour_HFA_film(ii+1) = 2/3*Y_vapour_HFA_surface(ii+1) + 1/3*Y_vapour_HFA_counter(ii+1); 
    Y_vapour_air_film(ii+1) = 1 - Y_vapour_HFA_film(ii+1) - Y_vapour_eth_film(ii+1);

    % Determine molar fractions in the film
    Y_M_tot_film = Y_vapour_HFA_film(ii+1)/M_HFA + Y_vapour_eth_film(ii+1)/M_eth + Y_vapour_air_film(ii+1)/M_air;
    X_vapour_eth_film(ii+1) = (Y_vapour_eth_film(ii+1)/M_eth)/Y_M_tot_film;
    X_vapour_HFA_film(ii+1) = (Y_vapour_HFA_film(ii+1)/M_HFA)/Y_M_tot_film;
    X_vapour_air_film(ii+1) = (Y_vapour_air_film(ii+1)/M_air)/Y_M_tot_film;
    
	% Calculate the evaporation fractions. This is just a ration of the
    % mass fractions of vapour at the surface due to teh partial pressures
    f_eth(ii+1) = max([min_mass_frac,min([(P_partial_eth(ii+1)*M_eth)/(P_partial_HFA(ii+1)*M_HFA + P_partial_eth(ii+1)*M_eth),1-min_mass_frac])]);
    f_HFA(ii+1) = max([min_mass_frac,min([(P_partial_HFA(ii+1)*M_HFA)/(P_partial_HFA(ii+1)*M_HFA + P_partial_eth(ii+1)*M_eth),1-min_mass_frac])]);
    
    %% THERMODYNAMIC PROPERTIES %%
    
    % Vapour density
    rho_vapour_eth_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_eth_temperature,imp_eth_density_vapour);
    if Y_vapour_eth_counter(ii+1) > 0
        rho_vapour_eth_counter(ii+1) = property_interpolating_function(T_counter(ii+1),imp_eth_temperature,imp_eth_density_vapour);
    end
    rho_vapour_HFA_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_HFA_temperature,imp_HFA_density_vapour);
    if Y_vapour_HFA_counter(ii+1) > 0
        rho_vapour_HFA_counter(ii+1) = property_interpolating_function(T_counter(ii+1),imp_HFA_temperature,imp_HFA_density_vapour);
    end
    rho_vapour_air_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_air_temperature,imp_air_density_vapour); 
    rho_vapour_air_counter(ii+1) = property_interpolating_function(T_counter(ii+1),imp_air_temperature,imp_air_density_vapour);

    % Liquid Density 
    rho_liquid_eth_mean(ii+1) = property_interpolating_function(T_droplet_mean(ii+1),imp_eth_temperature,imp_eth_density_liquid);
    rho_liquid_eth_surface(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_eth_temperature,imp_eth_density_liquid);
    rho_liquid_HFA_mean(ii+1) = property_interpolating_function(T_droplet_mean(ii+1),imp_HFA_temperature,imp_HFA_density_liquid);
    rho_liquid_HFA_surface(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_HFA_temperature,imp_HFA_density_liquid);
    
    % Vapour mixture density
    rho_vapour_film(ii+1) = (Y_vapour_air_film(ii+1)/rho_vapour_air_film(ii+1) + Y_vapour_eth_film(ii+1)/rho_vapour_eth_film(ii+1) + Y_vapour_HFA_film(ii+1)/rho_vapour_HFA_film(ii+1))^-1;

    % Liquid mixture density
    rho_liquid_mean(ii+1) = ((Y_liquid_eth_mean(ii+1)/rho_liquid_eth_mean(ii+1)) + (Y_liquid_HFA_mean(ii+1))/rho_liquid_HFA_mean(ii+1))^-1;
    rho_liquid_surface(ii+1) = ((Y_liquid_eth_surface(ii+1)/rho_liquid_eth_surface(ii+1)) + (Y_liquid_HFA_surface(ii+1))/rho_liquid_HFA_surface(ii+1))^-1;
    
	% Calculate mean density for iteration 1 for the liquid component if
    % ITC model is used
    if ii == 0 || strcmp(user_solver_method,'ITC') == 1
        rho_liquid_mean(ii+1) = ((Y_liquid_eth_mean(ii+1)/rho_liquid_eth_mean(ii+1)) + (Y_liquid_HFA_mean(ii+1))/rho_liquid_HFA_mean(ii+1))^-1;
    end
    
    % Vapour viscosity
    mu_vapour_eth_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_eth_temperature,imp_eth_mu_viscosity_vapour);
    mu_vapour_HFA_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_HFA_temperature,imp_HFA_mu_viscosity_vapour);
    mu_vapour_air_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_air_temperature,imp_air_mu_viscosity_vapour);
    
    % Liquid viscosity
    mu_liquid_eth_mean(ii+1) = property_interpolating_function(T_droplet_mean(ii+1),imp_eth_temperature,imp_eth_mu_viscosity_liquid);
    mu_liquid_eth_surface(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_eth_temperature,imp_eth_mu_viscosity_liquid);
    mu_liquid_HFA_mean(ii+1) = property_interpolating_function(T_droplet_mean(ii+1),imp_HFA_temperature,imp_HFA_mu_viscosity_liquid);
    mu_liquid_HFA_surface(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_HFA_temperature,imp_HFA_mu_viscosity_liquid);
    
    % Vapour kinematic viscosity
    nu_eth_vapour_film(ii+1) = mu_vapour_eth_film(ii+1)/rho_vapour_eth_film(ii+1);
    nu_HFA_vapour_film(ii+1) = mu_vapour_HFA_film(ii+1)/rho_vapour_HFA_film(ii+1);
    nu_air_vapour_film(ii+1) = mu_vapour_air_film(ii+1)/rho_vapour_air_film(ii+1);
    
    % Determine the viscosity of the film. First find the mixture for air
    % and HFA, then mix that with ethanol. We call the air-HFA the mix.
    % Get the mixture molar fractions
    X_mix_HFA = X_vapour_HFA_film(ii+1)/(X_vapour_HFA_film(ii+1)+X_vapour_air_film(ii+1));
    X_mix_air = X_vapour_air_film(ii+1)/(X_vapour_HFA_film(ii+1)+X_vapour_air_film(ii+1));
    
    % Get the mixture molar mass
    M_mix = X_mix_HFA*M_HFA + X_mix_air*M_air;
    
    % Get the mixture molar fraction
    X_mix = X_vapour_HFA_film(ii+1) + X_vapour_air_film(ii+1);
    
    % Get the  vapour ratios
    phi_HFA_air(ii+1) = sqrt(M_air/M_HFA);
    phi_air_HFA(ii+1) = sqrt(M_HFA/M_air);
    phi_mix_eth(ii+1) = sqrt(M_eth/M_mix);
    phi_eth_mix(ii+1) = sqrt(M_mix/M_eth);
    
    % Get the liquid ratios
    phi_HFA_eth(ii+1) = sqrt(M_eth/M_HFA);
    phi_eth_HFA(ii+1) = sqrt(M_HFA/M_eth);
    
    % Find viscosity in the mix of HFA and air
    mu_vapour_mix(ii+1) = (X_mix_HFA*mu_vapour_HFA_film(ii+1))/(X_mix_HFA+X_mix_air*phi_HFA_air(ii+1)) + (X_mix_air*mu_vapour_air_film(ii+1))/(X_mix_HFA*phi_air_HFA(ii+1)+X_mix_air);
    
    % Vapour mixture viscosity
    mu_vapour_film(ii+1) = (X_mix*mu_vapour_mix(ii+1))/(X_mix+X_vapour_eth_film(ii+1)*phi_mix_eth(ii+1)) + (X_vapour_eth_film(ii+1)*mu_vapour_eth_film(ii+1))/(X_mix*phi_eth_mix(ii+1)+X_vapour_eth_film(ii+1));
        
    % Vapour mixture kinematic viscosity 
    nu_vapour_film(ii+1) = mu_vapour_film(ii+1)/rho_vapour_film(ii+1);
    
    % Liquid mixture viscosity. Empirical relationn (Johnson 2018 RDD)
    % liquid viscosity of mixture of HFA134a and ethanol. Assume similar
    % relation for other propellants
    mu_liquid_mean(ii+1) = (X_liquid_HFA_mean(ii+1)*mu_liquid_HFA_mean(ii+1)^(-0.45) + X_liquid_eth_mean(ii+1)*mu_liquid_eth_mean(ii+1)^(-0.45))^(-1/0.45);
    mu_liquid_surface(ii+1) = (X_liquid_HFA_surface(ii+1)*mu_liquid_HFA_surface(ii+1)^(-0.45) + X_liquid_eth_surface(ii+1)*mu_liquid_eth_surface(ii+1)^(-0.45))^(-1/0.45);
    
    % Liquid Kinematic viscosity
    nu_liquid_mean(ii+1) = mu_liquid_mean(ii+1)/rho_liquid_mean(ii+1);
    nu_liquid_surface(ii+1) = mu_liquid_surface(ii+1)/(ii+1);
    
    % Vapour thermal conductivity
    lambda_vapour_eth_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_eth_temperature,imp_eth_conductivity_vapour);
    lambda_vapour_HFA_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_HFA_temperature,imp_HFA_conductivity_vapour);
    lambda_vapour_air_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_air_temperature,imp_air_conductivity_vapour);
    
    % Liquid thermal conductivity
    lambda_liquid_eth_mean(ii+1) = property_interpolating_function(T_droplet_mean(ii+1),imp_eth_temperature,imp_eth_conductivity_liquid);
    lambda_liquid_eth_surface(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_eth_temperature,imp_eth_conductivity_liquid); 
    lambda_liquid_HFA_mean(ii+1) = property_interpolating_function(T_droplet_mean(ii+1),imp_HFA_temperature,imp_HFA_conductivity_liquid);
    lambda_liquid_HFA_surface(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_HFA_temperature,imp_HFA_conductivity_liquid);  
    
    % Now do similar to viscosity to get the vapour thermal conductivity
    lambda_vapour_mix(ii+1) = (X_mix_HFA*lambda_vapour_HFA_film(ii+1))/(X_mix_HFA+X_mix_air*phi_HFA_air(ii+1)) + (X_mix_air*lambda_vapour_air_film(ii+1))/(X_mix_HFA*phi_air_HFA(ii+1)+X_mix_air);
    
    % Vapour mixture thermal conductivity
    lambda_vapour_film(ii+1) = (X_mix*lambda_vapour_mix(ii+1))/(X_mix+X_vapour_eth_film(ii+1)*phi_mix_eth(ii+1)) + (X_vapour_eth_film(ii+1)*lambda_vapour_eth_film(ii+1))/(X_mix*phi_eth_mix(ii+1)+X_vapour_eth_film(ii+1));
    
    % Thermal conductivity of liquid
    lambda_liquid_mean(ii+1) = (X_liquid_HFA_mean(ii+1)*lambda_liquid_HFA_mean(ii+1))/(X_liquid_HFA_mean(ii+1)+X_liquid_eth_mean(ii+1)*phi_HFA_eth(ii+1)) + (X_liquid_eth_mean(ii+1)*lambda_liquid_eth_mean(ii+1))/(X_liquid_HFA_mean(ii+1)*phi_eth_HFA(ii+1)+X_liquid_eth_mean(ii+1));
    lambda_liquid_surface(ii+1) = (X_liquid_HFA_surface(ii+1)*lambda_liquid_HFA_surface(ii+1))/(X_liquid_HFA_surface(ii+1)+X_liquid_eth_surface(ii+1)*phi_HFA_eth(ii+1)) + (X_liquid_eth_surface(ii+1)*lambda_liquid_eth_surface(ii+1))/(X_liquid_HFA_surface(ii+1)*phi_eth_HFA(ii+1)+X_liquid_eth_surface(ii+1));
   
    % Vapour specific heat 
    c_p_vapour_eth_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_eth_temperature,imp_eth_Cp_vapour);
    c_p_vapour_HFA_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_HFA_temperature,imp_HFA_Cp_vapour);
    c_p_vapour_air_film(ii+1) = property_interpolating_function(T_BL_v(ii+1),imp_air_temperature,imp_air_Cp_vapour);
    
    % Liquid specific heat
    c_p_liquid_eth_mean(ii+1) = property_interpolating_function(T_droplet_mean(ii+1),imp_eth_temperature,imp_eth_Cp_liquid);
    c_p_liquid_eth_surface(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_eth_temperature,imp_eth_Cp_liquid);
    c_p_liquid_HFA_mean(ii+1) = property_interpolating_function(T_droplet_mean(ii+1),imp_HFA_temperature,imp_HFA_Cp_liquid);
    c_p_liquid_HFA_surface(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_HFA_temperature,imp_HFA_Cp_liquid);
    
    % Vapour mixture specific heat
    c_p_vapour_film(ii+1) = Y_vapour_HFA_film(ii+1)*c_p_vapour_HFA_film(ii+1) + Y_vapour_eth_film(ii+1)*c_p_vapour_eth_film(ii+1) + Y_vapour_air_film(ii+1)*c_p_vapour_air_film(ii+1);
    
    % Liquid mixture specific heat 
    c_p_liquid_mean(ii+1) = Y_liquid_eth_mean(ii+1)*c_p_liquid_eth_mean(ii+1) + Y_liquid_HFA_mean(ii+1)*c_p_liquid_HFA_mean(ii+1);
    c_p_liquid_surface(ii+1) = Y_liquid_eth_surface(ii+1)*c_p_liquid_eth_surface(ii+1) + Y_liquid_HFA_surface(ii+1)*c_p_liquid_HFA_surface(ii+1);
    
    % Liquid surface tension
    surf_ten_eth(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_eth_temperature,imp_eth_surf_ten);
    surf_ten_HFA(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_HFA_temperature,imp_HFA_surf_ten);
    surf_ten_mix(ii+1) = X_liquid_HFA_surface(ii+1)*surf_ten_HFA(ii+1) + X_liquid_eth_surface(ii+1)*surf_ten_eth(ii+1);
    
    % Determine the heat of vapourisation for each component
    L_eth(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_eth_temperature,imp_eth_latent_heat);
    L_HFA(ii+1) = property_interpolating_function(T_droplet_surface(ii+1),imp_HFA_temperature,imp_HFA_latent_heat);
    
    % Get lennard-jones parameters for the vapour surface and counterflow
    eps_vapour_counter(ii+1) = eps_HFA*X_vapour_HFA_counter(ii+1) + eps_eth*X_vapour_eth_counter(ii+1) + eps_air*X_vapour_air_counter(ii+1);
    sigma_vapour_counter(ii+1) = sigma_HFA*X_vapour_HFA_counter(ii+1) + sigma_eth*X_vapour_eth_counter(ii+1) + sigma_air*X_vapour_air_counter(ii+1);
    eps_vapour_surface(ii+1) = eps_HFA*X_vapour_HFA_surface(ii+1) + eps_eth*X_vapour_eth_surface(ii+1) + eps_air*X_vapour_air_surface(ii+1);
    sigma_vapour_surface(ii+1) = sigma_HFA*X_vapour_HFA_surface(ii+1) + sigma_eth*X_vapour_eth_surface(ii+1) + sigma_air*X_vapour_air_surface(ii+1);
    
    % Get lennard-jones vapour mixture
    sigma_vapour_film(ii+1) = 0.5*(sigma_vapour_counter(ii+1) + sigma_vapour_surface(ii+1));
    
    % Get mixture molecular mass
    M_vs(ii+1) = 2*(1/M_vapour_counter(ii+1) + 1/M_vapour_surface(ii+1))^-1;
    
    % Find dimensionless collision integral
    omega_D(ii+1) = omega_integral(T_BL_v(ii+1),k_B,eps_vapour_counter(ii+1),eps_vapour_surface(ii+1));
    
    % Determine diffusion coefficient for the vapour phase
    D_vapour_film(ii+1) = ((3.03 - (0.98/(M_vs(ii+1))^0.5))*10^(-7)*T_BL_v(ii+1)^(3/2))/(p_amb_bar*(M_vs(ii+1))^0.5*(sigma_vapour_film(ii+1))^2*omega_D(ii+1));    
    
    %% DIMENSIONLESS NUMBERS %%
      
    % Calculate the spalding mass transfer number
    if (1 - Y_vapour_HFA_surface(ii+1) - Y_vapour_eth_surface(ii+1)) == 0
        B_M(ii+1) = 0;    
    else
        B_M(ii+1) = ((Y_vapour_HFA_surface(ii+1) + Y_vapour_eth_surface(ii+1)) - (Y_vapour_HFA_counter(ii+1) + Y_vapour_eth_counter(ii+1)))/(1 - Y_vapour_HFA_surface(ii+1) - Y_vapour_eth_surface(ii+1)); 
    end

    % Error check (Outwards mass transfer => B_M>0)
    if B_M(ii+1) <= 0
        warning('B_M < 0. Surface mass fraction > ambient mass fraction\n')
        B_M(ii+1) = 0;
    end

    % Determine relevant dimensionless numbers
    Sc(ii+1) = mu_vapour_film(ii+1)/(rho_vapour_film(ii+1)*D_vapour_film(ii+1));
    Pr(ii+1) = mu_vapour_film(ii+1)*c_p_vapour_film(ii+1)/lambda_vapour_film(ii+1); 
    Le(ii+1) = Sc(ii+1)/Pr(ii+1);
    
    % Calculate Re off droplet diameter (function doubles length scale)
    Re_d(ii+1) = Reynolds_number(r_droplet(ii+1), U_droplet(ii+1), U_counter(ii+1), rho_vapour_film(ii+1), mu_vapour_film(ii+1));
    
    % SCalculate droplet Weber number
    We_d(ii+1) = rho_liquid_mean(ii+1)*(U_droplet(ii+1)-U_counter(ii+1))^2*(2*r_droplet(ii+1))/surf_ten_mix(ii+1);
    
    % If B_M = 0. No outward evaporation. 
    if B_M(ii+1) == 0
        F_M(ii+1) = 1;
        Sh(ii+1) = 2*(1+((1+Re_d(ii+1)*Sc(ii+1))^(1/3)*max([1,Re_d(ii+1)^0.077])-1)/(2*F_M(ii+1)));
        F_T(ii+1) = 1;
        B_T(ii+1) = 0;
        Nu_d(ii+1) = 2*(1+((1+Re_d(ii+1)*Pr(ii+1))^(1/3)*max([1,Re_d(ii+1)^0.077])-1)/(2*F_T(ii+1)));
    
    % If non-zero B_M    
    else
        % Calculate the boundary layer correcting factor for mass F_M
        F_M(ii+1) = (1+B_M(ii+1))^(0.7)*(log(1+B_M(ii+1))/B_M(ii+1));
        
        % Calculate Sherwood number
        Sh(ii+1) = 2*(log(1+B_M(ii+1))/B_M(ii+1))*(1+((1+Re_d(ii+1)*Sc(ii+1))^(1/3)*max([1,Re_d(ii+1)^0.077])-1)/(2*F_M(ii+1)));
        
        % Nusselt number needs to be solved iteratively.
        % Reset counters.
        count = 0;
        
        % Kickstart iterative solver for Nu
        if ii == 0           
            B_T_old = 0.1;
        elseif ii > 0
            if B_M(ii) == 0
                B_T_old = 0.1;
            else
                B_T_old = B_T(ii);
            end
        end
        
        % Initiate loop by forcing tolerance above threshold
        tol = B_T_tol + 1;
        
        % Loop to find Nu
        while tol > B_T_tol
            
            % Exit with error if diverging or slow convergence
            count = count + 1;
            if count > 1E6
                error('stuck in loop ii = %d',ii)
            end
            
            % Find the Nu number, F_T coefficient and B_T
            F_T_new = (1+B_T_old)^(0.7)*(log(1+B_T_old)/B_T_old);
            Nu_new = 2*(log(1+B_T_old)/B_T_old)*(1+((1+Re_d(ii+1)*Pr(ii+1))^(1/3)*max([1,Re_d(ii+1)^0.077])-1)/(2*F_T_new));
            B_T_new = (1 + B_M(ii+1))^((c_p_liquid_mean(ii+1)*Sh(ii+1)*B_M(ii+1)*log(1+B_T_old))/(Le(ii+1)*c_p_vapour_film(ii+1)*Nu_new*B_T_old*log(1+B_M(ii+1))))-1; % Find new B_T
            
            % Check tolerance
            tol = abs(B_T_old - B_T_new);
            B_T_old = B_T_new;
        end
        
        % Update heat transfer dimensionless numbers
        B_T(ii+1) = B_T_new;
        F_T(ii+1) = F_T_new;
        Nu_d(ii+1) = Nu_new;
        
    end
    
    %% FIRST STEP INITIALISATION %%
    if ii == 0
        if strcmp(user_solver_method,'ETC')==1
                    
            % We are going to double the number of points so half them
            % intially
            tot_no_gridpoints = ceil(tot_no_gridpoints/2);
            
            % Define the spacing between each grid point using x^2
            index_spacing = zeros(1,tot_no_gridpoints);
            point = zeros(1,tot_no_gridpoints);
            for xx = 0:tot_no_gridpoints-1
                index_spacing(xx+1) = xx^2;
            end
            
            % Get the total length of the radial matrix
            vec = 1:sum(index_spacing)+1;
            
            % Select the points such that we have spacing increasing by x^2
            % gridpoint with each step
            for xx = 1:length(index_spacing)
                point(xx) = vec(end-sum(index_spacing(1:xx)));
            end
            
            % Make refinement at surface and select the index pts for the
            % discretisation
            point = [vec(end)-point(1:end-1)+1, flip(point)+length(vec)];
            
            % Get the overall length of the vector
            vec = 1:2*length(vec);
            
            % define the radial matrix
            r_location = linspace(0,r_droplet(1),length(vec));
            r_location = r_location(point);
            
            % Get the final number of points
            true_tot_no_gridpoints = length(r_location);
            
            % Postulate a value for the liquid diffusivity
            D_liquid_mean(1) = user_D_liquid_sweep(sensitivity_sweep_D)*D_liquid;
            lambda_liquid_mean(1) = user_K_liquid_sweep(sensitivity_sweep_K)*lambda_liquid_mean(1);
            
        	% Set initial variables for DE's
            R_0 = r_droplet(1);
            alpha_liquid_mean(1) = lambda_liquid_mean(1)/(c_p_liquid_mean(1)*rho_liquid_mean(1));
            alpha_liquid_0 = alpha_liquid_mean(1);
            T_liquid_0 = T_droplet_mean(1);
            r_s(1) = r_droplet(1)/R_0;
                        
            % Set initial matrices for DE's
            DE_Z_RHS(:,1) = (T_droplet_surface(1)-T_liquid_0)/T_liquid_0*ones(true_tot_no_gridpoints,1);%
            DE_Y_eth_RHS(:,1) = Y_liquid_eth_surface(1)*ones(true_tot_no_gridpoints,1);%
            DE_Y_HFA_RHS(:,1) = Y_liquid_HFA_surface(1)*ones(true_tot_no_gridpoints,1);%
            
            DE_Z_LHS(:,1) = (T_droplet_surface(1)-T_liquid_0)/T_liquid_0*ones(true_tot_no_gridpoints,1);%
            DE_Y_eth_LHS(:,1) = Y_liquid_eth_surface(1)*ones(true_tot_no_gridpoints,1);%
            DE_Y_HFA_LHS(:,1) = Y_liquid_HFA_surface(1)*ones(true_tot_no_gridpoints,1);%
            
        end
    end
    
    % Step forward in time
    ii = ii + 1;
    
    % Print interation number every 1000 iterations
    if mod(ii,1000) == 0
        fprintf('iteration #%f\n',ii)
    end
end

% Save data from simulation
filename_var = join([save_folder,user_propellant_name(sensitivity_sweep_HFA),'\D=',user_D_liquid_sweep(sensitivity_sweep_D),'\Yeth=',num2str(current_iter_initial_Y_eth_l)],'');
mkdir(filename_var);
save(join([filename_var,'\variables.mat'],''))


end
end
end
end
