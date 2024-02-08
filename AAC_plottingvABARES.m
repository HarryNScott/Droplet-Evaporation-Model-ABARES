close all
clc
clearvars -except DATA_MAT

%% FIGURE 1 droplet radius vs time %%

% Parameters: Ref Type & P model
% D = 1e-8;
% Yeth = 0.15
% REF = 134, 152 & 1234 (colour)
% P = normal and 1 (line type)

% Single case for ABARES speed
HFA_IMPORT = {'HFA152a'};
D_IMPORT = {'D=1000'};
Y_ETH_IMPORT = {'Yeth=0.15'};

% Reset counterts
count = 0;
leg_count = 0;

% Format figure
if count == 0
    figure
    hold on
    xlabel('time [s]')
    ylabel('R [m]')
    count = count+1;
end

% Sweep through the required refrigerants
for xx = 1:length(HFA_IMPORT)
    
    % Get the location of the refrigerants in the data cell
    for yy = 1:size(DATA_MAT,2)
        if strcmp(cell2mat(HFA_IMPORT(xx)),cell2mat(DATA_MAT{1,yy}))
            HFA_IMPORT_loc(xx) = yy;
        end
    end
    
    % get all the data for that refrigerant
    DATA_TOT = DATA_MAT{2,HFA_IMPORT_loc(xx)};
    
    % Get the location of the diffusivities in the data cell
    for yy = 1:length(D_IMPORT)
        for zz = 1:size(DATA_TOT,1)
            if strcmp(cell2mat(D_IMPORT(yy)),DATA_TOT{zz,1})
                D_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the location of the mass fractions in the data cell
    for yy = 1:length(Y_ETH_IMPORT)
        for zz = 1:size(DATA_TOT,2)
            if strcmp(cell2mat(Y_ETH_IMPORT(yy)),DATA_TOT{1,zz})
                Y_ETH_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the DATA for the current test
    for yy = 1:length(D_IMPORT_loc)
        for zz = 1:length(Y_ETH_IMPORT_loc)
            DATA_CURRENT = DATA_TOT{D_IMPORT_loc(yy),Y_ETH_IMPORT_loc(zz)};
            
            % Pull the radius and time data and plot it
            rad_temp = DATA_CURRENT{2,1};
            time_temp = DATA_CURRENT{2,13};
            
            % Crop data to 15% of the radius so we are consistent
            [min_rad,loc] = min(abs(rad_temp - rad_temp(1)*0.15));
            if loc == length(rad_temp)
                error('radius > 15%')
            end
            
            % plot the results and save teh legend
            plot(time_temp(1:loc),rad_temp(1:loc))
            leg_count = leg_count+1;
            legend_names(leg_count) = join([HFA_IMPORT(xx),D_IMPORT(yy),Y_ETH_IMPORT(zz)],' ');
            
        end
    end
end

legend(legend_names)
clearvars -except DATA_MAT

%P** TO COME

%% FIGURE 2 droplet velocity vs time %%

% Parameters: Ref Type & Yeth
% D = 1e-8;
% Yeth = 0.025, 0.15 (line type)
% REF = 134, 152 & 1234 (colour)
% P  = normal

% Single case for ABARES speed
HFA_IMPORT = {'HFA152a'};
D_IMPORT = {'D=1000'};
Y_ETH_IMPORT = {'Yeth=0.15'};

% Reset counterts
count = 0;
leg_count = 0;

% Format figure
if count == 0
    figure
    hold on
    xlabel('time [s]')
    ylabel('U [m/s]')
    count = count+1;
end

% Sweep through the required refrigerants
for xx = 1:length(HFA_IMPORT)
    
    % Get the location of the refrigerants in the data cell
    for yy = 1:size(DATA_MAT,2)
        if strcmp(cell2mat(HFA_IMPORT(xx)),cell2mat(DATA_MAT{1,yy}))
            HFA_IMPORT_loc(xx) = yy;
        end
    end
    
    % get all the data for that refrigerant
    DATA_TOT = DATA_MAT{2,HFA_IMPORT_loc(xx)};
    
    % Get the location of the diffusivities in the data cell
    for yy = 1:length(D_IMPORT)
        for zz = 1:size(DATA_TOT,1)
            if strcmp(cell2mat(D_IMPORT(yy)),DATA_TOT{zz,1})
                D_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the location of the mass fractions in the data cell
    for yy = 1:length(Y_ETH_IMPORT)
        for zz = 1:size(DATA_TOT,2)
            if strcmp(cell2mat(Y_ETH_IMPORT(yy)),DATA_TOT{1,zz})
                Y_ETH_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the DATA for the current test
    for yy = 1:length(D_IMPORT_loc)
        for zz = 1:length(Y_ETH_IMPORT_loc)
            DATA_CURRENT = DATA_TOT{D_IMPORT_loc(yy),Y_ETH_IMPORT_loc(zz)};
            
            % Pull the radius and time data and plot it
            rad_temp = DATA_CURRENT{2,1};
            time_temp = DATA_CURRENT{2,13};
            U_temp = DATA_CURRENT{2,2};
            
            % Crop data to 15% of the radius so we are consistent
            [min_rad,loc] = min(abs(rad_temp - rad_temp(1)*0.15));
            if loc == length(rad_temp)
                error('radius > 15%')
            end
            
            % plot the results and save teh legend
            plot(time_temp(1:loc),U_temp(1:loc))
            leg_count = leg_count+1;
            legend_names(leg_count) = join([HFA_IMPORT(xx),D_IMPORT(yy),Y_ETH_IMPORT(zz)],' ');
            
        end
    end
end

legend(legend_names)
clearvars -except DATA_MAT

%% FIGURE 3 droplet radius vs time (x 4 different Yeth) %%

% Parameters: Yeth & D
% D = 1e-1, 1e0, 1e1, 1e2, 5e-1, 5e0 (colour)
% Yeth = 0.025, 0.05, 0.15, 0.5 (plot)
% REF = 152

% Single case for ABARES speed
HFA_IMPORT = {'HFA152a'};
D_IMPORT = {'D=1000'};
Y_ETH_IMPORT = {'Yeth=0.15'};

% Reset counterts
count = 0;
leg_count = 0;

% Format figure
if count == 0
    figure
    subplot(2,2,1)
    xlabel('time [s]')
    ylabel('R [m]')
    hold on
    subplot(2,2,1)
    xlabel('time [s]')
    ylabel('R [m]')
    hold on
    subplot(2,2,1)
    xlabel('time [s]')
    ylabel('R [m]')
    hold on
    subplot(2,2,1)
    xlabel('time [s]')
    ylabel('R [m]')
    hold on
    count = count+1;
end

% Sweep through the required refrigerants
for xx = 1:length(HFA_IMPORT)
    
    % Get the location of the refrigerants in the data cell
    for yy = 1:size(DATA_MAT,2)
        if strcmp(cell2mat(HFA_IMPORT(xx)),cell2mat(DATA_MAT{1,yy}))
            HFA_IMPORT_loc(xx) = yy;
        end
    end
    
    % get all the data for that refrigerant
    DATA_TOT = DATA_MAT{2,HFA_IMPORT_loc(xx)};
    
    % Get the location of the diffusivities in the data cell
    for yy = 1:length(D_IMPORT)
        for zz = 1:size(DATA_TOT,1)
            if strcmp(cell2mat(D_IMPORT(yy)),DATA_TOT{zz,1})
                D_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the location of the mass fractions in the data cell
    for yy = 1:length(Y_ETH_IMPORT)
        for zz = 1:size(DATA_TOT,2)
            if strcmp(cell2mat(Y_ETH_IMPORT(yy)),DATA_TOT{1,zz})
                Y_ETH_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the DATA for the current test
    for yy = 1:length(D_IMPORT_loc)
        for zz = 1:length(Y_ETH_IMPORT_loc)
            DATA_CURRENT = DATA_TOT{D_IMPORT_loc(yy),Y_ETH_IMPORT_loc(zz)};
            
            % Pull the radius and time data and plot it
            rad_temp = DATA_CURRENT{2,1};
            time_temp = DATA_CURRENT{2,13};
            
            % Crop data to 15% of the radius so we are consistent
            [min_rad,loc] = min(abs(rad_temp - rad_temp(1)*0.15));
            if loc == length(rad_temp)
                error('radius > 15%')
            end
            
            % plot the results and save teh legend
            subplot(2,2,zz)
            plot(time_temp(1:loc),rad_temp(1:loc))
            legend_names(yy,zz) = join([HFA_IMPORT(xx),D_IMPORT(yy),Y_ETH_IMPORT(zz)],' ');
            legend(legend_names(:,zz))
        end
    end
end

clearvars -except DATA_MAT

%% FIGURE 4 droplet f_eth vs time %%
% Parameters: Yeth & D
% D = 1e-1, 1e0, 1e1 (colour)
% Y_eth = 0.025, 0.05, 0.15 (line type)
% REF = 152

% Single case for ABARES speed
HFA_IMPORT = {'HFA152a'};
D_IMPORT = {'D=1000'};
Y_ETH_IMPORT = {'Yeth=0.15'};

% Reset counterts
count = 0;
leg_count = 0;

% Format figure
if count == 0
    figure
    hold on
    xlabel('time [s]')
    ylabel('f_{eth}')
    count = count+1;
end

% Sweep through the required refrigerants
for xx = 1:length(HFA_IMPORT)
    
    % Get the location of the refrigerants in the data cell
    for yy = 1:size(DATA_MAT,2)
        if strcmp(cell2mat(HFA_IMPORT(xx)),cell2mat(DATA_MAT{1,yy}))
            HFA_IMPORT_loc(xx) = yy;
        end
    end
    
    % get all the data for that refrigerant
    DATA_TOT = DATA_MAT{2,HFA_IMPORT_loc(xx)};
    
    % Get the location of the diffusivities in the data cell
    for yy = 1:length(D_IMPORT)
        for zz = 1:size(DATA_TOT,1)
            if strcmp(cell2mat(D_IMPORT(yy)),DATA_TOT{zz,1})
                D_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the location of the mass fractions in the data cell
    for yy = 1:length(Y_ETH_IMPORT)
        for zz = 1:size(DATA_TOT,2)
            if strcmp(cell2mat(Y_ETH_IMPORT(yy)),DATA_TOT{1,zz})
                Y_ETH_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the DATA for the current test
    for yy = 1:length(D_IMPORT_loc)
        for zz = 1:length(Y_ETH_IMPORT_loc)
            DATA_CURRENT = DATA_TOT{D_IMPORT_loc(yy),Y_ETH_IMPORT_loc(zz)};
            
            % Pull the radius and time data and plot it
            rad_temp = DATA_CURRENT{2,1};
            time_temp = DATA_CURRENT{2,13};
            f_eth_temp = DATA_CURRENT{2,11};
            
            % Crop data to 15% of the radius so we are consistent
            [min_rad,loc] = min(abs(rad_temp - rad_temp(1)*0.15));
            if loc == length(rad_temp)
                error('radius > 15%')
            end
            
            % plot the results and save teh legend
            plot(time_temp(1:loc),f_eth_temp(1:loc))
            leg_count = leg_count+1;
            legend_names(leg_count) = join([HFA_IMPORT(xx),D_IMPORT(yy),Y_ETH_IMPORT(zz)],' ');
            
        end
    end
end

legend(legend_names)
clearvars -except DATA_MAT

%% FIGURE 5 droplet Y_eth_surface vs time %%
% Parameters: Yeth & D
% D = 1e-1, 1e0, 1e1 (colour)
% Y_eth = 0.025, 0.05, 0.15 (line type)
% REF = 152

% Single case for ABARES speed
HFA_IMPORT = {'HFA152a'};
D_IMPORT = {'D=1000'};
Y_ETH_IMPORT = {'Yeth=0.15'};

% Reset counterts
count = 0;
leg_count = 0;

% Format figure
if count == 0
    figure
    hold on
    xlabel('time [s]')
    ylabel('Y_{eth,surface}')
    count = count+1;
end

% Sweep through the required refrigerants
for xx = 1:length(HFA_IMPORT)
    
    % Get the location of the refrigerants in the data cell
    for yy = 1:size(DATA_MAT,2)
        if strcmp(cell2mat(HFA_IMPORT(xx)),cell2mat(DATA_MAT{1,yy}))
            HFA_IMPORT_loc(xx) = yy;
        end
    end
    
    % get all the data for that refrigerant
    DATA_TOT = DATA_MAT{2,HFA_IMPORT_loc(xx)};
    
    % Get the location of the diffusivities in the data cell
    for yy = 1:length(D_IMPORT)
        for zz = 1:size(DATA_TOT,1)
            if strcmp(cell2mat(D_IMPORT(yy)),DATA_TOT{zz,1})
                D_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the location of the mass fractions in the data cell
    for yy = 1:length(Y_ETH_IMPORT)
        for zz = 1:size(DATA_TOT,2)
            if strcmp(cell2mat(Y_ETH_IMPORT(yy)),DATA_TOT{1,zz})
                Y_ETH_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the DATA for the current test
    for yy = 1:length(D_IMPORT_loc)
        for zz = 1:length(Y_ETH_IMPORT_loc)
            DATA_CURRENT = DATA_TOT{D_IMPORT_loc(yy),Y_ETH_IMPORT_loc(zz)};
            
            % Pull the radius and time data and plot it
            rad_temp = DATA_CURRENT{2,1};
            time_temp = DATA_CURRENT{2,13};
            Y_eth_surf_temp = DATA_CURRENT{2,7};
            
            % Crop data to 15% of the radius so we are consistent
            [min_rad,loc] = min(abs(rad_temp - rad_temp(1)*0.15));
            if loc == length(rad_temp)
                error('radius > 15%')
            end
            
            % plot the results and save teh legend
            plot(time_temp(1:loc),Y_eth_surf_temp(1:loc))
            leg_count = leg_count+1;
            legend_names(leg_count) = join([HFA_IMPORT(xx),D_IMPORT(yy),Y_ETH_IMPORT(zz)],' ');
            
        end
    end
end

legend(legend_names)
clearvars -except DATA_MAT

%% FIGURE 6 droplet temperature vs time %%
% Parameters: Yeth & D
% D = 1e-1, 1e0, 1e1 (colour)
% Y_eth = 0.025, 0.05, 0.15 (line type)
% REF = 152

% Single case for ABARES speed
HFA_IMPORT = {'HFA152a'};
D_IMPORT = {'D=1000'};
Y_ETH_IMPORT = {'Yeth=0.15'};

% Reset counterts
count = 0;
leg_count = 0;

% Format figure
if count == 0
    figure
    hold on
    xlabel('time [s]')
    ylabel('T_{surface}')
    count = count+1;
end

% Sweep through the required refrigerants
for xx = 1:length(HFA_IMPORT)
    
    % Get the location of the refrigerants in the data cell
    for yy = 1:size(DATA_MAT,2)
        if strcmp(cell2mat(HFA_IMPORT(xx)),cell2mat(DATA_MAT{1,yy}))
            HFA_IMPORT_loc(xx) = yy;
        end
    end
    
    % get all the data for that refrigerant
    DATA_TOT = DATA_MAT{2,HFA_IMPORT_loc(xx)};
    
    % Get the location of the diffusivities in the data cell
    for yy = 1:length(D_IMPORT)
        for zz = 1:size(DATA_TOT,1)
            if strcmp(cell2mat(D_IMPORT(yy)),DATA_TOT{zz,1})
                D_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the location of the mass fractions in the data cell
    for yy = 1:length(Y_ETH_IMPORT)
        for zz = 1:size(DATA_TOT,2)
            if strcmp(cell2mat(Y_ETH_IMPORT(yy)),DATA_TOT{1,zz})
                Y_ETH_IMPORT_loc(yy) = zz;
            end
        end
    end
    
    % Get the DATA for the current test
    for yy = 1:length(D_IMPORT_loc)
        for zz = 1:length(Y_ETH_IMPORT_loc)
            DATA_CURRENT = DATA_TOT{D_IMPORT_loc(yy),Y_ETH_IMPORT_loc(zz)};
            
            % Pull the radius and time data and plot it
            rad_temp = DATA_CURRENT{2,1};
            time_temp = DATA_CURRENT{2,13};
            T_temp = DATA_CURRENT{2,3};
            
            % Crop data to 15% of the radius so we are consistent
            [min_rad,loc] = min(abs(rad_temp - rad_temp(1)*0.15));
            if loc == length(rad_temp)
                error('radius > 15%')
            end
            
            % plot the results and save teh legend
            plot(time_temp(1:loc),T_temp(1:loc))
            leg_count = leg_count+1;
            legend_names(leg_count) = join([HFA_IMPORT(xx),D_IMPORT(yy),Y_ETH_IMPORT(zz)],' ');
            
        end
    end
end

legend(legend_names)
clearvars -except DATA_MAT

%% FIGURE 7 droplet evap time vs D %%
% Parameters: Ref type & Yeth
% D = 1e-1, 1e0, 1e1, 1e2, 5e-1, 5e0 (x axis)
% Y_eth = 0.025, 0.05, 0.15 (line type)
% REF = 134, 152 & 1234 (colour)

% Single case for ABARES speed
HFA_IMPORT = {'HFA134a','HFA152a','HFA1234zeE'};
D_IMPORT = {'D=100','D=10','D=5','D=1'}; %,'D=0.5'}; %,'D=1e-1'};
Y_ETH_IMPORT = {'Yeth=0.025','Yeth=0.05','Yeth=0.15'};



