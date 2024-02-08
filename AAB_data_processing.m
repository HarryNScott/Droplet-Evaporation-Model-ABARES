close all
clear all
clc

% Location of results
overarching_folder_loc = 'H:\Data for AFMC\';

% Get the directory of Propellants
propellant_folders = dir(overarching_folder_loc);
propellant_folders = struct2cell(propellant_folders([propellant_folders.isdir]));
propellant_folders = propellant_folders(1,:);

% create logical matrix
log_mat_HFA = 0;

% Get the locations of the relevant propellant folders
for xx_a = 1:length(propellant_folders)
    
    temp = cell2mat(propellant_folders(xx_a));
    if length(temp) >= 3
        if strcmp(temp(1:3),'HFA') == 1
            log_mat_HFA(xx_a) = 1;
        else
            log_mat_HFA(xx_a) = 0;
        end
    else
        log_mat_HFA(xx_a) = 0;
    end
end
log_mat_HFA = logical(log_mat_HFA);
propellant_folders = propellant_folders(log_mat_HFA);

% Sweep through propellant folders
for ii_a = 1:length(propellant_folders)
    
    % Get the directory of Diffusivities
    diffusivity_folders = dir(join([overarching_folder_loc,cell2mat(propellant_folders(ii_a)),'\'],''));
    diffusivity_folders = struct2cell(diffusivity_folders([diffusivity_folders.isdir]));
    diffusivity_folders = diffusivity_folders(1,:);
    
    % create logical matrix
    log_mat_D = 0;
    
    % Get the locations of the relevant diffusivity folders
    for xx_a = 1:length(diffusivity_folders)
        
        temp = cell2mat(diffusivity_folders(xx_a));
        if length(temp) >= 2
            if strcmp(temp(1:2),'D=') == 1
                log_mat_D(xx_a) = 1;
            else
                log_mat_D(xx_a) = 0;
            end
        else
            log_mat_D(xx_a) = 0;
        end
    end
    log_mat_D = logical(log_mat_D);
    diffusivity_folders = diffusivity_folders(log_mat_D);
    
    % Sweep through all diffusivity folders
    for jj_a = 1:length(diffusivity_folders)
        
        % Get the directory of Mass Fractions
        mass_frac_folders = dir(join([overarching_folder_loc,cell2mat(propellant_folders(ii_a)),'\',cell2mat(diffusivity_folders(jj_a))],''));
        mass_frac_folders = struct2cell(mass_frac_folders([mass_frac_folders.isdir]));
        mass_frac_folders = mass_frac_folders(1,:);
        
        % create logical matrix
        log_mat_Y = 0;
        
        % Get the locations of the relevant mass fraction folders
        for xx_a = 1:length(mass_frac_folders)
            
            temp = cell2mat(mass_frac_folders(xx_a));
            if length(temp) >= 5
                if strcmp(temp(1:5),'Yeth=') == 1
                    log_mat_Y(xx_a) = 1;
                else
                    log_mat_Y(xx_a) = 0;
                end
            else
                log_mat_Y(xx_a) = 0;
            end
        end
        log_mat_Y = logical(log_mat_Y);
        mass_frac_folders = mass_frac_folders(log_mat_Y);
        
        % Sweep through all mass fraction folders
        for kk_a = 1:length(mass_frac_folders)
            
            % Check if variables exists. 
            variables_2_import = dir(join([overarching_folder_loc,cell2mat(propellant_folders(ii_a)),'\',cell2mat(diffusivity_folders(jj_a)),'\',cell2mat(mass_frac_folders(kk_a))],''));
            variables_2_import = struct2cell(variables_2_import);
            variables_2_import = variables_2_import(1,:);
            import_check = 0;
            
            % If so import the variable data file
            for ll = 1:length(variables_2_import)
                
                temp = cell2mat(variables_2_import(ll));
                if length(temp) >= 5
                    if strcmp(temp(1:5),'varia') == 1
                        import_check = 1;
                        break
                    end
                end
            end
            
            % Save away data into a cell
            if import_check == 1
                
                % Print filename and load data
                sprintf('%s\n',join([overarching_folder_loc,cell2mat(propellant_folders(ii_a)),'\',cell2mat(diffusivity_folders(jj_a)),'\',cell2mat(mass_frac_folders(kk_a)),'\',cell2mat(variables_2_import(ll))],''))
                load(join([overarching_folder_loc,cell2mat(propellant_folders(ii_a)),'\',cell2mat(diffusivity_folders(jj_a)),'\',cell2mat(mass_frac_folders(kk_a)),'\',cell2mat(variables_2_import(ll))],''))
                
                % Reset current cell used to save data 
                DATA_CELL_TEST = {0};
                
                % put propellant name in teh highest folder 
                DATA_MAT{1,ii_a} = propellant_folders(ii_a);
                
                % put names above data
                DATA_CELL_TEST{1,1} = {'r_droplet'};
                DATA_CELL_TEST{1,2} = {'U_dropet'};
                DATA_CELL_TEST{1,3} = {'T_droplet_surface'};
                DATA_CELL_TEST{1,4} = {'T_droplet_mean'};
                DATA_CELL_TEST{1,5} = {'m_droplet_eth'};
                DATA_CELL_TEST{1,6} = {'m_droplet_HFA'};
                DATA_CELL_TEST{1,7} = {'Y_liquid_eth_surface'};
                DATA_CELL_TEST{1,8} = {'Y_liquid_HFA_surface'};
                DATA_CELL_TEST{1,9} = {'Y_liquid_eth_mean'};
                DATA_CELL_TEST{1,10} = {'Y_liquid_HFA_mean'};
                DATA_CELL_TEST{1,11} = {'f_eth'};
                DATA_CELL_TEST{1,12} = {'f_HFA'};
                DATA_CELL_TEST{1,13} = {'time'};
                DATA_CELL_TEST{1,14} = {'user_propellant_name'};
                DATA_CELL_TEST{1,15} = {'D_liquid_mean'};
                DATA_CELL_TEST{1,16} = {'current_iter_initial_Y_eth_l'};
                
                % save data from a single test
                DATA_CELL_TEST{2,1} = r_droplet;
                DATA_CELL_TEST{2,2} = U_droplet;
                DATA_CELL_TEST{2,3} = T_droplet_surface;
                DATA_CELL_TEST{2,4} = T_droplet_mean;
                DATA_CELL_TEST{2,5} = m_droplet_eth;
                DATA_CELL_TEST{2,6} = m_droplet_HFA;
                DATA_CELL_TEST{2,7} = Y_liquid_eth_surface;
                DATA_CELL_TEST{2,8} = Y_liquid_HFA_surface;
                DATA_CELL_TEST{2,9} = Y_liquid_eth_mean;
                DATA_CELL_TEST{2,10} = Y_liquid_HFA_mean;
                DATA_CELL_TEST{2,11} = f_eth;
                DATA_CELL_TEST{2,12} = f_HFA;
                DATA_CELL_TEST{2,13} = time;
                DATA_CELL_TEST{2,14} = user_propellant_name;
                DATA_CELL_TEST{2,15} = D_liquid_mean;
                DATA_CELL_TEST{2,16} = current_iter_initial_Y_eth_l;
                
                % put headings and data in teh cell for teh particular
                % refrigerant
                TOT_DATA_MAT{1,1} = propellant_folders(ii_a);
                TOT_DATA_MAT{1,1+kk_a} = cell2mat(mass_frac_folders(kk_a));
                TOT_DATA_MAT{1+jj_a,1} = cell2mat(diffusivity_folders(jj_a));
                TOT_DATA_MAT{1+jj_a,1+kk_a} = DATA_CELL_TEST;
                
            end
            
            
        end
        
        % Save away the matrix for Y and D when we move to next propellant
        if jj_a == length(diffusivity_folders)
            DATA_MAT{2,ii_a} = TOT_DATA_MAT;
            TOT_DATA_MAT = {};
        end
        
    end
end

% clear all data
clearvars -except DATA_MAT