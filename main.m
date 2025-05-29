clear all
clc
close all
% ToDO!
% add gaussian integration points over the element so in your element 
addpath(genpath('visuals'))
addpath(genpath('solver'))
addpath(genpath('io'))
addpath(genpath('math'))
%%
Inputs = readCsvTables('Inputs.csv');

% -----------------------------------------------------
%               General structure properties
% -----------------------------------------------------
rho = Inputs.Properties.rho;
E= Inputs.Properties.E;
nu= Inputs.Properties.nu;

% -----------------------------------------------------
%                       Soil
% -----------------------------------------------------
k_s = Inputs.SoilProperties.slope; % N/m
L_soil = Inputs.SoilProperties.L_soil;
b_clamp = Inputs.SoilProperties.clamp;

% -----------------------------------------------------
%        Top Mass and Mass Moment of Inertia
% -----------------------------------------------------

M_lump = Inputs.TopMass.m;  % Lumped mass (kg)
I_lump = [Inputs.TopMass.i_xx, Inputs.TopMass.i_yy, Inputs.TopMass.i_zz]; % Mass moment of inertia [I_x, I_y, I_z] (kg·m²)
% Define the lumped mass position
X_offset = Inputs.TopMass.x_offset; % X offset (m)
Y_offset = Inputs.TopMass.y_offset; % Y offset (m)
Z_offset = Inputs.TopMass.z_offset; % Z offset (m)
b_TopMass = Inputs.TopMass.consdier_topmass;

% ----------------------------------------------------
%     General script 
% ----------------------------------------------------
spacing = Inputs.General.spacing; % m
FAK = Inputs.General.FAK; % from meter to mili meter
b_Write_stfile = Inputs.General.b_export_st_file;

% ----------------------------------------------------
%     	Can geometry 
% ----------------------------------------------------
geometry = [Inputs.Beam.sec_no,Inputs.Beam.z_bot, ...
    Inputs.Beam.z_top, Inputs.Beam.d_bot, ... 
    Inputs.Beam.d_top, Inputs.Beam.thick];       % Top

geometry = geometry/FAK;
comb_st  = table(geometry(:,1),geometry(:,2),geometry(:,3), geometry(:,4),geometry(:,5),geometry(:,6),...
                'Variablenames', ... 
                {'sec',  'zTop','zBot',   'dTop',   'dBotm',   'Thick'});
            
% constructuring the Timoshenko general beam structure
results = make_timoshenko_beam(comb_st, spacing, rho, E, nu);    

% visualize the structure with its boundary conditions
visualize_structure(comb_st, b_TopMass, ... 
    L_soil, max(Inputs.Beam.d_bot)/FAK, b_clamp)

%% Formatting and writing an st file
% Display results in a table
if b_Write_stfile
    write_st_file('tower_hossein.st', results);
    disp(['Structure mass : ', num2str(trapz(results.H,results.m)/FAK), ' [Ton]'])
end

%% Euler FEM 
[eigenvectors, eigenvalues] = euler_beam_solver(results,b_clamp,b_TopMass , ... 
    X_offset,Y_offset,Z_offset,M_lump,I_lump, ... 
    L_soil,k_s);

frequencies = sqrt(diag(eigenvalues)) / (2 * pi);  % Convert to Hz
% Sort frequencies in ascending order
frequencies = sort(frequencies);

% Output first few natural frequencies
fprintf('First 18 natural frequencies (Hz):\n');
for i=1:18
    disp(['Mode ',num2str(i),' ',num2str(frequencies(i)), ' Hz']);
end

