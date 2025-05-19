clear all
clc
close all
global rho
% ToDO!
% add gaussian integration points over the element so in your element 
addpath(genpath('visuals'))
addpath(genpath('solver'))
addpath(genpath('io'))
%%

% -----------------------------------------------------
%               General structure properties
% -----------------------------------------------------
rho = 7850;
E= 2.1e11;
nu= 0.3;

% -----------------------------------------------------
%                       Soil
% -----------------------------------------------------
k_s = 1e7; % N/m
L_soil = 26.1;
b_clamp = true;

% -----------------------------------------------------
%        Top Mass and Mass Moment of Inertia
% -----------------------------------------------------

M_lump = 50000;  % Lumped mass (kg)
I_lump = [0, 0, 0]; % Mass moment of inertia [I_x, I_y, I_z] (kg·m²)
% Define the lumped mass position
X_offset = 160; % X offset (m)
Y_offset = 0.5; % Y offset (m)
Z_offset = 0.0; % Z offset (m)
b_TopMass = true;

% ----------------------------------------------------
%     General script 
% ----------------------------------------------------
spacing = 0.1; % m
FAK = 1000; % from meter to mili meter
b_Write_stfile = false;

% ----------------------------------------------------
%     	Can geometry 
% ----------------------------------------------------
geometry = [1	0	    4800	4200	4200	76        % Bottom
            2	4800	6800	4200	4200	60
            3	6800	8800	4200	4200	52
            4	8800	10800	4200	4200	46
            5	10800	12800	4200	4200	42
            6	12800	14800	4200	4200	39
            7	14800	16800	4200	4200	37
            8	16800	19800	4200	4200	35
            9	19800	22800	4200	4200	33
            10	22800	25800	4200	4200	31
            11	25800	28800	4200	4200	30
            12	28800	31800	4200	4200	29
            13	31800	34800	4200	4200	28
            14	34800	37800	4200	4200	27
            15	37800	42800	4200	4200	26
            16	42800	45800	4200	4200	24
            17	45800	48800	4200	4200	23
            18	48800	51800	4200	4200	22
            19	51800	61800	4200	3900	22
            20	61800	70800	3900	3630	20
            21	70800	81800	3630	3300	19
            22	81800	83800	3300	3240	20
            23	83800	85800	3240	3180	21
            24	85800	87800	3180	3120	23];       % Top


% geometry = [1	0	        150000	    8000	8000	20];

geometry = geometry/FAK;
comb_st  = table(geometry(:,1),geometry(:,2),geometry(:,3), geometry(:,4),geometry(:,5),geometry(:,6),...
                'Variablenames', ... 
                {'sec',  'zTop','zBot',   'dTop',   'dBotm',   'Thick'});
            
% constructuring the Timoshenko general beam structure
results = make_timoshenko_beam(comb_st, spacing, rho, E, nu);    

% visualize the structure with its boundary conditions
visualize_structure(comb_st, b_TopMass, b_clamp)

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

%% functions


function [mass,mass_Over_L, cog,Larr] = get_section_properties(zTop, zBot,dTop, dBot, thick)
    global rho
    % Calculate section length (convert from mm to m)
    dH = 0.1;
    L = (zBot-zTop); 
    dL = ceil(L/dH);
    Larr = linspace(zTop,zBot,dL+1);
%     Larr = Larr(2:end);
    mass_can = calculate_can_mass(L,dBot,dTop,thick);
    disp(['Overall can mass : ' , num2str(mass_can), ' [kg]']);
    for i =1:length(Larr)-1
        
        dBotLoc = interp1([zTop,zBot], [dTop,dBot ],Larr(i))/2; % radius
        dTopLoc = interp1([zTop,zBot], [dTop,dBot ],Larr(i+1))/2;% radius
        
        V_outer = (1/3) * pi* dH * (dBotLoc^2 + dBotLoc*dTopLoc + dTopLoc^2);
        V_in = (1/3) * pi* dH * ((dBotLoc-thick)^2 + (dBotLoc-thick)*(dTopLoc-thick) + (dTopLoc-thick)^2);

        % Calculate average outer and inner diameters (convert from mm to m)
    % %     D_outer_avg = (dBot + dTop) / 2; % Average diameter in meters
    % %     D_inner_avg = D_outer_avg - 2 * (thick ); % Inner diameter

        % Cross-sectional area of a hollow cylinder
    % %     cross_sectional_area = pi / 4 * (D_outer_avg^2 - D_inner_avg^2);

        % Volume and mass calculation
    % %     volume = cross_sectional_area * L;
        volume = V_outer - V_in;
        mass = volume * rho;
        mass_Over_L(i) = mass/dH;
        cog(i) = get_cog(zTop, zBot,dTop, dBot);
    end
    
    mass_Over_L(end+1) = mass_Over_L(end);
    
    disp(['Back integrated mass from the discritized can : ', num2str(trapz(Larr, mass_Over_L)), ' [kg]'])
end

function cog = get_cog(zTop, zBot,dTop, dBot)
    h = zBot - zTop;
    cog = (h/3 )*((2*dBot+3)/(dBot+dTop)) ;
    
end

function mmi = get_mass_moment_of_inertia(cog, mass)
    mmi = cog^2*mass;
end

function [A,Din,Dout] = area_calculation(comb_st,can,hCog)
    Dout = interp1([comb_st.zTop(can), comb_st.zBot(can)], ...
                   [comb_st.dTop(can), comb_st.dBotm(can)],... 
                   hCog);
    Din  = interp1([comb_st.zTop(can), comb_st.zBot(can)], ... 
                   [comb_st.dTop(can), comb_st.dBotm(can)], ... 
                   hCog)- 2*comb_st.Thick(can);
   
               
    A = pi.*(Dout.^2)./4 - pi.*(Din.^2)./4;
    
end

function [Ix, Iy] = calculate_area_moment_of_inertia(Dout,Din)
    Ix  = pi .* (Dout.^4.0 - Din.^4.0) ./ 64.0;
    Iy = Ix;
end




function drawCenterOfGravitySymbol(radius, centerX, centerY)
    % Function to draw a center of gravity symbol
    % Inputs:
    % - radius: radius of the circle
    % - centerX: x-coordinate of the circle center
    % - centerY: y-coordinate of the circle center
    
    % Define angles for the circle's quarters
    theta1 = linspace(0, pi/2, 100);       % First quadrant
    theta2 = linspace(pi/2, pi, 100);      % Second quadrant
    theta3 = linspace(pi, 3*pi/2, 100);    % Third quadrant
    theta4 = linspace(3*pi/2, 2*pi, 100);  % Fourth quadrant
    
    % Compute the coordinates of each quarter
    x1 = radius * cos(theta1) + centerX;
    y1 = radius * sin(theta1) + centerY;
    x2 = radius * cos(theta2) + centerX;
    y2 = radius * sin(theta2) + centerY;
    x3 = radius * cos(theta3) + centerX;
    y3 = radius * sin(theta3) + centerY;
    x4 = radius * cos(theta4) + centerX;
    y4 = radius * sin(theta4) + centerY;
    
    % Draw the quarters
    fill([centerX, x1, centerX], [centerY, y1, centerY], 'k'); % First quadrant (dark)
    hold on;
    fill([centerX, x2, centerX], [centerY, y2, centerY], 'w'); % Second quadrant (white)
    fill([centerX, x3, centerX], [centerY, y3, centerY], 'k'); % Third quadrant (dark)
    fill([centerX, x4, centerX], [centerY, y4, centerY], 'w'); % Fourth quadrant (white)
    
    % Add circle boundary
    theta = linspace(0, 2*pi, 200);
    x = radius * cos(theta) + centerX;
    y = radius * sin(theta) + centerY;
    plot(x, y, 'k', 'LineWidth', 1.5); hold on
    yline(centerY, 'r--', 'LineWidth', 1.5); % Dashed red line at x = 0
    
    % Set axis properties
    axis equal;
    grid on;
    xlabel('x-axis');
    ylabel('y-axis');
    hold off;
end



function frequencies = calcFixedFixedBeamFreq(E, rho, R_o, R_i, L)
    % Function to calculate the first 4 natural frequencies of a fixed-fixed hollow cylindrical beam
    % Inputs:
    % - E: Young's modulus (Pa)
    % - rho: Density of the material (kg/m^3)
    % - R_o: Outer radius of the hollow cylinder (m)
    % - R_i: Inner radius of the hollow cylinder (m)
    % - L: Length of the beam (m)
    % Output:
    % - frequencies: Vector containing the first 4 natural frequencies (Hz)
    
    % Mode constants for fixed-fixed beam
    beta = [4.730, 7.853, 10.996, 14.137]; % First 4 mode constants
    
    % Second moment of area (I) and mass per unit length (m)
    I = (pi / 4) * (R_o^4 - R_i^4);               % Second moment of area
    m = rho * pi * (R_o^2 - R_i^2);               % Mass per unit length
    
    % Precompute the factor
    factor = sqrt(E / rho) * sqrt(I / m) / (2 * pi * L^2);
    
    % Calculate the frequencies
    frequencies = beta.^2 * factor;
end



function mass = calculate_can_mass(L,dBot,dTop,thick)
    global rho
    rBot = dBot/2;
    rTop = dTop/2;
        V_outer = (1/3) * pi* L * (rBot^2 + rBot*rTop + rTop^2);
        V_in = (1/3) * pi* L * ((rBot-thick)^2 + (rBot-thick)*(rTop-thick) + (rTop-thick)^2);
        dV = V_outer - V_in;
        mass = dV * rho;
end


function results = make_timoshenko_beam(comb_st, spacing, rho, E, nu)

    Zint= linspace(min(comb_st.zTop), ... 
                   max(comb_st.zBot),...
                   floor((max(comb_st.zBot)-min(comb_st.zTop))/spacing));
    initialization_array = zeros(1,length(Zint));

    Zarr = initialization_array;
    Marr = initialization_array;
    Aarr = initialization_array;
    Ixarr = initialization_array;
    Iyarr = initialization_array;
    rxarr = initialization_array;
    ryarr = initialization_array;
    Garr = initialization_array;
    Jarr = initialization_array;
    HcogArr = initialization_array;
    MassArr = initialization_array;
    flag = false;

    for i =1:length(Zint)
        hTop = comb_st.zTop;
        hBot = comb_st.zBot;
        idx1 = find(Zint(i)>=hTop,1,'last');
        data_can.zTop(1) = hTop(idx1);
        data_can.zBot(1) = hBot(idx1);
        data_can.dBotm(1) = comb_st.dBotm(idx1);
        data_can.dTop(1) =  comb_st.dTop(idx1);
        data_can.Thick(1) =  comb_st.Thick(idx1);

        [A,Din,Dout] = area_calculation(data_can,1,Zint(i));
         mL = A *rho;
        [Ix, Iy] = calculate_area_moment_of_inertia(Dout,Din);
         G = E / (2 * (1 + nu));
         rix = (Ix ./ A) .^ 0.5;
         riy = (Iy ./ A) .^ 0.5;  

        Zarr(i) = Zint(i);
        Marr(i) = mL;
        Aarr(i)  = A;
        Ixarr(i) = Ix;
        Iyarr(i) = Iy;
        rxarr(i) = rix;
        ryarr(i) = riy;
    end

    zero_array  = repmat(0, length(Zarr), 1);

    results = table(Zarr', ...
                    Marr', ...
                    zero_array,zero_array, ...
                    rxarr',ryarr', zero_array,zero_array, ...
                    repmat(E, length(Zarr), 1), ... 
                    repmat(G, length(Zarr), 1), ... 
                    Ixarr', Iyarr', ...
                    2*Ixarr', ...
                    repmat(2*(1+nu)/(4+3*nu), length(Zarr), 1), repmat(2*(1+nu)/(4+3*nu), length(Zarr), 1), ...
                    Aarr', ...
                    zero_array, zero_array , zero_array , ...
                    'VariableNames', {'H', 'm', 'x_cg', 'y_cg', 'ri_x', 'ri_y', 'x_sh', 'y_sh', ...
                                      'E', 'G', 'I_x', 'I_y', 'J', 'Kx', 'Ky', 'A', 'pitch', 'xe', 'ye'});
end
