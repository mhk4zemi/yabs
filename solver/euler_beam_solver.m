function [eigenvectors, eigenvalues]=euler_beam_solver(results,b_Clamp,b_TopMass, ...
                                                        X_offset,Y_offset,Z_offset,M_lump,I_lump, ... 
                                                        L_soil,k_s)

H = results.H;      % Node positions
m_kg = results.m;   % Mass per element
E_new = results.E;      % Young's modulus (constant)
Ix =results.I_x;     % Moment of inertia about x-axis (same as y in this case)
A = results.A;

[~, lumped_node] = min(abs(H - X_offset));
% Compute the DOFs associated with the node
node_dof = 2 * lumped_node - 1; % Displacement DOF
rotation_dof = 2 * lumped_node; % Rotation DOF


n_segments = length(H) - 1;
dof = 2 * (n_segments + 1);  % Degrees of freedom for each node (displacement, rotation)
K_global = zeros(dof, dof);
M_global = zeros(dof, dof);

x = linspace(0, length(H), n_segments); % N+1 nodes


% Assemble element matrices and add to global matrices
for i = 1:n_segments
    L = H(i+1) - H(i);  % Length of the element
    rho = m_kg(i) / (A(i));  % Density (kg/m^3)
    
    % Element stiffness matrix (for a beam element with bending)
    k = E_new(i) * Ix(i) / L^3 * [12 6*L -12 6*L;
                            6*L 4*L^2 -6*L 2*L^2;
                            -12 -6*L 12 -6*L;
                            6*L 2*L^2 -6*L 4*L^2];
    
    % Element mass matrix (consistent mass matrix for a beam element)
    m_elem = rho * A(i) * L / 420 * [156 22*L 54 -13*L;
                                     22*L 4*L^2 13*L -3*L^2;
                                     54 13*L 156 -22*L;
                                     -13*L -3*L^2 -22*L 4*L^2];
    
    % Assign to global matrices (two degrees of freedom per node)
    K_global(2*i-1:2*i+2, 2*i-1:2*i+2) = K_global(2*i-1:2*i+2, 2*i-1:2*i+2) + k;
    M_global(2*i-1:2*i+2, 2*i-1:2*i+2) = M_global(2*i-1:2*i+2, 2*i-1:2*i+2) + m_elem;
end

if b_TopMass
    % adding the lumped mass and MMI
    M_global(node_dof, node_dof) = M_global(node_dof, node_dof) + M_lump;
    % Update rotational inertia contributions
    M_global(rotation_dof, rotation_dof) = M_global(rotation_dof, rotation_dof) + I_lump(1);
    % Include Coupling Effects from Mass Offset (Rotational Contributions)
    % The offset creates additional moment contributions due to mass inertia coupling
    % These terms arise from m*r^2 effects in the rotational equations
    M_global(rotation_dof, rotation_dof) = M_global(node_rot_dof, node_rot_dof) + M_lump * (Y_offset^2 + Z_offset^2);

    % Cross coupling terms: Offset creates interaction between rotation and displacement DOFs
    M_global(node_dof, rotation_dof) = M_global(node_disp_dof, node_rot_dof) - M_lump * Y_offset;
    M_global(rotation_dof, node_dof) = M_global(node_rot_dof, node_disp_dof) - M_lump * Y_offset;
end


if b_Clamp
    % Apply boundary conditions (fixed at the first node, zero displacement and rotation)
    K_global(1:2, :) = [];  % Remove rows and columns corresponding to the first node's DOFs
    K_global(:, 1:2) = [];
    M_global(1:2, :) = [];
    M_global(:, 1:2) = [];
else
    % Find number of affected segments
    num_segments_soil = find(cumsum(H(2:end) - H(1:end-1)) <= L_soil, 1, 'last');
    idx = find(cumsum(H(2:end) - H(1:end-1)) <= L_soil);
    dx= diff(H(idx));
    % Add soil stiffness to the global stiffness matrix
    for i = 1:num_segments_soil
        K_global(2*i-1, 2*i-1) = K_global(2*i-1, 2*i-1) + k_s*dx(1);
    end
end

% Solve eigenvalue problem for natural frequencies
[eigenvectors, eigenvalues] = eig(K_global, M_global);
