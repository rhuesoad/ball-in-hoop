%% Validate the cascade linearization for all rolling modes.
clear; close all; clc;
addpath(genpath(pwd));
params = ball_hoop_params();

modes = {
    'rolling_out',        0,    'stable';
    'rolling_in_outside', pi,   'inverted (unstable in open loop)';
    'rolling_in_inside',  0,    'stable';
};

for k = 1:size(modes, 1)
    mode_name = modes{k, 1};
    psi_eq    = modes{k, 2};
    nature    = modes{k, 3};
    
    fprintf('\n=== Mode: %s ===\n', mode_name);
    fprintf('Equilibrium     : psi_eq = %.4f rad (%.1f deg) [%s]\n', ...
            psi_eq, rad2deg(psi_eq), nature);
    
    [A, B, C_out, D] = linearize_system(psi_eq, params, mode_name);
    
    % --- Eigenvalues ---
    eigA = eig(A);
    fprintf('Eigenvalues     : ');
    disp(eigA.');
    
    % --- Controllability ---
    Co = ctrb(A, B);
    fprintf('Controllability rank : %d / 4\n', rank(Co));
    
    % --- Cascade signature ---
    fprintf('B(2,1) = %.6f (should be 1.0)\n', B(2));
end