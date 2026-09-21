function K_traj = tvlqr_design(t_traj, x_traj, u_traj, Q, R, Qf, params)
%TVLQR_DESIGN  Compute the time-varying LQR gain along a reference trajectory.
%
%   K_traj = tvlqr_design(t_traj, x_traj, u_traj, Q, R, Qf, params)
%
%   Linearizes the cascade dynamics along the trajectory (x*, u*) and
%   solves the differential Riccati equation backwards in time:
%       -dS/dt = A'*S + S*A - S*B*inv(R)*B'*S + Q,  S(Tf) = Qf
%   to produce the time-varying gain K(t) = inv(R)*B'(t)*S(t).
%
%   Inputs
%   ------
%   t_traj : time grid (N x 1)                                          [s]
%   x_traj : reference reduced state trajectory (N x 4)
%            [theta, theta_dot, psi, psi_dot], [rad, rad/s, rad, rad/s]
%   u_traj : reference control trajectory (N x 1)                       [rad/s^2]
%   Q, Qf  : LQR weight matrices (4x4), one entry per state component squared
%   R      : LQR control weight (scalar), per (rad/s^2) squared
%   params : physical parameters
%
%   Outputs
%   -------
%   K_traj : time-varying gain (N x 4) — K(t_k) for each knot k
%            [(rad/s^2)/state-unit, per column]
%
%   Reference: Gurtner & Zemanek (2017), eqs. 26-30; Bryson & Ho (1975),
%   ch. 5 (neighboring extremals).

    N   = length(t_traj);
    n_x = 4;
    
    %% --- Linearize along the trajectory ---
    A_traj = zeros(n_x, n_x, N);
    B_traj = zeros(n_x, 1, N);
    
    for k = 1:N
        [A_traj(:,:,k), B_traj(:,:,k)] = linearize_cascade( ...
            x_traj(k,:).', u_traj(k), params);
    end
    
    %% --- Solve differential Riccati equation backwards ---
    % Integrate from Tf (S=Qf) to t=0, requesting output exactly at the
    % trajectory knots (flipud(t_traj) is a decreasing sequence, which
    % ode45 accepts directly as the output times for backward integration).
    riccati_rhs = @(t, S_vec) riccati_dynamics(t, S_vec, t_traj, A_traj, B_traj, R, Q, n_x);
    [~, S_resampled] = ode45(riccati_rhs, flipud(t_traj), Qf(:));
    S_resampled = flipud(S_resampled);   % back to forward time order, matching t_traj

    S_traj = zeros(n_x, n_x, N);
    for k = 1:N
        S_traj(:,:,k) = reshape(S_resampled(k,:), n_x, n_x);
    end
    
    %% --- Compute K(t) = R^-1 * B'(t) * S(t) ---
    K_traj = zeros(N, n_x);
    for k = 1:N
        K_traj(k, :) = (B_traj(:,:,k).' * S_traj(:,:,k)) / R;
    end
end


%% ====================================================================
%  HELPERS
%% ====================================================================

function dS_vec = riccati_dynamics(t, S_vec, t_traj, A_traj, B_traj, R, Q, n_x)
%RICCATI_DYNAMICS  RHS of the differential Riccati equation.
%   -dS/dt = A'S + SA - SBR^-1B'S + Q
%   We integrate backwards, so dS/dt = -(A'S + SA - SBR^-1B'S + Q)
%   But ode45 integrates forward; we reverse time externally by giving
%   it a decreasing tspan.
%
%   Inputs
%   ------
%   t      : current time [s]
%   S_vec  : Riccati matrix S(t), vectorized (n_x^2 x 1)
%   t_traj, A_traj, B_traj : trajectory grid [s] and linearized matrices, see tvlqr_design
%   R, Q   : LQR weight matrices, see tvlqr_design
%   n_x    : number of reduced states (4)
%
%   Outputs
%   -------
%   dS_vec : d/dt of S(t), vectorized (n_x^2 x 1)


    S = reshape(S_vec, n_x, n_x);
    
    % Interpolate A(t) and B(t) at the current time
    A = interp_matrix(t, t_traj, A_traj);
    B = interp_matrix(t, t_traj, B_traj);
    
    dS = -(A.' * S + S * A - S * B * (B.' * S) / R + Q);
    
    % We integrate forward over a backward-going tspan, so the sign is correct
    % as written above.
    dS_vec = dS(:);
end


function M_t = interp_matrix(t, t_grid, M_traj)
%INTERP_MATRIX  Linear interpolation of a stack of matrices in time.
%
%   Inputs
%   ------
%   t      : query time [s]
%   t_grid : time grid the stack is sampled at (N x 1) [s]
%   M_traj : matrix stack (n x m x N)
%
%   Outputs
%   -------
%   M_t : interpolated matrix (n x m) at time t
    if t <= t_grid(1)
        M_t = M_traj(:,:,1);
    elseif t >= t_grid(end)
        M_t = M_traj(:,:,end);
    else
        % Find index
        idx_low = find(t_grid <= t, 1, 'last');
        idx_high = idx_low + 1;
        alpha = (t - t_grid(idx_low)) / (t_grid(idx_high) - t_grid(idx_low));
        M_t = (1-alpha) * M_traj(:,:,idx_low) + alpha * M_traj(:,:,idx_high);
    end
end


function [A, B] = linearize_cascade(x_eq, u_eq, params, mode)
%LINEARIZE_CASCADE  Jacobians of the cascade dynamics at (x_eq, u_eq).
%
%   Thin wrapper around linearize_system.m (docs/AUDIT.md item 1.7): the
%   cascade Jacobian only depends on psi (the dynamics are already
%   linear in theta, theta_dot and u, so their partial derivatives don't
%   depend on the operating point's theta_dot/u_eq values), so this is
%   exactly linearize_system evaluated at psi_eq = x_eq(3). An earlier
%   version hand-rederived M(2,2)/M(2,1)/C(2,1)/C(2,2) here independently
%   of rolling_matrices.m and got both wrong: it used the hoop contact
%   radius (R_o/R_b) where the ball's own-spin terms need the ball-centre
%   orbit radius (R_eff/R_b), and a "+1" where the coupling term needs a
%   "-1". Delegating to linearize_system.m (which already goes through
%   hoop_geometry.m -> rolling_matrices.m) removes the second,
%   independent copy of this derivation entirely.
%
%   mode defaults to 'rolling_out', the only mode this was ever called
%   with (tvlqr_design.m itself has no mode parameter yet).
%
%   Inputs
%   ------
%   x_eq   : reduced state at the operating point [theta, theta_dot, psi, psi_dot]
%            [rad, rad/s, rad, rad/s]
%   u_eq   : control at the operating point [rad/s^2]
%   params : physical parameters
%   mode   : 'rolling_out' | 'rolling_in_outside' | 'rolling_in_inside' (optional)
%
%   Outputs
%   -------
%   A : 4x4 state matrix
%   B : 4x1 input matrix
    if nargin < 4
        mode = 'rolling_out';
    end
    psi_eq = x_eq(3);
    [A, B, ~, ~] = linearize_system(psi_eq, params, mode);
end