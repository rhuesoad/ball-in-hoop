%% Verification of stabilisability and detectability, task T1
% Runs over the two relevant equilibria of tab: equilibria.

modes = {'rolling_out', 'rolling_in_inside', 'rolling_in_outside'};

for m = 1:numel(modes)
    mode = modes{m};

    % --- Build A and B at the equilibrium of this mode -------------------
    % params must supply m_b, g, R_eff, and the entries M21, M22, C21, C22
    % of the rolling-mode matrices for the mode considered.
    p    = get_params(mode);
    psi_eq = p.psi_eq;                       % 0 or pi

    k_g = p.m_b * p.g * p.R_eff * cos(psi_eq);

    A = [ 0   1              0            0           ;
          0   0              0            0           ;
          0   0              0            1           ;
          0  -p.C21/p.M22   -k_g/p.M22   -p.C22/p.M22 ];

    B = [ 0 ; 1 ; 0 ; -p.M21/p.M22 ];

    n = size(A,1);

    % --- Controllability, and stabilisability if it fails ---------------
    Co     = ctrb(A,B);
    rank_C = rank(Co);

    if rank_C == n
        stab = true;                          % controllable, hence stabilisable
    else
        % PBH test restricted to the unstable and marginally stable modes.
        lam  = eig(A);
        stab = true;
        for k = 1:numel(lam)
            if real(lam(k)) >= 0
                if rank([A - lam(k)*eye(n), B]) < n
                    stab = false;
                end
            end
        end
    end

    % --- Detectability of the pair (A, Q^{1/2}) -------------------------
    Qh     = sqrtm(p.Q);                      % Q = diag(...) from eq: bryson
    Ob     = obsv(A,Qh);
    rank_O = rank(Ob);

    if rank_O == n
        det_ok = true;
    else
        lam    = eig(A);
        det_ok = true;
        for k = 1:numel(lam)
            if real(lam(k)) >= 0
                if rank([A - lam(k)*eye(n); Qh]) < n
                    det_ok = false;
                end
            end
        end
    end

    % --- Report ---------------------------------------------------------
    fprintf('--- %s (psi_eq = %.0f deg) ---\n', mode, rad2deg(psi_eq));
    fprintf('  rank(ctrb) = %d / %d   -> stabilisable : %d\n', rank_C, n, stab);
    fprintf('  rank(obsv) = %d / %d   -> detectable   : %d\n', rank_O, n, det_ok);
    fprintf('  eig(A)     = %s\n', mat2str(round(eig(A).',4)));
end