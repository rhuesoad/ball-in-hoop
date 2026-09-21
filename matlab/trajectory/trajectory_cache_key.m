function key = trajectory_cache_key(x0_red, xf_red, mode, params, plan_opts)
%TRAJECTORY_CACHE_KEY  Deterministic hash of every input plan_trajectory_casadi.m's
%                      solve actually depends on, for cache/trajectories/<key>.mat.
%
%   key = trajectory_cache_key(x0_red, xf_red, mode, params, plan_opts)
%
%   Inputs
%   ------
%   x0_red, xf_red : boundary reduced states (4x1) [rad;rad/s;rad;rad/s]
%   mode           : dynamic mode string
%   params         : physical parameters (ball_hoop_params()); every
%                    numeric field is included, so replanning after any
%                    physical-parameter edit misses the cache rather than
%                    silently reusing a stale trajectory
%   plan_opts      : the SAME options struct handed to
%                    plan_trajectory_casadi.m -- every field is hashed,
%                    whatever it is called
%
%   Outputs
%   -------
%   key : 32-character hex MD5 digest, safe to use as a filename
%
%   WHY plan_opts IS TAKEN WHOLE rather than field by field. This function
%   used to list the solve inputs positionally (K, Tf_min, Tf_max, u_max,
%   contact_margin, psi_bounds), which meant every new planner option had
%   to be remembered in two places. It was not: constrain_theta_f reached
%   validate_scenario.m and docs/MODEL.md but neither the cache key nor
%   plan_trajectory_casadi.m's actual argument, so a scenario setting it
%   would have been handed a cached trajectory planned without it and
%   never known. Hashing the struct the solver is given removes the class
%   of bug rather than the instance: a field that reaches the planner
%   reaches the key by construction.
%
%   Note: uses java.security.MessageDigest (standard MATLAB, part of the
%   JVM every desktop MATLAB ships with). Octave's availability of the
%   java package is platform-dependent -- this is the one piece of the
%   scenario architecture not guaranteed to run under Octave.

    numeric_payload = [flatten_struct(params); flatten_struct(plan_opts); ...
                       x0_red(:); xf_red(:)];
    payload_bytes   = [typecast(numeric_payload.', 'uint8'), uint8(mode), ...
                       source_bytes(solve_dependencies())];

    md5    = java.security.MessageDigest.getInstance('MD5');
    digest = md5.digest(payload_bytes);   % returns Java byte[], comes back as int8
    key    = sprintf('%02x', typecast(digest, 'uint8'));
end


function files = solve_dependencies()
%SOLVE_DEPENDENCIES  The source files whose contents determine what the
%                    planner returns for a given set of inputs.
%
%   WHY THE CODE IS PART OF THE KEY. Hashing the physical parameters
%   catches "the plant changed"; it does not catch "the equations changed".
%   Editing rolling_matrices.m, hoop_geometry.m or the planner itself
%   leaves every hashed input identical, so a trajectory solved against
%   the OLD equations is handed back and used as though it came from the
%   new ones -- and it is used to produce headline results.
%
%   This is not hypothetical. T4b_flying_ball.m's docstring records a
%   successful catch; replanning it from a cold cache at the same commit
%   lands the ball on the inner hoop's outer surface instead
%   (rolling_in_outside at t = 1.2799 s), i.e. the published result came
%   from a cached trajectory the committed source no longer reproduces.
    files = { ...
        'plan_trajectory_casadi.m', ...
        'cascade_dynamics_reduced.m', ...
        'rolling_matrices.m', ...
        'hoop_geometry.m', ...
        'contact_forces.m'};
end


function bytes = source_bytes(files)
%SOURCE_BYTES  Raw contents of each dependency, concatenated.
%
%   Fails loudly on a missing file: a dependency silently dropped from the
%   key is exactly the failure this function exists to prevent, so it must
%   not degrade quietly into "hash what you could find".
    parts = cell(1, numel(files));
    for i = 1:numel(files)
        path = which(files{i});
        if isempty(path)
            error('trajectory_cache_key:missing_dependency', ...
                ['Cannot hash ''%s'': not on the MATLAB path. The trajectory cache key ' ...
                 'must cover every file the solve depends on, so this is fatal rather ' ...
                 'than skipped.'], files{i});
        end
        fid = fopen(path, 'r');
        if fid < 0
            error('trajectory_cache_key:unreadable_dependency', ...
                'Cannot read ''%s'' for hashing.', path);
        end
        parts{i} = fread(fid, Inf, '*uint8').';
        fclose(fid);
    end
    bytes = [parts{:}];
end


function v = flatten_struct(s)
%FLATTEN_STRUCT  Every field of s as one double column, field order fixed
%                by sorting the names so it cannot depend on the order the
%                caller happened to assign them in.
%
%   The field NAMES are hashed alongside their values: without them,
%   renaming an option or swapping two same-valued fields would collide.
    names = sort(fieldnames(s));
    parts = cell(numel(names), 1);
    for i = 1:numel(names)
        value = s.(names{i});
        if isstruct(value)
            % Nested option structs -- opts.liftoff is one -- are flattened
            % recursively rather than skipped. Skipping them would let two
            % genuinely different plans (different release windows, say)
            % share a cache entry, which is exactly the failure this
            % function exists to prevent.
            value = flatten_struct(value);
        elseif ischar(value)
            value = double(value);
        elseif islogical(value)
            value = double(value);
        end
        parts{i} = [double(names{i}).'; double(value(:))];
    end
    v = vertcat(parts{:});
end
