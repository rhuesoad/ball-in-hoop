# The `run` struct schema

Single source of truth for the field names and units of the `run` struct produced
by `build_run.m`, persisted by `save_run.m`/`load_run.m`, and consumed by every
`plot_*` function in `results/visualization/`. `run_scenario.m` builds one for
every scenario, open-loop or closed-loop (docs/MODEL.md sec. 9); an open-loop
scenario's `run` has `ref_psi_rad`/`ref_psi_dot_rad_s`/`error_rad` all-NaN
(there is no reference to track) and `u_rad_s2` all-zero (no controller).
`run_scenario.m` also adds fields beyond this schema for some scenario kinds
(`eig_cl`, `solver_info`, `mode_sequence_check`) — see docs/MODEL.md sec. 9.7,
the source of truth for those rather than a duplicate list here. `source`,
`calibration`, and `ode_segments` (below) are the simulation-vs-hardware
fields — see docs/MODEL.md sec. 10 and docs/VALIDATION.md for how they're
used, not just what they are.

Unit suffix convention: `_Nm` (N·m), `_rad_s2` (rad/s²), `_rad_s` (rad/s), `_rad`
(rad), `_s` (s). A field with no suffix is unitless or not a physical quantity
(strings, structs, timestamps).

## Fields

| Field         | Type              | Units    | Description |
|---------------|-------------------|----------|-------------|
| `name`        | char              | —        | Scenario name (used for filenames and plot titles) |
| `source`      | char              | —        | `'simulation'` (set by `build_run.m`) or `'experiment'` (set by `import_experiment.m`); every consumer that treats the two differently switches on this field, never on which function produced the struct |
| `t_s`         | `N×1 double`      | s        | Time vector. `t_s(1) = 0` is defined by an explicit event (docs/MODEL.md sec. 10.3) — never by eye |
| `ode_segments` | cell array of ODE solution structs | — | **Simulation only.** One continuous solution struct per `scenario`-run mode interval (same segmentation as `sol.mode_intervals`, dynamics/ball_hoop_ode.m), each usable directly with `deval`. This is what lets `compare_runs.m` resample the simulation onto the experiment's time grid without ever calling `interp1` on discrete samples (docs/MODEL.md sec. 10.4) |
| `calibration` | struct            | —        | **Experiment only.** Everything needed to map raw sensor output to the state vector — see docs/MODEL.md sec. 10.2 for the field list |
| `x`           | `N×6 double`      | mixed    | State trajectory, one row per `t_s(i)`. Columns match the notation contract's state vector exactly — see table below |
| `tau_Nm`      | `N×1 double`      | N·m      | Applied motor torque, already saturated to `±cfg.tau_max` (the only saturation point is `motor_inner_loop.m`) |
| `u_rad_s2`    | `N×1 double`      | rad/s²   | Controller acceleration setpoint (`u = theta_ddot_ref`) at each `t_s(i)`, before the `u -> tau` conversion |
| `ref_psi_rad` | `N×1 double`      | rad      | Reference `psi_ref(t)` evaluated along `t_s` |
| `ref_psi_dot_rad_s` | `N×1 double` | rad/s    | Reference `psi_dot_ref(t)` evaluated along `t_s` |
| `error_rad`   | `N×1 double`      | rad      | Tracking error, `ref_psi_rad - x(:,5)` |
| `z_ctrl`      | `N×n_ctrl double` | mixed    | Controller internal state trajectory (e.g. the PID integral term); meaning depends on `ctrl_type` |
| `scenario`    | struct            | —        | Full scenario struct as defined in `scenarios/**/*.m` (docs/MODEL.md sec. 9), kept for reproducibility |
| `params`      | struct            | —        | Physical parameters used for the run (see `ball_hoop_params.m`), kept for reproducibility |
| `cfg`         | struct            | —        | Numerical settings used for the run (see `sim_config.m`), kept for reproducibility |
| `timestamp`   | `datetime`        | —        | Creation time, used to build the `save_run.m` filename |

### `x` columns (notation contract, `ball_hoop_params.m`/dynamics files)

| Column | Symbol      | Units | Description |
|--------|-------------|-------|-------------|
| 1      | `r`         | m     | Ball-centre radial distance from the hoop axis |
| 2      | `r_dot`     | m/s   | Radial velocity |
| 3      | `theta`     | rad   | Hoop (motor) angle |
| 4      | `theta_dot` | rad/s | Hoop angular velocity |
| 5      | `psi`       | rad   | Ball angular position, inertial frame, from the downward vertical |
| 6      | `psi_dot`   | rad/s | Ball angular velocity |

## Producers / consumers that must obey this schema

- `results/visualization/build_run.m` — the only place a `run` struct is assembled.
- `results/io/save_run.m` / `results/io/load_run.m` — persist/restore it unchanged.
- `results/visualization/plot_closed_loop_results.m` — reads `t_s`, `x`, `tau_Nm`,
  `ref_psi_rad`, `error_rad`, `cfg.tau_max`, `name`.
- `results/visualization/plot_runs.m` — reads the same fields across multiple
  saved runs for comparison plots.
- `results/visualization/plot_run.m` — the single entry point (Phase 2 target
  item 6); dispatches to the lower-level `plot_*` helpers in
  `results/visualization/plots/`, unpacking `run` into the raw arrays those
  helpers already expect.

## History

Renamed from an earlier, undocumented schema (`t`, `X`, `tau`, `ref_psi`,
`ref_psi_dot`, `error`, `Z_ctrl`) that had no unit suffixes and was never written
down in one place — every consumer had to infer the field meanings from
`build_run.m`'s source. `u_rad_s2` is new: the controller's acceleration setpoint
was computed at every step but never logged, only the post-conversion torque was.

`source`, `ode_segments`, and `calibration` (Phase 6) exist so a simulation run
and an experiment run are the *same struct* — the explicit design goal being
that no plotting or metric function ever needs an "experiment version." See
docs/MODEL.md sec. 10.
