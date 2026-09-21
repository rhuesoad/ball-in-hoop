# MATLAB — model, control and simulation

## Requirements

- MATLAB (tested on R2025b) with the Control System Toolbox
- [CasADi](https://web.casadi.org/) 3.7.2 on the path, for the trajectory
  planner used by T3 and T4

## Layout

| Folder | Content |
|--------|---------|
| `config/` | Physical parameters (`ball_hoop_params.m`) and solver settings (`sim_config.m`) |
| `dynamics/` | Hybrid model: rolling on each hoop, free flight, impacts, motor inner loop |
| `control/` | LQR, TVLQR and the T4 controller |
| `trajectory/` | Direct-collocation planner, and the export of plans to the bench |
| `scenarios/` | One function per scenario, open loop and T1 to T4 |
| `analysis/` | Run I/O, plotting, animation, verification scripts |
| `results/` | One example run per scenario, and the run format (`RUN_SCHEMA.md`) |

## Running

From this folder:

```matlab
addpath(genpath(pwd))
T1_stabilization_10deg()      % any scenario, by name
run_all_scenarios()           % all 21
record_animations()           % one MP4 per saved run, into results/animations
```

Each run is saved under `results/runs/<date>/`.

## Sending a plan to the bench

T3 and T4 follow a trajectory planned here. `export_tvlqr_for_bench` and
`export_t4_for_bench` write it where `python/tasks/` reads it.
