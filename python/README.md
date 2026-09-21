# Python — bench

Runs on the Raspberry Pi connected to the ODrive and the camera.

## Requirements

Python 3 with `numpy`, `scipy`, `opencv-python`, `odrive`, and `picamera2`
(Raspberry Pi only).

## Layout

| Folder | Content |
|--------|---------|
| `common/` | Ball detection, control loop and logging, estimators, hoop geometry |
| `experiments/` | Identification experiments E1 to E4 |
| `tasks/` | Control tasks T1 to T4: one config and one runner each, and the plans exported from MATLAB |
| `tools/` | Camera calibration and capture, ODrive calibration and tuning, one-off bench checks |

| Experiment | Measures |
|------------|----------|
| E1 | Motor and hoop inertia and friction |
| E2 | Ball friction, from free oscillations |
| E3 | The whole acquisition chain, from a hoop step |
| E4 | Static friction between ball and track |

## Running

From this folder, as modules:

```bash
python3 -m tasks.t1_run --psi-init 45
T1_HOOP=inner python3 -m tasks.t1_run --psi-init 45
python3 -m experiments.e2_run
python3 -m tools.odrive_calibration
```

Each run writes a `.npz` capture.
