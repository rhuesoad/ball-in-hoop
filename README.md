# Ball in a Double Hoop

Code, CAD and documents of the master's thesis *Flying Ball in a Hoop — Design
of New Didactic Devices for Teaching of Control Engineering* (Rafaël Hueso
Adelantado), electromechanical engineering (mechatronics), ULB, 2025–2026.

Supervisors: Prof. Michel Kinnaert and Laurent Catoire.

A ball rolls inside two concentric hoops, turned by a single motor. The inner
hoop has a 90° gap. Four control tasks of increasing difficulty are posed on it:

| Task | Goal | Simulation | Bench |
|------|------|------------|-------|
| T1 | Bring the ball back to the bottom and hold it there | works | works, both hoops |
| T2 | Tracking | works | works, both hoops |
| T3 | Loop the loop | works | the ball leaves the track near the top (95–104°) |
| T4 | Throw the ball through the gap and catch it inside | works | not ported yet |

The thesis (`docs/thesis.pdf`) gives the model, the controllers and the results.

## Layout

The project works with a simulation, made on **Matlab**, and real hardware built with 
3D printing, and controlled with a Raspberry Pi controller and **Python**. 

| Folder | Content |
|--------|---------|
| `matlab/` | Model, controllers, trajectory planning, scenarios, and one example run per scenario |
| `python/` | Code that runs on the Raspberry Pi and drives the bench |
| `cad/` | OpenSCAD source of every printed part, one file per part |
| `analysis/` | Notebook that produces the report figures from the bench captures |
| `docs/` | Thesis and defence |

Each folder has its own README.

## Hardware

ODrive S1 driver, D5065 motor, AMT102-V encoder, Raspberry Pi with a Pi Camera,
and the printed parts in `cad/`.

## Where to pick it up

- **T4 on the bench.** The simulation works because the hoop is steered during
  the flight so that the gap faces the ball when it crosses the inner hoop
  (`matlab/scenarios/closed_loop/T4/T4_flying_ball.m`). The bench code in
  `python/tasks/t4_*.py` predates this and still loads the old plan.
- **T3 on the bench.** The two kept trials show where the ball leaves the track.
- **The model tests the ball centre only** when deciding whether it passes the
  gap, not the whole ball, so a clean mode sequence in simulation does not
  guarantee the ball clears the edge.
- ** Hardware **: There is a play on the connection piece. Try to start from there
   to try to get better results on the looping itself. 

## License

No license yet. `cad/lib/nutsnbolts.scad` is a third-party library by Johannes
Kneer under GPLv3; its notice is kept in the file.
