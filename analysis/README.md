# Analysis — report figures

`bench_analysis.ipynb` turns the bench captures into the figures of the thesis.

## Requirements

Python 3 with `numpy`, `scipy`, `pandas`, `matplotlib`, `h5py`, Jupyter, and a
LaTeX distribution with `siunitx` (the figures are typeset with LaTeX).

## Running

Start Jupyter from this folder: the paths are relative to it. The notebook reads

- `npz/`, the bench captures,
- `../matlab/results/runs/`, the simulated runs it compares against,

and writes the figures to `figures/rapport/`.

It runs top to bottom: data first, then experiments E1 to E3 and the inner
loop, then tasks T1 to T4.
