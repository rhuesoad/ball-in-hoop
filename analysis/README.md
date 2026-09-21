# Analysis — report figures

`bench_analysis.ipynb` turns the bench captures into the figures of the thesis.

## Requirements

Python 3 with `numpy`, `scipy`, `pandas`, `matplotlib`, `h5py`, Jupyter, and a
LaTeX distribution with `siunitx` (the figures are typeset with LaTeX).

## Running

Start Jupyter from this folder: all paths are relative to it. The notebook will read:

- `npz/`, the bench captures, in a "npz" format, 
- `../matlab/results/runs/`, the simulated runs it compares against,

and writes the figures to `figures/rapport/`.

It runs top to bottom: 
the data is first extracted, analyzed and sorted, the experiments are then used to extract the identification parameters, the inner loop tuning results are then analysed, and finally the four tasks T1 to t4. 
