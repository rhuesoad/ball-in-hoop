function run = load_run(filepath)
%LOAD_RUN  Loads a simulation run struct from a .mat file.
%
%   run = load_run(filepath)
%
%   Symmetric counterpart to save_run: loads the run struct that save_run
%   wrote. The variable name inside the .mat file is always 'run'.
%
%   Input
%   -----
%   filepath : full or relative path to the .mat file
%
%   Output
%   ------
%   run : run struct obeying results/RUN_SCHEMA.md
%         (see analysis/io/build_run.m, which assembles it)

data = load(filepath, 'run');
run  = data.run;
end
