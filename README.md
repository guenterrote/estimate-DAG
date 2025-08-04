# estimate-DAG
Estimating the size of a DAG by following random paths

The python programs
estimating-DAG...py perform experiments for two
types of examples: the permutation DAG and Klee-Minty cubes.
The size parameter n can be given on the command line.
It is preferable to process these files with the pypy interpreter
instead of the standard python interpreter, which may
have a "MemoryError" on larger instances.

The results are written to the standard output, as well as to files of the form result...py,
for further processing.
The python script analyze.py reads these result files and creates the
tables in LaTeX format
for the paper "Estimating the size of a DAG by following random paths"
in the subdirectory tables.
The file
estimating-DAG-full-tables.pdf
(generated from 
estimating-DAG-full-tables.tex) in that directory
contains the full results for the permutation DAG for n=7,8,...,12,
and for Klee-Minty cubes of dimension n=8,10,12,...,18.