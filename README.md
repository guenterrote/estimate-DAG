# estimate-DAG
Estimating the number of nodes and paths in a DAG by following a random path

This repository contains the programs and experimental results for the paper
"Estimating the Number of Nodes and Paths in a DAG by Following a Random Path"
by Günter Rote.
The full version of the paper is in the file fullpaper.pdf.

The python programs
`estimating-DAG-permutations.py` and `estimating_DAG_klee_minty.py`
perform experiments for two
types of examples: the permutation DAG and Klee-Minty cubes.
The size parameter n can be given on the command line, for example
```
    pypy3 estimating_DAG_klee_minty.py 10
```
(`estimating-DAG-permutations.py` imports some functions from
`estimating_DAG_klee_minty.py`, so both files must be in the same directory.)

The following methods are tested:
- Algorithm P (Path Sampling), 10^7 samples;
- Importance Sampling (IS) by outdegree, 10^7 samples;
- for Klee-Minty cubes only: Importance Sampling by the size of the
  reachable set (IS*) and by size/indegree (IS**), 10^7 samples;
- Bucket Sampling (BS) by outdegree, 10^6 samples;
- Crowd Sampling (CS) with budget B=10, 10^6 samples.

For Algorithm P and the variants of Importance Sampling,
the programs also compute the exact standard deviation of the estimate.
For computing the exact variance, the programs need a lot of memory.
It is preferable to run these programs with the pypy interpreter
instead of the standard python interpreter, which may
generate a "MemoryError" on larger instances.

The results are logged to the standard output (see the files
`L-...log`; the last part of the filename refers to the computer on
which the program was run),
and they are written to files of the form `results-PERM-n.py` and `results-KM-n.py`,
for further processing.
The python script `analyze.py` reads these result files and creates
tables in LaTeX format
for the paper
in the subdirectory `tables`.
The file
`estimating-DAG-full-tables.pdf`
(generated from 
`estimating-DAG-full-tables.tex`) in that directory
contains the full results for the permutation DAG for n=7,8,...,12,
and for Klee-Minty cubes of dimension n=8,10,12,...,18.

The method can also estimate the number of source-sink paths, but no
experiments have been carried out in this direction.
