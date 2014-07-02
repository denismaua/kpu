k-Neighborhood Policy Updating for Influence Diagrams
=====================================================

A Python package that implements the *k-neighborhood local search algorithms for selecting strategies in (limited memory) influence diagrams* described in

>    D.D. Maua and F.G. Cozman (2012), Speeding Up k-Neighborhood Local Search in Limited Memory Influence Diagrams. In Proceedings of the Seventh European Conference on Probabilistic Graphical Models (to appear).
    
If you use this code, please cite the publication above.

LICENSE
-------
    
Copyright by Denis D. Maua (2014)

THIS CODE IS PROVIDED "AS-IS". USE AT YOUR OWN RISK.

This package is released so that others can reproduce the experiments in the paper, and potentially use the algorithms in further work. You can contact-me for simple questions, but please do not expect to get real support from me.

INSTALLATION
------------

Simply clone the repository. It is highly advisable to run the scripts with the Pypy interpreter instead of the standard Python intepreter.
Check [pypy.org](http://pypy.org "Pypy") for info on downloading and installing PyPy.

USAGE
-----

*** If you have Pypy installed, replace `python` with `pypy` in the sequel.

Every runnable script returns its usage if run without arguments (unless the function implemented takes no arguments). Try `python spu3.py` to get usage on running SPU with dominance pruning.

The diagrams used in the experiments in the paper are in the compressed file experiments.zip. To reproduce the results, uncompress that file and the script run-chain.py with each algorithm, e.g.:
   
>   `python run-chain.py SPUhull experiments/exp-100-10-?.p`
   
will run SPU (1-neighborhood search) with dominance pruning on ten chain-like diagrams with 100 decision variables and 10-ary variables (each file is a pick dump of a diagram object).

New chain-like diagrams can be generated using the script `generate-chain.py`. Comparison among different algorithms on randomly generated chain-like diagrams can be obtained running the script `chainID.py`
