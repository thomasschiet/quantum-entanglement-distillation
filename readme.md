# Distillation project
The implementations of the methods described in the notes are found in the `/lib` folder. For example the code for purity enforcing SDPs is found in `/lib/purity.jl`. The programmes are then run in the root folder.

## PPT relaxations

## Epsilon net
The SDPs and methods to build the epsilon net are in `/lib/constructEpsnet.jl`. It can be run using `epsilonNet.jl`, but some things still need to be implemented:
- Combine the found Choi states for Alice and Bob and compute output fidelity
- Parallelize code to run it on Amazon

## k-extentions
To run a single extension the function `PPTprogrammeNoTwirling1Ext` is available and can be seen in `lib/k-ext.jl`. It takes about 5-10 minutes to run. 

## Purity enforcing SDPs

## Seesaw methods
The code for the seesaw method is currently written in MATLAB, but will be ported to Julia soon.

## TODO:
- Finish purity constrains and loop and make a script to run it
- Finish epsilon net
- Clean up some more old code
