## In this file I want to draw attention to the fact, that the code for B-NEM evolved from the CellNetOptimizer method from Saez-Rodriguez et al. 2009 or to be more precise from the R-package CellNOptR (https://github.com/cellnopt/CellNOptR).

## That is also why B-NEM requires the CNO R-package to run and still uses objects and some routines from CellNOptR. E.g. the CNOlist S4-class. But several functions have been rewritten (simulateStatesRecursive function to simulate the steady state of a Boolean network or disjunctive normal form is completely new).

## References:
  
## https://github.com/cellnopt/CellNOptR

## Saez-Rodriguez, Julio, Alexopoulos, Leonidas G, Epperlein, Jonathan, Samaga, Regina, Lauffenburger, Douglas A, Klamt, Steffen,
## & Sorger, Peter K. 2009. Discrete logic modelling as a means to link protein signalling networks with functional analysis of
## mammalian signal transduction. Mol Syst Biol, 5, 331.
