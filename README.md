# Accompanying material for the Custom Orthogonal Weight functions (COWs) paper

The material in this repository is intended to supplement the paper: ["Custom Orthogonal Weight functions (COWs) for Event Classification", Dembinski, Kenzie, Langenbruch and Schmelling](arXiV.org)

Provided here are a few helper modules / classes which can be used to extract sWeights and COWs:

- `examples.py` demonstrates how each of the classes provided works

- `SWeighter.py` - this provides a class which implements sWeights in 6 different ways (depending on users desire)
  1. Using the "summation" method from the original sPlot paper (referred to as Variant B in our paper)
  2. Using the "integration" method rederived in our paper, originally by Schmelling, (referred to as Variant A in our paper)
  3. Using the "refit" method, i.e. taking the covariance matrix of a yield only fit (referred to as Variant Ci in our paper)
  4. Using the "subhess" method, i.e. taking the sub-covariance matrix for the yields (referred to as Variant Cii in our paper)
  5. Using the implementation in ROOT's TSPlot (this we believe should be equivalent to Variant B but is more susceptible to numerical differences)
  6. Using the implementation in RooStat's SPlot (we found this identical to Variant B ("summation") above in all the cases we tried)

- `CovarianceCorrect.py` - this provides a function which implements a correction to the covariance matrix based on the fit result. This follows the original work done by Langenbruch in arXiv:1911.01303.

- `Cow.py` - this provides a class which implements COWs in a variety of ways (depending on users desire). It expects to be passed:
  - gs - the signal function for the discrimant variable
	- gb - the background function(s) for the discriminant variable (can pass orders of polynomials here if desired)
	- Im - the variance function for the discriminant variance (can also pass 1 and it will be set of uniform)
	- obs - one can instead or additionally pass the observed distribution in the discriminant variable which will be used for the variance function instead. In this case you must pass a two element tuple giving the bin entries and bin edges for the observed dataset (the same as what `np.histogram(data)` would return)


