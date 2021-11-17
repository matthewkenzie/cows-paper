# Custom Orthogonal Weight functions (COWs) - Accompanying Material

The material in this repository is intended to supplement the paper: ["Custom Orthogonal Weight functions (COWs) for Event Classification", Dembinski, Kenzie, Langenbruch and Schmelling](arXiV.org)

Provided here are a few helper modules / classes which can be used to extract sWeights and COWs:

- `examples.py` demonstrates how each of the classes provided works and will:
  1. Generate some toy data
  
  ![toy](https://user-images.githubusercontent.com/1140576/142237277-0485e6e7-8ccf-489a-affd-6b81028ed5c3.png)

  2. Fit the toy data in the discriminanting variable to get an estimate of the discriminating variable pdfs

  3. Run the "summation" sWeights method (using the `SWeight` class provided below) ![sws](https://user-images.githubusercontent.com/1140576/142237391-0b37f428-5668-4602-98bb-097fdaae62e8.png)
  4. Run the COW method with variance function of unity, I(m)=1, (using the `Cow` class provided below) ![cows](https://user-images.githubusercontent.com/1140576/142237453-8c3dfa2b-b38d-4e22-96d8-30f31f61d1c8.png) 

  5. Fit the weighted distributions and correct the covariance using the `cov_correct` function provided below ![tfit](https://user-images.githubusercontent.com/1140576/142237505-11032b1c-b6fa-47dc-9a0e-e965210fdf6b.png)

- `SWeighter.py` - this provides a class which implements sWeights in 6 different ways (depending on users desire)
  1. Using the "summation" method from the original sPlot paper (referred to as Variant B in our paper)
  2. Using the "integration" method rederived in our paper, originally by Schmelling, (referred to as Variant A in our paper)
  3. Using the "refit" method, i.e. taking the covariance matrix of a yield only fit (referred to as Variant Ci in our paper)
  4. Using the "subhess" method, i.e. taking the sub-covariance matrix for the yields (referred to as Variant Cii in our paper)
  5. Using the implementation in ROOT's TSPlot (this we believe should be equivalent to Variant B but is more susceptible to numerical differences)
  6. Using the implementation in RooStat's SPlot (we found this identical to Variant B ("summation") above in all the cases we tried)

- `CovarianceCorrect.py` - this provides a function which implements a correction to the covariance matrix based on the fit result. This follows the original work done by Langenbruch in [arXiv:1911.01303](https://arxiv.org/abs/1911.01303).

- `Cow.py` - this provides a class which implements COWs in a variety of ways (depending on users desire). It expects to be passed:
  - gs - the signal function for the discrimant variable
  - gb - the background function(s) for the discriminant variable (can pass orders of polynomials here if desired)
  - Im - the variance function for the discriminant variance (can also pass 1 and it will be set of uniform)
  - obs - one can instead or additionally pass the observed distribution in the discriminant variable which will be used for the variance function instead. In this case you must pass a two element tuple giving the bin entries and bin edges for the observed dataset (the same as what `np.histogram(data)` would return)


