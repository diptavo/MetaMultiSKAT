### About

MetaMultiSKAT provides a meta-analysis framework for performing rare-variant based tests of associations with multiple phenotypes across multiple studies. It is a kernel regression framework extension of [MultiSKAT](https://github.com/diptavo/MultiSKAT). This is an initial working version. A number of functionalities will be added soon.


### MetaMultiSKAT v(1.0)

This package contains the R-codes/functions (including an example dataset) to carry out the MetaMultiSKAT tests. The package is still under developement. Any questions or bug-reports should be mailed to diptavo@umich.edu

### Main Functions

**Meta.MultiSKAT.base**: Performs MetaMultiSKAT tests with a given Sigma_s, Sigma_P and Sigma_G structure

**Meta.MultiSKAT.wResample**: Performs MetaMultiSKAT tests with a given Sigma_s, Sigma_P and Sigma_G structure and generates null p-values to be used in omnibus tests

**minP.Meta**: Minimum p-value based aggregate tests using copula

**Meta.Hom**: Performs MetaMultiSKAT test of association against the alternative that all the studies have similar effects

**Meta.Het**: Performs MetaMultiSKAT test of association against the alternative that all the studies have heterogeneous/uncorrelated effects

**Meta.Com**: Combination test of Meta.Het and Meta.Hom
