# NK

[![Build Status](https://travis-ci.org/glesica/NK.jl.svg?branch=master)](https://travis-ci.org/glesica/NK.jl)

A Julia library for conducting evolutionary experiments using the NK family of
fitness landscapes.

## Tests

Tests can be run with the included script. There are (or will be) several test
suites. First, a set of unit tests that exercise basic functionality and check
for sane results. For the most part, these may be thought of as regression
tests. Second, there are functional tests that replicate experiments drawn
from the literature and check that the results approximately match the
published results.

The unit tests may be run with `./runtests.sh unit`. The functional tests may
be run with `./runtests.sh SUITE` where `SUITE` is one of the following:

  * `kauffman` - experiments from [1]
  * `nowak` - experiments from [2]

Note that these tests may take quite a long time to finish even on very
powerful workstations.

## Bibliography

  1. Kauffman, S. A. (1993). The Origins of Order: Self-Organization and
     Selection in Evolution. New York, New York, USA: Oxford University Press.
  2. Nowak, S., & Krug, J. (2015). Analysis of adaptive walks on NK fitness
     landscapes with different interaction schemes. Journal of Statistical
     Mechanics: Theory and Experiment, 2015(6), P06014.
     doi:10.1088/1742-5468/2015/06/P06014

## TODO

  * Look at implementing something similar to Newman and Englehart.
  * Make a list of mutation algorithms to implement.
  * Check to see if neutral degree drops off exp. with increasing fitness.

