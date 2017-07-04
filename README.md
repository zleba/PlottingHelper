# PlottingHelper
Plotting utilities to makes plotting in ROOT a bit easier

Includes several methods for positioning of titles, labels and defining their sizes.
Furthermore, package contains functions for creation of grid of pads.
There is also an testing version of the automatic legend.

### Compilation

To use this package in script, you first need to complile it by calling
```
$ make
```

This shell command creates two shared libraries, one called `libPlottingHelper.so` using `g++` compiler.
Another called `plottingHelper_C.so` using compiler in ROOT.
For the time being the `g++` version was found much faster and therefore is recommended for both Macros, compiled Macros and ROOT programs compiled by `g++`.

### Using in Macros

To use it put these lines to the beginning of your program:
```
#include "plottingHelper.h"
R__LOAD_LIBRARY(libPlottingHelper.so) //necessary only for ROOT Macros
using namespace PlottingHelper;//pollute the namespace!
```

In case of some problems one can alternatively load the library created by ROOT:
```
#include "plottingHelper.h"
R__LOAD_LIBRARY(plottingHelper_C.so)
using namespace PlottingHelper;
```
### Using in programs compiled by g++ with ROOT

For programs complied by g++, the library must be added to the Makefile by the standard way:
```
-Wl,-rpath,. -L. -lPlottingHelper
```
where both dots "." can be replaced by the acutall location of the dynamic library.
The `R__LOAD_LIBRARY` does nothing for compiled programs and therefore can be kept.


### Documentation

The code is documented by doxygen, in future the tutorial section is planned as well.

https://codedocs.xyz/zleba/PlottingHelper/namespacePlottingHelper.html

[![Documentation](https://codedocs.xyz/zleba/PlottingHelper.svg)](https://codedocs.xyz/zleba/PlottingHelper/)
