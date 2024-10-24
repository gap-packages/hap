[![Build Status](https://github.com/gap-packages/hap/workflows/CI/badge.svg?branch=master)](https://github.com/gap-packages/hap/actions?query=workflow%3ACI+branch%3Amaster)
[![Code Coverage](https://codecov.io/github/gap-packages/hap/coverage.svg?branch=master&token=)](https://codecov.io/gh/gap-packages/hap)

# The HAP package for group cohomology and related functors

(C) Graham Ellis

HAP is a package for some calculations in elementary algebraic topology and
the cohomology of groups. The initial focus of the library was on computations 
related to the cohomology of finite and infinite groups, with particular 
emphasis on integral coefficients. The focus has since broadened to include 
Steenrod algebras of finite groups, Bredon homology, cohomology of simplicial 
groups, and general computations in algebraic topology relating to finite 
CW-complexes, covering spaces, knots, knotted surfaces, and topics such as
 persitent homology arising in topological data analysis. A full description 
of the library can be found at https://gap-packages.github.io/hap/ .

HAP is written almost entirely in the GAP language. A few of HAP's (optional)
functions rely on Polymake software (available from
www.math.tu-berlin.de/polymake/) which should be installed so that it starts
using the command "polymake" at the Linux prompt. A few (optional) 
functions rely on Singular software (available from 
https://www.singular.uni-kl.de/) which should be installed so that it starts 
using the command "singular" at the Linux prompt. A few (oprional) functions
rely on Graphviz software (available from http://www.graphviz.org/).

Please send your bug reports to graham.ellis(at)nuigalway.ie .


## Installing HAP

On a Linux machine with GAP (and optionally Polymake) installed, the HAP
library can be loaded as follows:

* First download the file hap1.66.tar.gz to the subdirectory "pkg/" of GAP. (If
you don't have access to this, then create a directory "pkg" in your home
directory and download the file there.)

* Change to directory "pkg/" and type "gunzip hap1.66.tar.gz" followed by
"tar -xvf hap1.66.tar" .

* Start GAP. (If you have created "pkg" in your home directory then start GAP
with the command "gap -l 'path/homedir;' "   where path/homedir is the path to
your home directory.)

* In GAP type " LoadPackage("HAP"); " .

* Help on HAP can be found on the HAP home page (a version of which is
included in directory "pkg/Hap1.66/www" of this distribution).

* Performance can be significantly improved by using a compiled version of the
HAP library. A compiled version can be created by the following steps.

1. Change to the directory "pkg/Hap1.66/" .
2. Edit the file "compile" so that: PKGDIR is equal to the path to the
directory "pkg" where your GAP packages are stored; GACDIR is equal to the
path to the directory where the GAP compiler "gac" is stored.
3. Type "./compile" .

The next time HAP is loaded a compiled version will be loaded.

* Should you want to return to an uncompiled version, change to the directory
"pkg/Hap1.66/" and type "./uncompile".
