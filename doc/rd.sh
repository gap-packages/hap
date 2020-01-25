#!/bin/bash
read version < ./version;
read date < ./date;

# A script to update the version number in various pages.

file=$(more<<STOP
[![Build Status](https://travis-ci.org/gap-packages/hap.svg?branch=master)](https://travis-ci.org/gap-packages/hap)
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
of the library can be found at http://hamilton.nuigalway.ie/Hap/www/ .

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

* First download the file hap1.25.tar.gz to the subdirectory "pkg/" of GAP. (If
you don't have access to this, then create a directory "pkg" in your home
directory and download the file there.)

* Change to directory "pkg/" and type "gunzip hap$version.tar.gz" followed by
"tar -xvf hap$version.tar" .

* Start GAP. (If you have created "pkg" in your home directory then start GAP
with the command "gap -l 'path/homedir;' "   where path/homedir is the path to
your home directory.)

* In GAP type " LoadPackage("HAP"); " .

* Help on HAP can be found on the HAP home page (a version of which is
included in directory "pkg/Hap$version/www" of this distribution).

* Performance can be significantly improved by using a compiled version of the
HAP library. A compiled version can be created by the following steps.

1. Change to the directory "pkg/Hap1.25/" .
2. Edit the file "compile" so that: PKGDIR is equal to the path to the
directory "pkg" where your GAP packages are stored; GACDIR is equal to the
path to the directory where the GAP compiler "gac" is stored.
3. Type "./compile" .

The next time HAP is loaded a compiled version will be loaded.

* Should you want to return to an uncompiled version, change to the directory
"pkg/Hap$version/" and type "./uncompile".
STOP
);

echo "$file" >./README.md;

file=$(more<<STOP
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<!-- --------------------HEADER------------------------ -->
  <title>HAP - content</title>
  <meta http-equiv="Content-Type"
 content="text/html; charset=iso-8859-1">
  <link rel="Stylesheet" href="../HAPstyles.css" type="text/css"
 media="screen,print">
  <style type="text/css"></style>
</head>
<body id="content" style="background-image: url(../background.png);">
<div id="content">
<p>HAP is a homological
algebra library for use with the
GAP computer algebra system, and is still under development.
 The current version <a href="../../hap$version.tar.gz">$version</a> was released on
$date.</p>


<p>
The initial focus of the library was on computations related to the cohomology of finite and infinite groups, with particular emphasis on integral coefficients. The focus has since broadened to include Steenrod algebras of finite groups, Bredon homology, cohomology of simplicial groups, and general computations in algebraic topology relating to finite CW-complexes, covering spaces, knots, knotted surfaces, and topics such as persitent homology arising in topological data analysis. 


 This tutorial <a href="../SideLinks/About/aboutContents.html" target="_top">table of contents</a> gives an indication of  topics covered. This <a href="https://global.oup.com/academic/product/an-invitation-to-computational-homotopy-9780198832980" target="_top">book</a> provides a detailed overview.</p>

<p> A list of some papers which refer to HAP can be found <a href="hap.pdf">here</a>.



<p>Some functions rely on the deposited packages <a href=../../hapcryst.tar.gz>HAPcryst</a> and <a href="../../polymaking.tar.gz">Polymaking</a> and on the package <a href="../../hapcocyclic.tar.gz">HAPcocyclic</a>. </p>

<p> <a href="../copyright/copyrightIndex.html" target="_top"> <small>GPL
Copyright</small></a> </p>

</div>
</body>
</html>
STOP
);

echo "$file" > ./www/home/content.html;

file=$(more<<STOP
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<!-- --------------------HEADER------------------------ -->
  <title>HAP - content download</title>
  <meta http-equiv="Content-Type"
 content="text/html; charset=iso-8859-1">
  <link rel="Stylesheet" href="../HAPstyles.css" type="text/css"
 media="screen,print">
  <style type="text/css"></style>
</head>
<body id="content" style="background-image: url(../background.png);">
<div id="content">
<p></p>
<h3>Download Instructions</h3>
<p>
On a Linux machine with GAP installed, the HAP library can be loaded as
follows:<br>
</p>
<ul>
  <li>First download the file <a
 href="http://hamilton.nuigalway.ie/Hap/hap$version.tar.gz">hap$version.tar.gz</a>
which contains the most recent development version of HAP to the
subdirectory
"pkg/" of GAP. If you don't have access to this subdirectory,
then create a directory "pkg" in your home directory and download the
file there. (if for some reason you'd prefer to download the most
recent non-development version of HAP then download the file&nbsp; <a
 href="http://hamilton.nuigalway.ie/Hap/hap$version.tar.gz">hap$version.tar.gz</a>
instead.)<br>
  </li>
</ul>
<ul>
  <li>Change to directory "pkg/" and type "gunzip
hap$version.tar.gz"
followed
by "tar
-xvf hap$version.tar" .</li>
</ul>
<ul>
  <li>Start GAP. (If you have created "pkg" in your home
directory then
start GAP with the command "gap -l 'path/homedir;' " &nbsp; where
path/homedir is the path to your home directory.)</li>
</ul>
<ul>
  <li>In GAP type " LoadPackage("HAP"); "<br>
  </li>
</ul>
<ul>
  <li>Help on HAP can be found on the <a
 href="http://hamilton.nuigalway.ie/Hap/www/index.html" target="index">HAP
home page</a> (a version of which is included in directory
"pkg/Hap$version/www" of the distribution).</li>
</ul>
<ul>
  <li>A few of HAP's (optional) functions rely on Polymake
software
(available from www.math.tu-berlin.de/polymake/) which should be
installed so that it
starts using the command "polymake" at the Linux prompt. A few
(optional) functions rely on Graphviz software (available from
http://www.graphviz.org/).<br>
  </li>
</ul>
<ul>
  <li>Performance can be improved by using a compiled
version of the HAP library. The following steps will produce a compiled
version. <br>
(1) Change to the directory "pkg/Hap$version/" .<br>
(2) Edit the file "compile" so that: PKGDIR is equal to the path to the<br>
directory "pkg" where your GAP packages are stored; GACDIR is equal to
the<br>
path to the directory where the GAP compiler "gac" is stored.<br>
(3) Type "./compile" .<br>
The next time HAP is loaded a compiled version will be loaded.</li>
</ul>
<ul>
  <li>Should you want to return to an uncompiled version, change
to the directory<br>
"pkg/Hap$version/" and type "./uncompile".</li>
</ul>
</div>
</body>
</html>
STOP
);

echo "$file" > ./www/download/downloadContent.html;

      
