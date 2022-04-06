#############################################################################
##
##  HAPPRIME - polynomials.gd
##  Functions, Operations and Methods to extend GAP's polynomials
##  Paul Smith
##
##  Copyright (C)  2008
##  Paul Smith
##  National University of Ireland Galway
##
##  This file is part of HAPprime. 
## 
##  HAPprime is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
## 
##  HAPprime is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
## 
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <https://www.gnu.org/licenses/>.
## 
##
#############################################################################


#######################################################################
#
# Declarations for Operations in polynomials.gi
#
DeclareAttribute("TermsOfPolynomial", IsPolynomial);
DeclareAttribute("UnivariateMonomialsOfMonomial", IsPolynomial);
DeclareAttribute("IndeterminateAndExponentOfUnivariateMonomial", IsPolynomial);
DeclareAttribute("IndeterminatesOfPolynomial", IsPolynomial);

DeclareOperation("ReduceIdeal",[IsHomogeneousList and IsRationalFunctionCollection, IsMonomialOrdering]);
DeclareOperation("ReducedPolynomialRingPresentation", [IsPolynomialRing, IsHomogeneousList and IsRationalFunctionCollection, IsHomogeneousList]);
DeclareOperation("ReducedPolynomialRingPresentationMap", [IsPolynomialRing, IsHomogeneousList and IsRationalFunctionCollection, IsHomogeneousList]);

DeclareGlobalFunction("HAPPRIME_SwitchPolynomialIndeterminates");
DeclareGlobalFunction("HAPPRIME_MapPolynomialIndeterminates");
DeclareGlobalFunction("HAPPRIME_CombineIndeterminateMaps");
DeclareGlobalFunction("HAPPRIME_SingularGroebnerBasis");
DeclareGlobalFunction("HAPPRIME_SingularReducedGroebnerBasis");
