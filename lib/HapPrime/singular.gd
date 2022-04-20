#############################################################################
##
##  HAPPRIME - singular.gd
##  Functions, Operations and Methods to interface with singular
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

DeclareAttribute("SingularGroebnerBasis", IsPolynomialRingIdeal); 
DeclareAttribute("SingularReducedGroebnerBasis", IsPolynomialRingIdeal); 

DeclareOperation("SingularSetNormalFormIdeal", [IsPolynomialRingIdeal]); 
DeclareOperation("SingularSetNormalFormIdealNC", [IsPolynomialRingIdeal]); 
DeclareOperation("SingularPolynomialNormalForm", 
  [IsPolynomial, IsPolynomialRingIdeal]); 



