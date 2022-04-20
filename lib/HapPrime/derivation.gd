#############################################################################
##
##  HAPPRIME - derivation.gd
##  Functions, Operations and Methods to implement derivations
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



#####################################################################
##  <#GAPDoc Label="HAPDerivationFilter_DTmanDerivationNODOC">
##  <ManSection>
##  <Filt Name="IsHAPDerivation" Arg="O" Type="Category"/>
##  <Returns>
##  <K>true</K> if the object is a <K>HAPDerivation</K>, or 
##  <K>false</K> otherwise
##  </Returns>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareCategory("IsHAPDerivation", IsObject);
# Note this also defines the function IsHAPDerivation
#####################################################################

#####################################################################
##  <#GAPDoc Label="HAPDerivationFamilyAttr_DTmanDerivationNODOC">
##  <ManSection>
##  <Fam Name="HAPDerivationFamily"/>
##  <Description>
##  The family to which <K>HAPDerivation</K> objects belong.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareAttribute("HAPDerivationFamily", IsFamily);
# Note this also defines the function HAPDerivationFamily
#####################################################################


#######################################################################
#
# Declarations for Operations in derivation.gi
#
DeclareOperation("HAPDerivation", [IsPolynomialRing, IsHomogeneousList, IsHomogeneousList and IsRationalFunctionCollection]);
DeclareOperation("HAPDerivationNC", [IsPolynomialRing, IsHomogeneousList, IsHomogeneousList and IsRationalFunctionCollection]);

DeclareAttribute("DerivationRing", IsHAPDerivation);
DeclareAttribute("DerivationRelations", IsHAPDerivation);
DeclareAttribute("DerivationImages", IsHAPDerivation);

DeclareOperation("PolynomialToRModuleRep", [IsPolynomialRing, IsPolynomial]);
DeclareOperation("ImageOfDerivation", [IsHAPDerivation, IsPolynomial]);
DeclareOperation("KernelOfDerivation", [IsHAPDerivation, IsHomogeneousList]);
DeclareOperation("HomologyOfDerivation", [IsHAPDerivation, IsHomogeneousList]);

DeclareGlobalFunction("HAPPRIME_SModule");