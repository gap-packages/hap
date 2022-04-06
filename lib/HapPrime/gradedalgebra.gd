#############################################################################
##
##  HAPPRIME - gradedalgebra.gd
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
##  <#GAPDoc Label="GradedAlgebraPresentationFilter_DTmanGradedAlgebraNODOC">
##  <ManSection>
##  <Filt Name="IsGradedAlgebraPresentation" Arg="O" Type="Category"/>
##  <Returns>
##  <K>true</K> if the object is a <K>GradedAlgebraPresentation</K>, or 
##  <K>false</K> otherwise
##  </Returns>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareCategory("IsGradedAlgebraPresentation", IsObject);
# Note this also defines the function IsGradedAlgebraPresentation
#####################################################################

#####################################################################
##  <#GAPDoc Label="GradedAlgebraPresentationFamilyAttr_DTmanGradedAlgebraNODOC">
##  <ManSection>
##  <Fam Name="GradedAlgebraPresentationFamily"/>
##  <Description>
##  The family to which <K>GradedAlgebraPresentation</K> objects belong.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareAttribute("GradedAlgebraPresentationFamily", IsFamily);
# Note this also defines the function GradedAlgebraPresentationFamily
#####################################################################


#######################################################################
#
# Declarations for Operations in derivation.gi
#
DeclareOperation("GradedAlgebraPresentation", [IsPolynomialRing, IsHomogeneousList, IsHomogeneousList]);
DeclareOperation("GradedAlgebraPresentationNC", [IsPolynomialRing, IsHomogeneousList, IsHomogeneousList]);

DeclareAttribute("BaseRing", IsGradedAlgebraPresentation);
DeclareAttribute("CoefficientsRing", IsGradedAlgebraPresentation);
DeclareAttribute("IndeterminatesOfGradedAlgebraPresentation", IsGradedAlgebraPresentation);
DeclareAttribute("PresentationIdeal", IsGradedAlgebraPresentation);
DeclareAttribute("GeneratorsOfPresentationIdeal", IsGradedAlgebraPresentation);
DeclareAttribute("IndeterminateDegrees", IsGradedAlgebraPresentation);

DeclareOperation("TensorProductOp", [IsGradedAlgebraPresentation, IsGradedAlgebraPresentation]);

DeclareOperation("DegreeOfRepresentative", [IsGradedAlgebraPresentation, IsPolynomial]);
DeclareOperation("AreIsomorphicGradedAlgebras", [IsGradedAlgebraPresentation, IsGradedAlgebraPresentation]);
DeclareOperation("IsAssociatedGradedRing", [IsGradedAlgebraPresentation, IsGradedAlgebraPresentation]);
DeclareAttribute("MaximumDegreeForPresentation", IsGradedAlgebraPresentation);
DeclareOperation("SubspaceDimensionDegree", [IsGradedAlgebraPresentation, IsPosInt]);
DeclareOperation("SubspaceBasisRepsByDegree", [IsGradedAlgebraPresentation, IsPosInt]);
    
DeclareOperation("CoefficientsOfPoincareSeries", [IsGradedAlgebraPresentation, IsPosInt]);
DeclareAttribute("HilbertPoincareSeries", IsGradedAlgebraPresentation);
  
DeclareAttribute("HAPPRIME_HilbertSeries", IsGradedAlgebraPresentation);

DeclareGlobalFunction("HAPPRIME_SwitchGradedAlgebraRing");
DeclareGlobalFunction("HAPPRIME_GradedAlgebraPresentationAvoidingIndeterminates");
