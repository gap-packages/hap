#############################################################################
##
##  HAPPRIME - ringhomomorphism.gd
##  Functions, Operations and Methods to implement ring homomorphisms
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
##  <#GAPDoc Label="HAPRingHomomorphismFilter_DTmanRingHomomorphismNODOC">
##  <ManSection>
##  <Filt Name="IsHAPRingHomomorphism" Arg="O" Type="Category"/>
##  <Returns>
##  <K>true</K> if the object is a <K>HAPRingHomomorphism</K>, or 
##  <K>false</K> otherwise
##  </Returns>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareCategory("IsHAPRingHomomorphism", IsObject);
# Note this also defines the function IsHAPRingHomomorphism
#####################################################################

#####################################################################
##  <#GAPDoc Label="HAPRingHomomorphismFamilyAttr_DTmanRingHomomorphismNODOC">
##  <ManSection>
##  <Fam Name="HAPRingHomomorphismFamily"/>
##  <Description>
##  The family to which <K>HAPRingHomomorphism</K> objects belong.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareAttribute("HAPRingHomomorphismFamily", IsFamily);
# Note this also defines the function HAPRingHomomorphismFamily
#####################################################################


#######################################################################
#
# Declarations for Operations in ringhomomorphism.gi
#
DeclareOperation("HAPRingToSubringHomomorphism", 
  [IsPolynomialRing, IsHomogeneousList, 
  IsHomogeneousList and IsRationalFunctionCollection]);
DeclareOperation("HAPSubringToRingHomomorphism", 
  [IsHomogeneousList and IsRationalFunctionCollection, 
  IsHomogeneousList, IsPolynomialRing]);
DeclareOperation("HAPSubringToRingHomomorphism", 
  [IsHomogeneousList and IsRationalFunctionCollection, 
  IsPolynomialRing, IsHomogeneousList]);
DeclareOperation("HAPRingHomomorphismByIndeterminateMap", 
  [IsPolynomialRing, IsHomogeneousList, IsPolynomialRing]);
DeclareOperation("HAPRingReductionHomomorphism", 
  [IsPolynomialRing, IsHomogeneousList, IsHomogeneousList]);
DeclareOperation("HAPRingReductionHomomorphism", 
  [IsHAPRingHomomorphism, IsHomogeneousList]);
DeclareOperation("HAPZeroRingHomomorphism", 
  [IsPolynomialRing, IsHomogeneousList]);
DeclareOperation("CompositionRingHomomorphism", 
  [IsHAPRingHomomorphism, IsHAPRingHomomorphism]);
DeclareAttribute("InverseRingHomomorphism", IsHAPRingHomomorphism);

DeclareAttribute("SourcePolynomialRing", IsHAPRingHomomorphism);
DeclareAttribute("SourceGenerators", IsHAPRingHomomorphism);
DeclareAttribute("SourceRelations", IsHAPRingHomomorphism);
DeclareAttribute("ImagePolynomialRing", IsHAPRingHomomorphism);
DeclareAttribute("ImageGenerators", IsHAPRingHomomorphism);
DeclareAttribute("ImageRelations", IsHAPRingHomomorphism);

DeclareOperation("ImageOfRingHomomorphism", 
  [IsHAPRingHomomorphism, IsHomogeneousList and IsRationalFunctionCollection]);
DeclareOperation("PreimageOfRingHomomorphism", 
  [IsHAPRingHomomorphism, IsHomogeneousList and IsRationalFunctionCollection]);

DeclareGlobalFunction("HAPPRIME_MakeEliminationOrdering");
DeclareGlobalFunction("HAPPRIME_RingHomomorphismsAreComposable");
  
#####################################################################
##  <#GAPDoc Label="SourceGenerators_DTmanRingHomomorphism_Dat">
##  <ManSection>
##  <Attr Name="SourceGenerators" Arg="phi"/>
##  <Returns>
##    List
##  </Returns>
##  <Description>
##    A list of generators for the source ring <M>R/I</M> of the ring 
##    homomorphism.
##    <A>phi</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
##  <#GAPDoc Label="SourceRelations_DTmanRingHomomorphism_Dat">
##  <ManSection>
##  <Attr Name="SourceRelations" Arg="phi"/>
##  <Returns>
##    List
##  </Returns>
##  <Description>
##    A list of the relations that generate the ideal <M>I</M> of in the source 
##    ring of the ring homomorphism <A>phi</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
##  <#GAPDoc Label="SourcePolynomialRing_DTmanRingHomomorphism_Dat">
##  <ManSection>
##  <Attr Name="SourcePolynomialRing" Arg="phi"/>
##  <Returns>
##    <K>PolynomialRing</K>
##  </Returns>
##  <Description>
##    Returns the polynomial ring which contains the source ring of the 
##    ring homomorphism <A>phi</A>.
##    Polynomials to be mapped by <A>phi</A> must be in this ring.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
##  <#GAPDoc Label="ImageGenerators_DTmanRingHomomorphism_Dat">
##  <ManSection>
##  <Attr Name="ImageGenerators" Arg="phi"/>
##  <Returns>
##    List
##  </Returns>
##  <Description>
##    A list of generators for the image ring <M>S/J</M> of the ring 
##    homomorphism <A>phi</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
##  <#GAPDoc Label="ImageRelations_DTmanRingHomomorphism_Dat">
##  <ManSection>
##  <Attr Name="ImageRelations" Arg="phi"/>
##  <Returns>
##    List
##  </Returns>
##  <Description>
##    A list of the relations that generate the ideal <M>J</M> of in the image 
##    ring of the ring homomorphism <A>phi</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
##  <#GAPDoc Label="ImagePolynomialRing_DTmanRingHomomorphism_Dat">
##  <ManSection>
##  <Attr Name="ImagePolynomialRing" Arg="phi"/>
##  <Returns>
##    <K>PolynomialRing</K>
##  </Returns>
##  <Description>
##    Returns the polynomial ring which contains the image of the 
##    ring homomorphism <A>phi></A>. All polynomials mapped by <A>phi</A> will 
##    be in this ring.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
