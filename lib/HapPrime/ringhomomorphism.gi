#############################################################################
##
##  HAPPRIME - ringhomomorphism.gi
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
##  <#GAPDoc Label="IsHAPSubringToRingHomomorphismRep_DTmanRingHomomorphismNODOC">
##  <ManSection>
##  <Filt Name="IsHAPSubringToRingHomomorphismRep" Arg="O" Type="Representation"/>
##  <Description>
##  Returns <K>true</K> if the object is in the for a <K>HAPRingHomomorphism</K> 
##  from a subring to a ring, or <K>false</K> otherwise
##  </Returns>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareRepresentation(
  "IsHAPRingToSubringHomomorphismRep",
  IsComponentObjectRep and IsAttributeStoringRep and IsHAPRingHomomorphism,
  [""]
);
# Note this also defines the IsHAPRingHomomorphismGeneralRep filter
#####################################################################

#####################################################################
##  <#GAPDoc Label="IsHAPRingToSubringHomomorphismRep_DTmanRingHomomorphismNODOC">
##  <ManSection>
##  <Filt Name="IsHAPRingToSubringHomomorphismRep" Arg="O" Type="Representation"/>
##  <Description>
##  Returns <K>true</K> if the object is in the for a <K>HAPRingHomomorphism</K> 
##  from a ring to a subring, or <K>false</K> otherwise
##  </Returns>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareRepresentation(
  "IsHAPSubringToRingHomomorphismRep",
  IsComponentObjectRep and IsAttributeStoringRep and IsHAPRingHomomorphism,
  ["eliminationideal"]
);
# Note this also defines the IsHAPRingHomomorphismGeneralRep filter
#####################################################################

#####################################################################
##  <#GAPDoc Label="HAPRingHomomorphismIndeterminateMapRep_DTmanRingHomomorphismNODOC">
##  <ManSection>
##  <Filt Name="HAPRingHomomorphismIndeterminateMapRep" Arg="O" 
##    Type="Representation"/>
##  <Description>
##  Returns <K>true</K> if the object is in the indeterminate map representation 
##  used for a <K>HAPRingHomomorphism</K>, or <K>false</K> otherwise
##  </Returns>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareRepresentation(
  "IsHAPRingHomomorphismIndeterminateMapRep",
  IsComponentObjectRep and IsAttributeStoringRep and IsHAPRingHomomorphism,
  ["Rindetnums", "Sindetnums"]
);
# Note this also defines the IsHAPRingHomomorphismIndeterminateMapRep filter
#####################################################################


#####################################################################
##  <#GAPDoc Label="IsHAPRingReductionHomomorphismRep_DTmanRingHomomorphismNODOC">
##  <ManSection>
##  <Filt Name="IsHAPRingReductionHomomorphismRep" Arg="O" 
##    Type="Representation"/>
##  <Description>
##  Returns <K>true</K> if the object is in the representation 
##  used for a <K>IsHAPRingReductionHomomorphismRep</K>, or <K>false</K> otherwise
##  </Returns>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareRepresentation(
  "IsHAPRingReductionHomomorphismRep",
  IsComponentObjectRep and IsAttributeStoringRep and IsHAPRingHomomorphism,
  ["eliminationideal", "ringindetmap"]
);
# Note this also defines the IsHAPRingHomomorphismIndeterminateMapRep filter
#####################################################################


#####################################################################
##  <#GAPDoc Label="IsHAPZeroRingHomomorphismRep_DTmanRingHomomorphismNODOC">
##  <ManSection>
##  <Filt Name="IsHAPZeroRingHomomorphismRep" Arg="O" 
##    Type="Representation"/>
##  <Description>
##  Returns <K>true</K> if the object is in the representation 
##  used for a <K>IsHAPZeroRingHomomorphismRep</K>, or <K>false</K> otherwise
##  </Returns>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareRepresentation(
  "IsHAPZeroRingHomomorphismRep",
  IsComponentObjectRep and IsAttributeStoringRep and IsHAPRingHomomorphism,
  []
);
# Note this also defines the IsHAPRingHomomorphismIndeterminateMapRep filter
#####################################################################


#####################################################################
##  <#GAPDoc Label="ViewObj_DTmanRingHomomorphismNODOC">
##  <ManSection>
##  <Meth Name="ViewObj" Arg="phi" Label="for HAPRingHomomorphism"/>
##
##  <Description>
##  Prints a short description of the ring homomorphism <A>phi</A>. This is the 
##  usual description printed by &GAP;.
##  </Description>
##  </ManSection>
##  <Log><![CDATA[
##  gap> View(d);
##  <Ring homomorphism>
##  ]]></Log>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  ViewObj,
  "for HAPRingHomomorphism",
  [IsHAPRingHomomorphism],
  function(obj)
    Print("<Ring homomorphism>");
  end
);
#####################################################################

#####################################################################
##  <#GAPDoc Label="PrintObj_DTmanRingHomomorphismNODOC">
##  <ManSection>
##  <Meth Name="PrintObj" Arg="phi" Label="for HAPRingHomomorphism"/>
##
##  <Description>
##  Prints a detailed description of the ring homomorphism <A>phi</A>.
##  </Description>
##  </ManSection>
##  <Log><![CDATA[
##  ]]></Log>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  PrintObj,
  "for IsHAPRingToSubringHomomorphismRep",
  [IsHAPRingToSubringHomomorphismRep],
  function(obj)
    Print("HAPRingToSubringHomomorphism(", SourcePolynomialRing(obj), ", ",
    SourceRelations(obj), ", ", ImageGenerators(obj), ")");
  end
);
#####################################################################
InstallMethod(
  PrintObj,
  "for IsHAPRingHomomorphismIndeterminateMapRep",
  [IsHAPRingHomomorphismIndeterminateMapRep],
  function(obj)
    Print("HAPRingHomomorphismByIndeterminateMap(", SourcePolynomialRing(obj), ", ",
    SourceRelations(obj), ", ", ImagePolynomialRing(obj), ")");
  end
);
#####################################################################
InstallMethod(
  PrintObj,
  "for IsHAPRingHomomorphism",
  [IsHAPRingHomomorphism],
  function(obj)
    # this one also does for HAPRingReductionHomomorphism
    Print("HAPSubringToRingHomomorphism(", SourceGenerators(obj), ", ",
    ImagePolynomialRing(obj), ", ", ImageRelations(obj), ")");
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="Display_DTmanRingHomomorphismNODOC">
##  <ManSection>
##  <Meth Name="Display" Arg="phi" Label="for HAPRingHomomorphism"/>
##
##  <Description>
##  Displays the ring homomorphism <A>phi</A> in a human-readable form. 
##  </Description>
##  </ManSection>
##  <Log><![CDATA[
##  ]]></Log>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  Display,
  "for HAPRingHomomorphism",
  [IsHAPRingHomomorphism],
  function(obj)
    local i;
    Print("Ring homomorphism\n");
    for i in [1..Length(SourceGenerators(obj))] do
      Print("  ", SourceGenerators(obj)[i], " -> ", ImageGenerators(obj)[i], "\n");
    od; 
    if not IsEmpty(SourceRelations(obj)) then
      Print("with relations\n");
      Print("  ", SourceRelations(obj), "\n");
    fi;
    if not IsEmpty(ImageRelations(obj)) then
      if IsEmpty(SourceRelations(obj)) then
        Print("with relations\n");
      else
        Print("and\n");
      fi;
      Print("  ", ImageRelations(obj), "\n");
    fi;
  end
  );
#####################################################################



#####################################################################
##  <#GAPDoc Label="HAPRingToSubringHomomorphism_DTmanRingHomomorphism_Con">
##  <ManSection>
##  <Oper Name="HAPRingToSubringHomomorphism" Arg="Rring, Rrels, Simages"/>
##
##  <Returns>
##    <K>HAPRingHomomorphism</K>
##  </Returns>
##  <Description>
##    Creates a <K>HAPRingHomomorphism</K> which represents the mapping 
##    <M>R/I \to S/J</M>. In this form, <A>Rring</A> a polynomial ring <M>R</M>
##    and <A>Rrels</A> an ideal <M>I</M> in that ring. The image of the indeterminates
##    of <A>R</A> under this mapping are given in <A>Simages</A> and generate
##    the ring <M>S</M>, while the relations <A>Rrels</A> are mapped to 
##    generate <M>J</M>. The ring <M>S</M> may be a subring of the full polynomial
##    ring in its indeterminates.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(HAPRingToSubringHomomorphism,
  [IsPolynomialRing, IsHomogeneousList, 
  IsHomogeneousList and IsRationalFunctionCollection],
  function(Rring, Rrels, Simages)
    local Sring, Rindets, Sindets, phi;

    # Simages must correspond to the indeterminates of Rring
    if Length(IndeterminatesOfPolynomialRing(Rring)) <> Length(Simages) then
      Error("the <Simages> list must be the same length as the number if indeterminates in <Rring> list");
    fi;

    # Check Rrels are in Rring
    if not ForAll(Rrels, i->i in Rring) then
      Error("<Rgens> must be in <Rring>");
    fi;
    # Sring has the same coefficents ring as Rring, but the indeterminates 
    # from Simages
    Sring := PolynomialRing(CoefficientsRing(Rring),
      Reversed(Set(Flat(List(Simages, IndeterminatesOfPolynomial)))));
    # Get the indeterminates in Rgens and Simages and check that they're 
    # distinct
    Rindets := IndeterminatesOfPolynomialRing(Rring);
    Sindets := IndeterminatesOfPolynomialRing(Sring);
    if not IsEmpty(Intersection(Rindets, Sindets)) then
      Error("<Simages> must be from a different ring to <Rring>");
    fi;

    # Make sure that the relations are a Groebner basis 
    # Set the term ordering and then find the Groebner Basis
    SetTermOrdering(Rring, "dp");
    Rrels := SingularReducedGroebnerBasis(Ideal(Rring, Rrels));

    # Create the object 
    phi := Objectify( 
      NewType(NewFamily("HAPRingHomomorphismFamily"), 
        IsHAPRingHomomorphism and IsHAPRingToSubringHomomorphismRep), 
        rec());

    # And remember the generators and so on
    SetSourcePolynomialRing(phi, Rring);
    SetSourceGenerators(phi, Rindets);
    SetImagePolynomialRing(phi, Sring);
    SetImageGenerators(phi, Simages);
    SetImageRelations(phi, ImageOfRingHomomorphism(phi, Rrels));
    # set the source relations after finding the image of them, otherwise the
    # image will simply be killed by the source relations
    SetSourceRelations(phi, Rrels);

    return phi;
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="HAPSubringToRingHomomorphism_DTmanRingHomomorphism_Con">
##  <ManSection>
##  <Heading>HAPSubringToRingHomomorphism</Heading>
##  <Oper Name="HAPSubringToRingHomomorphism" Arg="Rgens, Rrels, Sring" 
##    Label="for relations defined at source"/>
##  <Oper Name="HAPSubringToRingHomomorphism" Arg="Rgens, Sring, Srels" 
##    Label="for relations defined at image"/>
##
##  <Returns>
##    <K>HAPRingHomomorphism</K>
##  </Returns>
##  <Description>
##    Creates a <K>HAPRingHomomorphism</K> which represents the mapping 
##    <M>R/I \to S/J</M>. The ring <M>R</M> is generated by a set of 
##    polynomials <A>Rgens</A> (so <M>R</M> may be a subring of the full 
##    polynomial ring in its indeterminates). The images of <A>Rgens</A> under 
##    the mapping are the indeterminates of the polynomial ring given in 
##    <A>Sring</A>. The ideals can be specified either as a set of relations 
##    <A>Srels</A> in the target ring <M>S</M>, or as a set of relations 
##    <A>Rrels</A> in the source ring. In this second case, <A>Rrels</A> can be
##    polynomials in the full polynomial ring, in which case the ideal <M>I</M>
##    is the intersection of the ideal they generate in the full ring with the
##    subring generated by <A>Rgens</A>. In both cases, the specified ideal is
##    mapped with the homomorphism (or its inverse) to find the corresponding 
##    ideal in the other ring.
##    <P/>
##    This ring homomorphism uses Gröbner bases to perform the mapping, and 
##    the time taken to calculate the basis in this function can be influenced
##    by the choice of monomial ordering. See 
##    <Ref Subsect="RingHomEliminationOrdering"/> for more details.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(HAPSubringToRingHomomorphism,
  [IsHomogeneousList and IsRationalFunctionCollection, 
  IsHomogeneousList, IsPolynomialRing],
  function(Rgens, Rrels, Sring)
    local Rring, Rindets, Sindets, elimring, rels, phi;

    # Simages must correspond to the images of Rgens
    if Length(Rgens) <> Length(IndeterminatesOfPolynomialRing(Sring)) then
      Error("the <Rgens> list must be the same length as the number of indeterminates in <Sring>");
    fi;

    # Find Rring 
    # We may have relations in Rrels that aren't generated by Rgens (and 
    # may even feature more indeterminates then Rgens have), but that is OK
    # in some cases, so we shan't check for that. We'll just find the smallest
    # polynomial ring that contains both Rgens and Rrels.
    Rring := PolynomialRing(CoefficientsRing(Sring),
      Reversed(Set(Flat(List(Concatenation(
        Rgens, Rrels), IndeterminatesOfPolynomial)))));

    # Get the indeterminates in Rgens and Simages and check that they're 
    # distinct
    Rindets := IndeterminatesOfPolynomialRing(Rring);
    Sindets := IndeterminatesOfPolynomialRing(Sring);
    if not IsEmpty(Intersection(Rindets, Sindets)) then
      Error("<Rgens> and must be from a different ring from <Sring>");
    fi;

    # Now make the elimination ordering that we want
    elimring := HAPPRIME_MakeEliminationOrdering(
      CoefficientsRing(Sring), Rindets, Sindets);

    # Add the map relations to the input relations
    # and find a Groebner basis. With the elimination ordering this will
    # also find the relations in the target ring
    rels := Concatenation(Rrels, Rgens - Sindets);
    # Set the base ring if it is not the same or it doesn't have 
    # the same indeterminate order
    if elimring <> SingularBaseRing or
      IndeterminatesOfPolynomialRing(elimring) <> 
        IndeterminatesOfPolynomialRing(SingularBaseRing) then
      SingularSetBaseRing(elimring);
    fi;
    # now do GroebnerBasis
    rels := SingularReducedGroebnerBasis(Ideal(elimring, rels));

    # Create the object 
    phi := Objectify( 
      NewType(NewFamily("HAPRingHomomorphismFamily"), 
        IsHAPRingHomomorphism and IsHAPSubringToRingHomomorphismRep), 
        rec(eliminationideal := Ideal(elimring, rels)));

    # And remember the generators and so on
    SetSourcePolynomialRing(phi, Rring);
    SetSourceGenerators(phi, Rgens);
    SetSourceRelations(phi, Rrels);
    SetImageGenerators(phi, Sindets);
    SetImagePolynomialRing(phi, Sring);
    SetImageRelations(phi, Filtered(rels, i->i in Sring));

    return phi;
  end
);
#####################################################################
InstallMethod(HAPSubringToRingHomomorphism,
  [IsHomogeneousList and IsRationalFunctionCollection, 
  IsPolynomialRing, IsHomogeneousList],
  function(Rgens, Sring, Srels)
    local Rring, Rindets, Sindets, elimring, rels, ring, phi, invphi;

    # Simages must correspond to the images of Rgens
    if Length(Rgens) <> Length(IndeterminatesOfPolynomialRing(Sring)) then
      Error("the <Rgens> list must be the same length as number of indeterminates in <Sring>");
    fi;

    # Find Rring 
    Rring := PolynomialRing(CoefficientsRing(Sring),
      Reversed(Set(Flat(List(Rgens, IndeterminatesOfPolynomial)))));
    # Check that Srels are in S
    if not ForAll(Srels, i->i in Sring) then
      Error("<Srels> must be in <Sring>");
    fi;
    # Get the indeterminates in Rgens and Simages and check that they're 
    # distinct
    Rindets := IndeterminatesOfPolynomialRing(Rring);
    Sindets := IndeterminatesOfPolynomialRing(Sring);
    if not IsEmpty(Intersection(Rindets, Sindets)) then
      Error("<Rgens> must be from a differnt ring from <Sring>");
    fi;

    # Now make the elimination ordering that we want
    elimring := HAPPRIME_MakeEliminationOrdering(
      CoefficientsRing(Sring), Rindets, Sindets);

    # Add the map relations to the input relations
    # and find a Groebner basis. With the elimination ordering this will
    # also find the relations in the target ring
    rels := Concatenation(Srels, Rgens - Sindets);
    # Set the base ring if it is not the same or it doesn't have 
    # the same indeterminate order
    if elimring <> SingularBaseRing or
      IndeterminatesOfPolynomialRing(elimring) <> 
        IndeterminatesOfPolynomialRing(SingularBaseRing) then
      SingularSetBaseRing(elimring);
    fi;
    # now do GroebnerBasis
    rels := SingularReducedGroebnerBasis(Ideal(elimring, rels));

    # Create the object 
    phi := Objectify( 
      NewType(NewFamily("HAPRingHomomorphismFamily"), 
        IsHAPRingHomomorphism and IsHAPSubringToRingHomomorphismRep), 
        rec(eliminationideal := Ideal(elimring, rels)));

    # And remember the generators and so on
    SetImageGenerators(phi, Sindets);
    SetImagePolynomialRing(phi, Sring);
    # Srels may not have been a GroebnerBasis, but it is now
    SetImageRelations(phi, Filtered(rels, i->i in Sring));
    SetSourcePolynomialRing(phi, Rring);
    SetSourceGenerators(phi, Rgens);
    # Map the image relations back to Rring to get the source relations
    # Create the inverse by hand rather than using InverseRingHomomorphism
    # since phi isn't finished yet, the (partial) inverse would be stored
    # as an attribute otherwise
    invphi := HAPRingToSubringHomomorphism(Sring, [], Rgens);
    SetSourceRelations(phi, ImageOfRingHomomorphism(invphi, ImageRelations(phi)));

    return phi;
  end
);
#####################################################################

#####################################################################
##  <#GAPDoc Label="HAPRingHomomorphismByIndeterminateMap_DTmanRingHomomorphism_Con">
##  <ManSection>
##  <Oper Name="HAPRingHomomorphismByIndeterminateMap" Arg="R, Rrels, S"/>
##
##  <Returns>
##    <K>HAPRingHomomorphism</K>
##  </Returns>
##  <Description>
##    Creates a <K>HAPRingHomomorphism</K> which represents the map
##    <M>R/I \to S/J </M> which is a simple relabelling of indeterminates: the 
##    image of the <M>i</M>th indeterminate of <A>R</A> under the mapping is 
##    taken to be the <M>i</M>th indeterminate of <A>S</A>. The ideal <M>I</M>
##    is generated by <A>Rrels</A> and are mapped using the homomorphism to
##    generate <M>J</M>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(HAPRingHomomorphismByIndeterminateMap,
  [IsPolynomialRing, IsHomogeneousList, IsPolynomialRing],
  function(R, Rrels, S)
    local phi, Rindetnums, Sindetnums;
     
    # Check the inputs
    if CoefficientsRing(R) <> CoefficientsRing(S) then
      Error("<R> and <S> must have the same coefficient ring");
    fi;
    if Length(IndeterminatesOfPolynomialRing(R)) <> 
      Length(IndeterminatesOfPolynomialRing(S)) then 
      Error("<R> and <S> must have the same number of indeterminates");
    fi;
    if not ForAll(Rrels, i->i in R) then
      Error("<Rgens> must be in <R>");
    fi;

    Rindetnums := List(IndeterminatesOfPolynomialRing(R), 
      IndeterminateNumberOfUnivariateRationalFunction);
    Sindetnums := List(IndeterminatesOfPolynomialRing(S), 
      IndeterminateNumberOfUnivariateRationalFunction);

    #################
    # Create the object 
    phi := Objectify( 
      NewType(NewFamily("HAPRingHomomorphismFamily"), 
        IsHAPRingHomomorphism and IsHAPRingHomomorphismIndeterminateMapRep), 
        rec(Rindetnums := Rindetnums, Sindetnums := Sindetnums));

    # And remember the generators and so on
    SetSourcePolynomialRing(phi, R);
    SetSourceGenerators(phi, IndeterminatesOfPolynomialRing(R));
    SetSourceRelations(phi, Rrels);
    SetImagePolynomialRing(phi, S);
    SetImageGenerators(phi, IndeterminatesOfPolynomialRing(S));
    # Use the object itself to compute the image relations
    SetImageRelations(phi, ImageOfRingHomomorphism(phi, Rrels));

    return phi;
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="HAPRingReductionHomomorphism_DTmanRingHomomorphism_Con">
##  <ManSection>
##  <Oper Name="HAPRingReductionHomomorphism" Arg="R, Rrels[, avoid]"
##    Label="for ring presentation"/>
##  <Oper Name="HAPRingReductionHomomorphism" Arg="phi[, avoid]"
##    Label="for image of ring homomorphism"/>
##  <Returns>
##    <K>HAPRingHomomorphism</K>
##  </Returns>
##  <Description>
##    For a polynomial ring <A>R</A> and ideal <M>I</M> generated by 
##    <A>Rrels</A>, this function finds an isomorphic ring in fewer 
##    indeterminates (or the same number, if this is not possible). This new 
##    ring will avoid the indeterminates of <A>R</A> and any further 
##    indeterminates listed in <A>avoid</A>. The function returns the map 
##    between <M>R/I</M> and the new ring.
##    <P/>
##    In the second form, this function reduces the target ring of the 
##    ring homomorphism <A>phi</A> and returns the map between this
##    and the reduced ring. This map will also avoid the indeterminates in the
##    source ring of <A>phi</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallOtherMethod(HAPRingReductionHomomorphism,
  [IsPolynomialRing, IsHomogeneousList],
  function(R, Rrels)
    return HAPRingReductionHomomorphism(R, Rrels, []);
  end
);    
#####################################################################
InstallMethod(HAPRingReductionHomomorphism,
  [IsPolynomialRing, IsHomogeneousList, IsHomogeneousList],
  function(R, Rrels, avoid)
    local ideal, allindets, indets, removedindets, eliminationrelations, indet, 
    i, t, t2, unimons, relation, elimring, elimorder, len, poly, newring, 
    ringindetmap, images, phi, remainingR;
     
    # Check the inputs
    if not ForAll(Rrels, i->i in R) then
      Error("<Rrels> must be in <R>");
    fi;

    # Check for unit in the ideal - if so, we return a trivial ring
    if Length(Rrels) = 1 and IsOne(Rrels[1]) then
      return HAPZeroRingHomomorphism(R, Rrels);
    fi;

    ideal := ShallowCopy(Rrels);;
    # indeterminates will be removed from indets and added to removedindets
    # as we do our reduction 
    indets := ShallowCopy(IndeterminatesOfPolynomialRing(R));
    allindets := ShallowCopy(IndeterminatesOfPolynomialRing(R));
    removedindets := [];
    eliminationrelations := [];
    elimring := PolynomialRing(GF(2), allindets);
    elimorder := MonomialLexOrdering(allindets);


    # Are any of the leading terms of the ideal single indeterminates?
    # If so then we can (after reduction) remove those indeterminates and those 
    # terms of the ideal
    repeat
      indet := false;
      # Find a relation in the ideal that involves a solitary indeterminate
      for i in [1..Length(ideal)] do
        for t in TermsOfPolynomial(ideal[i]) do
          unimons := UnivariateMonomialsOfMonomial(t[1]);
          if Length(unimons) = 1 
            and unimons[1] = 
              IndeterminateOfUnivariateRationalFunction(unimons[1]) then
            # This term involves a solitary indeterminate. 
            # Check that no other term in this relation also involves this 
            # indeterminate
            indet := unimons[1];
            for t2 in TermsOfPolynomial(ideal[i] - t[1]*t[2]) do
              if IsOne(DenominatorOfRationalFunction(t2[1] / indet)) then
                indet := false;
                break;
              fi;
            od;
            if indet <> false then
              # We're OK - no other term in this relation involves this
              # indeterminate, so we can go on to remove it
              break;
            fi;
          fi;
        od;
        # Have we found a solitary indeterminate?
        if indet <> false then
          relation := Remove(ideal, i);
          Add(eliminationrelations, relation);
          break;
        fi;
      od;

      if indet <> false then
        # Remove this indeterminate from the indets list and put it
        # onto the removedindets list instead
        Remove(indets, Position(indets, indet));
        Add(removedindets, indet);

        # The elimination order has the removed indeterminates first 
        elimring := PolynomialRing(GF(2), allindets);
        elimorder := MonomialLexOrdering(Concatenation(removedindets, indets));
        SetTermOrdering(elimring, elimorder);

        # And make sure that our relation has a unit coefficient
        relation := relation / 
          LeadingCoefficientOfPolynomial(relation, elimorder);
        # need to set the base ring again since I've changed the term ordering
        SingularSetBaseRing(elimring); 
        SingularSetNormalFormIdealNC(Ideal(elimring, [relation]));

        # Now reduce all the other relations in the ideal with this one, which
        # will get rid of this indeterminate
        i := 1;
        len := Length(ideal);
        while i <= len do
          poly := SingularPolynomialNormalForm(ideal[i]);
          if IsZero(poly) then
            ideal[i] := ideal[len];
            Unbind(ideal[len]);
            len := len - 1;
          else
            ideal[i] := poly;
            i := i + 1;
          fi;
        od;
      fi;
    until indet = false;
  
    # Tidy up the elimination relations so that they will all
    # give something in the remaining indeterminates
    eliminationrelations := SingularReducedGroebnerBasis(
      Ideal(elimring, eliminationrelations));

    # make sure the remaining indets are a Groebner Basis
    remainingR := PolynomialRing(CoefficientsRing(R), indets);
    ideal := SingularReducedGroebnerBasis(Ideal(remainingR, ideal));

    # Now make the new ring and map to that ring
    newring := PolynomialRing(CoefficientsRing(R), Length(indets),
      Concatenation(avoid, IndeterminatesOfPolynomialRing(R)));
    ringindetmap := HAPRingHomomorphismByIndeterminateMap(
      remainingR, ideal, newring);

    # Finally, work out what the images of our original indeterminates are
    images := [];  
    SingularSetBaseRing(elimring); 
    SingularSetNormalFormIdealNC(Ideal(elimring, eliminationrelations));
    for i in IndeterminatesOfPolynomialRing(R) do
      if i in indets then
        Add(images, 
          IndeterminatesOfPolynomialRing(newring)[Position(indets, i)]);  
      else
        Add(images, 
          ImageOfRingHomomorphism(ringindetmap, 
            SingularPolynomialNormalForm(i)));
      fi;
    od;

    ################
    # Create the object 
    phi := Objectify( 
      NewType(NewFamily("HAPRingHomomorphismFamily"), 
        IsHAPRingHomomorphism and IsHAPRingReductionHomomorphismRep), 
        rec(eliminationideal := Ideal(elimring, eliminationrelations),
          ringindetmap := ringindetmap));

    # And remember the generators and so on
    SetSourcePolynomialRing(phi, R);
    SetSourceGenerators(phi, IndeterminatesOfPolynomialRing(R));
    SetSourceRelations(phi, Rrels);
    SetImagePolynomialRing(phi, newring);
    SetImageGenerators(phi, images);
    SetImageRelations(phi, ImageRelations(ringindetmap));

    return phi;
  end
);
#####################################################################
InstallOtherMethod(HAPRingReductionHomomorphism,
  [IsHAPRingHomomorphism],
  function(phi)
    return HAPRingReductionHomomorphism(phi, []);
  end
);
#####################################################################
InstallMethod(HAPRingReductionHomomorphism,
  [IsHAPRingHomomorphism, IsHomogeneousList],
  function(phi, avoid)
    return HAPRingReductionHomomorphism(
      ImagePolynomialRing(phi), ImageRelations(phi), 
      Concatenation(
        IndeterminatesOfPolynomialRing(SourcePolynomialRing(phi)), avoid));
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="HAPZeroRingHomomorphism_DTmanRingHomomorphism_Con">
##  <ManSection>
##  <Oper Name="HAPZeroRingHomomorphism" Arg="R, Rrels"/>
##
##  <Returns>
##    <K>HAPRingHomomorphism</K>
##  </Returns>
##  <Description>
##    Creates a <K>HAPRingHomomorphism</K> which maps from the ring <M>R</M>, 
##    with an ideal generated by <A>Rrels</A>, into the trival ring generated
##    by zero.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(HAPZeroRingHomomorphism,
  [IsPolynomialRing, IsHomogeneousList],
  function(R, Rrels)
    local Sring, Rindets, Sindets, elimring, rels, ring, phi;

    if not ForAll(Rrels, i->i in R) then
      Error("<Rgens> must be in <R>");
    fi;

    # Create the object 
    phi := Objectify( 
      NewType(NewFamily("HAPRingHomomorphismFamily"), 
        IsHAPRingHomomorphism and IsHAPZeroRingHomomorphismRep), 
        rec());

    # And remember the generators and so on
    SetSourcePolynomialRing(phi, R);
    SetSourceGenerators(phi, IndeterminatesOfPolynomialRing(R));
    SetSourceRelations(phi, Rrels);
    SetImageGenerators(phi, ListWithIdenticalEntries(
      Length(IndeterminatesOfPolynomialRing(R)), Zero(R)));
    SetImagePolynomialRing(phi, Ring(Zero(R)));
    SetImageRelations(phi, []);

    return phi;
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="InverseRingHomomorphism_DTmanRingHomomorphism_Con">
##  <ManSection>
##  <Attr Name="InverseRingHomomorphism" Arg="phi"/>
##
##  <Returns>
##  <K>HAPRingHomomorphism</K>
##  </Returns>
##  <Description>
##  Returns (as a ring homomorphism) the inverse of the ring homomorphism 
##  <A>phi</A>.
##  <P/>
##  If the inverse homomorphism requires an elimination Gröbner basis to 
##  perform the mapping (for example when computing the inverse of a 
##  <K>HAPRingHomomorphism</K> constructed with 
##  <Ref Func="HAPRingToSubringHomomorphism"/>) then the ordering can be 
##  specified using the options stack. See 
##  <Ref Subsect="RingHomEliminationOrdering"/> for more details.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(InverseRingHomomorphism, 
  [IsHAPSubringToRingHomomorphismRep],  
  function(phi)
    local inverse;
    inverse := HAPRingToSubringHomomorphism(
      ImagePolynomialRing(phi), ImageRelations(phi), SourceGenerators(phi));
    # We can also remember the inverse!
    SetInverseRingHomomorphism(inverse, phi);
    return inverse;  
  end
);
#####################################################################
InstallMethod(InverseRingHomomorphism, 
  [IsHAPRingToSubringHomomorphismRep],
  function(phi)
    local inverse;
    inverse := HAPSubringToRingHomomorphism(
      ImageGenerators(phi), ImageRelations(phi), SourcePolynomialRing(phi));
    # We can also remember the inverse!
    SetInverseRingHomomorphism(inverse, phi);
    return inverse;  
  end
);
#####################################################################
InstallMethod(InverseRingHomomorphism, 
  [IsHAPRingReductionHomomorphismRep],
  function(phi)
    local inverse;
    inverse := InverseRingHomomorphism(phi!.ringindetmap);
    # We can also remember the inverse!
    SetInverseRingHomomorphism(inverse, phi);
    return inverse;  
  end
);
#####################################################################
InstallMethod(InverseRingHomomorphism, 
  [IsHAPRingHomomorphismIndeterminateMapRep],
  function(phi)
    local inverse;

    # Just swap over source and target
    inverse := Objectify( 
      NewType(NewFamily("HAPRingHomomorphismFamily"), 
        IsHAPRingHomomorphism and IsHAPRingHomomorphismIndeterminateMapRep), 
        rec(Rindetnums := phi!.Sindetnums, Sindetnums := phi!.Rindetnums));

    # And remember the generators and so on
    SetSourcePolynomialRing(inverse, ImagePolynomialRing(phi));
    SetSourceGenerators(inverse, ImageGenerators(phi));
    SetSourceRelations(inverse, ImageRelations(phi));
    SetImagePolynomialRing(inverse, SourcePolynomialRing(phi));
    SetImageGenerators(inverse, SourceGenerators(phi));
    SetImageRelations(inverse, SourceRelations(phi));
    # We can also remember the inverse!
    SetInverseRingHomomorphism(inverse, phi);

    return inverse;
  end
);
#####################################################################
InstallMethod(InverseRingHomomorphism, 
  [IsHAPZeroRingHomomorphismRep],
  function(phi)
    return fail;
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="CompositionRingHomomorphism_DTmanRingHomomorphism_Con">
##  <ManSection>
##  <Oper Name="CompositionRingHomomorphism" Arg="phiA, phiB"/>
##
##  <Returns>
##    <K>HAPRingHomomorphism</K>
##  </Returns>
##  <Description>
##    Returns the ring homomorphism that is the composition of the ring 
##    homomorphisms <A>phiA</A> and <A>phiB</A>. The source ring of <A>phiB</A>
##    must be in the image ring of <A>phiA</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(CompositionRingHomomorphism, 
  [IsHAPRingToSubringHomomorphismRep, IsHAPRingHomomorphism],
  function(phiA, phiB)

    # Check that phiA and phiB are compatible
    HAPPRIME_RingHomomorphismsAreComposable(phiA, phiB);

    # The indeterminates of the source of phiA map to ImageGenerators(phiA).
    # Map those on through phiB to find the image of the composed homomorphism
    # and simply reuse the source relations (which will be mapped through)
    return HAPRingToSubringHomomorphism(
      SourcePolynomialRing(phiA), 
      SourceRelations(phiA), 
      List(ImageGenerators(phiA), i->ImageOfRingHomomorphism(phiB, i)));
  end
);
#####################################################################
InstallMethod(CompositionRingHomomorphism, 
  [IsHAPSubringToRingHomomorphismRep, IsHAPRingToSubringHomomorphismRep],
  function(phiA, phiB)

    Error("Can't compose a HAPSubringToRingHomomorphism with a HAPRingToSubringHomomorphism");
  end
);
#####################################################################
InstallMethod(CompositionRingHomomorphism, 
  [IsHAPSubringToRingHomomorphismRep, IsHAPRingHomomorphism],
  function(phiA, phiB)
    local gens;

    # Check that phiA and phiB are compatible
    HAPPRIME_RingHomomorphismsAreComposable(phiA, phiB);

    # Map the image indeterminates of phiB through the two homomorphisms
    gens := List(
      IndeterminatesOfPolynomialRing(ImagePolynomialRing(phiB)),
      i->PreimageOfRingHomomorphism(phiA, PreimageOfRingHomomorphism(phiB, i)));

    # and reuse the relations from the image of phiB
    return HAPSubringToRingHomomorphism(
      gens, 
      SourceRelations(phiA), 
      ImagePolynomialRing(phiB));
  end
);
#####################################################################
InstallMethod(CompositionRingHomomorphism, 
  [IsHAPRingHomomorphismIndeterminateMapRep, IsHAPRingHomomorphismIndeterminateMapRep],
  function(phiA, phiB)

    # Check that phiA and phiB are compatible
    HAPPRIME_RingHomomorphismsAreComposable(phiA, phiB);

    # The indeterminates of the source of phiA map to the indeterminates of
    # ImagePolynomialRing(phiA).
    # Map those on through phiB to find the image of the composed homomorphism
    # and simply reuse the source relations (which will be mapped through)
    return HAPRingHomomorphismByIndeterminateMap(
      SourcePolynomialRing(phiA), 
      SourceRelations(phiA), 
      PolynomialRing(CoefficientsRing(SourcePolynomialRing(phiA)),
        List(ImageGenerators(phiA), i->ImageOfRingHomomorphism(phiB, i))));
  end
);
#####################################################################
InstallMethod(CompositionRingHomomorphism, 
  [IsHAPRingHomomorphism, IsHAPZeroRingHomomorphismRep],
  function(phiA, phiB)

    # Check that phiA and phiB are compatible
    HAPPRIME_RingHomomorphismsAreComposable(phiA, phiB);

    return HAPZeroRingHomomorphism(
      SourcePolynomialRing(phiA), 
      SourceRelations(phiA));
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="ImageOfRingHomomorphism_DTmanRingHomomorphism_Gen">
##  <ManSection>
##  <Heading>ImageOfRingHomomorphism</Heading>
##  <Oper Name="ImageOfRingHomomorphism" Arg="phi, poly" 
##    Label="for one polynomial"/>
##  <Oper Name="ImageOfRingHomomorphism" Arg="phi, coll" 
##    Label="for collection of polynomials"/>
##
##  <Returns>
##    Polynomial or list
##  </Returns>
##  <Description>
##    Returns the image of the polynomial <A>poly</A> under the ring 
##    homomorphism <A>phi</A>. The input must be an element(s) of the 
##    source ring of <A>phi</A> (see <Ref Attr="SourcePolynomialRing"/>).
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(ImageOfRingHomomorphism, 
  [IsHAPRingToSubringHomomorphismRep, 
    IsHomogeneousList and IsRationalFunctionCollection],
  function(phi, coll)
    local imcoll, p, i, im, t, term, um, iexp;

    # Check the input is valid
    for p in coll do
      if not p in SourcePolynomialRing(phi) then
        Error("polynomials must be from the source ring of <phi>");
      fi;
    od;
        
    # We should always have SourceRelations unless we are constructing
    # one of these
    if HasSourceRelations(phi) and (not IsEmpty(SourceRelations(phi))) then
      SingularSetNormalFormIdealNC(
        Ideal(SourcePolynomialRing(phi), SourceRelations(phi)));

      imcoll := List(coll, SingularPolynomialNormalForm);
    else
      imcoll := ShallowCopy(coll);
    fi;
    
    for i in [1..Length(imcoll)] do
      # Now map to the image ring
      im := Zero(imcoll[1]);
      for t in TermsOfPolynomial(imcoll[i]) do
        term := t[2];
        for um in UnivariateMonomialsOfMonomial(t[1]) do
          iexp := IndeterminateAndExponentOfUnivariateMonomial(um);
          term := term * ImageGenerators(phi)[
            Position(SourceGenerators(phi), iexp[1])]^iexp[2];
        od;
        im := im + term;
      od;  
      imcoll[i] := im;
    od;
    return imcoll;
  end
);
#####################################################################
InstallMethod(ImageOfRingHomomorphism, 
  [IsHAPSubringToRingHomomorphismRep, 
    IsHomogeneousList and IsRationalFunctionCollection],
  function(phi, coll)
    local i, imcoll, p;

     SingularSetNormalFormIdealNC(phi!.eliminationideal);
 
     imcoll := [];
     for p in coll do
       if not p in SourcePolynomialRing(phi) then
         Error("polynomials must be from the source ring of <phi>");
       fi;
       i := SingularPolynomialNormalForm(p);
       Add(imcoll, i);
     od;
     return imcoll;
  end
);
#####################################################################
InstallMethod(ImageOfRingHomomorphism, 
  [IsHAPRingHomomorphismIndeterminateMapRep, IsHomogeneousList and IsRationalFunctionCollection],
  function(phi, coll)
    local newpolys, fam, p, extrep, extrepnew, i, e, j;

    #newpolys := EmptyPlist(Length(coll));
newpolys:=[];
    fam := FamilyObj(coll[1]);
    for p in coll do
      if not p in SourcePolynomialRing(phi) then
        Error("<polys> must be a polynomial or a list of polynomials in the source ring of <phi>");
      fi;

      extrep := ExtRepPolynomialRatFun(p);
      #extrepnew := EmptyPlist(Length(extrep));
extrepnew := [];
      i := 1;
      while i < Length(extrep) do
        e := ShallowCopy(extrep[i]);
        if not IsEmpty(e) then
          j := 1;
          repeat
            # Swap the indeterminate number
            e[j] := phi!.Sindetnums[Position(phi!.Rindetnums, e[j])];
            j := j+2;
          until j > Length(e);
        fi;
        Add(extrepnew, e);
        Add(extrepnew, extrep[i+1]);
        i := i+2;
      od;
      Add(newpolys, PolynomialByExtRep(fam, extrepnew));
    od;

    return newpolys;
  end
);
#####################################################################
InstallMethod(ImageOfRingHomomorphism, 
  [IsHAPRingReductionHomomorphismRep, IsHomogeneousList and IsRationalFunctionCollection],
  function(phi, coll)

    # Reduce the polynomials with the eliminationrelations and elimorder
    # and then map using the indeterminate map
    SingularSetNormalFormIdealNC(phi!.eliminationideal);
    return List(coll,
      p->ImageOfRingHomomorphism(phi!.ringindetmap, 
        SingularPolynomialNormalForm(p)));
  end
);
#####################################################################
InstallMethod(ImageOfRingHomomorphism, 
  [IsHAPZeroRingHomomorphismRep, IsHomogeneousList and IsRationalFunctionCollection],
  function(phi, coll)
    return List(coll, i->Zero(i));
  end
);
#####################################################################
InstallOtherMethod(ImageOfRingHomomorphism, 
  [IsHAPRingHomomorphism, IsEmpty],
  function(phi, coll)
    return [];
  end
);
#####################################################################
InstallOtherMethod(ImageOfRingHomomorphism, 
  [IsHAPRingHomomorphism, IsRationalFunction],
  function(phi, poly)
    return ImageOfRingHomomorphism(phi, [poly])[1];
  end
);
#####################################################################



#####################################################################
##  <#GAPDoc Label="PreimageOfRingHomomorphism_DTmanRingHomomorphism_Gen">
##  <ManSection>
##  <Heading>PreimageOfRingHomomorphism</Heading>
##  <Oper Name="PreimageOfRingHomomorphism" Arg="phi, poly" 
##    Label="for one polynomial"/>
##  <Oper Name="PreimageOfRingHomomorphism" Arg="phi, coll" 
##    Label="for collection of polynomials"/>
##
##  <Returns>
##    Polynomial or list
##  </Returns>
##  <Description>
##    Returns the preimage of the polynomial <A>poly</A> under the ring 
##    homomorphism <A>phi</A>. The input must be an element(s) of the 
##    image ring of <A>phi</A> (see <Ref Attr="ImagePolynomialRing"/>).
##    This function is a synonym for 
##    <C>ImageOfRingHomomorphism(InverseRingHomomorphism(phi), poly)</C>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(PreimageOfRingHomomorphism, 
  [IsHAPRingHomomorphism, IsHomogeneousList and IsRationalFunctionCollection],
  function(phi, coll)
    return ImageOfRingHomomorphism(InverseRingHomomorphism(phi), coll);
  end
);
#####################################################################
InstallOtherMethod(PreimageOfRingHomomorphism, 
  [IsHAPRingHomomorphism, IsRationalFunction],
  function(phi, poly)
    return ImageOfRingHomomorphism(InverseRingHomomorphism(phi), poly);
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="HAPPRIME_ShuffleRandomSource_Int">
##  <ManSection>
##  <Var Name="HAPPRIME_ShuffleRandomSource"/>
##
##  <Description>
##  The random source for shuffling the list in 
##  <Ref Func="HAPPRIME_MakeEliminationOrdering"/>. This is a source
##  of type <Ref Func="IsMersenneTwister" BookName="ref"/> which is initialised
##  with a seed of 1 when the package is loaded. For a more random random seed,
##  use <C>HAPPRIME_ShuffleRandomSource := RandomSource(IsMersenneTwister, Runtime());</C>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
HAPPRIME_ShuffleRandomSource := RandomSource(IsMersenneTwister, 1);
#####################################################################
##  <#GAPDoc Label="HAPPRIME_MakeEliminationOrdering_DTmanRingHomomorphism_Int">
##  <ManSection>
##  <Oper Name="HAPPRIME_MakeEliminationOrdering" Arg="coeff, Rindets, Sindets"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  Returns a list <C>[jointring, ord]</C> which defines an ordering to
##  be used by the <Ref Sect="singular" Ref="singular"/> package when
##  performing variable elimination. The indeterminates in the list 
##  <A>Rindets</A> are guaranteed to be greater than the indeterminates in 
##  <A>Sindets</A>, and both have coefficient field <A>coeff</A>.
##  <P/>
##  The precise ordering is determined by the options 
##  <Ref Chap="Options Stack" BookName="ref"/> <C>EliminationIndexOrder</C> and
##  <C>EliminationBlockOrdering</C>.
##  See <Ref Sect="RingHom_EliminationOrdering"/> for details of possible
##  orderings.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_MakeEliminationOrdering,
  function(coeff, Rindets, Sindets)
    local jointring, str, rs, start, ord, i, ShuffleList;

    ########################################
    ShuffleList := function(rs, L)
      local l, Lcopy, LShuffled;
      l := Length(L);
      Lcopy := ShallowCopy(L);
      LShuffled := [];
      repeat
        Add(LShuffled, Remove(Lcopy, 
          Random(rs, 1, l)));
        l := l -1;
      until l = 1;
      Add(LShuffled, Lcopy[1]);
      return LShuffled;
    end;
    ########################################

    # This is the default index order
    jointring := PolynomialRing(coeff, Concatenation(Rindets, Sindets));
    # Order the joint ring with the required index order
    str := ValueOption("EliminationIndexOrder");
    if str = "Forward" or str = fail then 
      # use the default
    elif str = "Reverse" then 
      jointring := PolynomialRing(coeff, 
        Concatenation(Reversed(Rindets), Reversed(Sindets)));
    elif str  = "Shuffle" then 
      jointring := PolynomialRing(coeff, 
        Concatenation(ShuffleList(HAPPRIME_ShuffleRandomSource, Rindets), 
          ShuffleList(HAPPRIME_ShuffleRandomSource, Sindets)));
    elif Length(str) > 7 and str{[1..7]} = "Shuffle" then 
      # This is the LexShuffleIndexXX case. Extract the XX
      rs := Int(str{[8..Length(str)]});
      if rs = fail then
        Error("unrecognised integer at the end of ", 
          ValueOption("EliminationOrdering"));
      else
        rs := RandomSource(IsMersenneTwister, rs);
        jointring := PolynomialRing(coeff, 
          Concatenation(ShuffleList(rs, Rindets), ShuffleList(rs, Sindets)));
      fi;
    else
      Info(InfoWarning, 1, "option EliminationIndexOrder:=", str, 
        " is not a valid option. Using the default Forward order.");
    fi;

    # Make the Singular ordering
    # This is the default ordering
    ord := ["lp", Length(Rindets), "Dp", Length(Sindets)];
    # See if we have any other requests
    str := ValueOption("EliminationBlockOrdering");
    if str = fail then
      # use the default
    elif Length(str) < 3 then 
      Info(InfoWarning, 1, "option EliminationBlockOrdering:=", str, 
        " is not a valid option. Using the default LexGrlex order.");
    else
      # Find the ordering for eliminated variables
      start := 1; # Where does the next option start?
      for i in [1, 3] do
        if Length(str) >= (start+2) and str{[start..(start+2)]} = "Lex" then
          ord[i] := "lp";
          start := 4;
        elif Length(str) >= (start+4) and str{[start..(start+4)]} = "Grlex" then
          ord[i] := "Dp";
          start := 6;
        elif Length(str) >= (start+6) and str{[start..(start+6)]} = "Grevlex" then
          ord[i] := "dp";
          start := 8;
        elif i = 1 then
          Info(InfoWarning, 1, "option EliminationBlockOrdering:=", str, 
            " is not a valid option. Using the default LexGrlex order.");
        fi;
      od;
    fi;
    SetTermOrdering(jointring, ord);
    return jointring;
  end
);
#####################################################################



#####################################################################
##  <#GAPDoc Label="HAPPRIME_RingHomomorphismsAreComposable_DTmanRingHomomorphism_Int">
##  <ManSection>
##  <Oper Name="HAPPRIME_RingHomomorphismsAreComposable" Arg="phiA, phiB"/>
##
##  <Returns>
##    nothing
##  </Returns>
##  <Description>
##  Checks that the ring homomorphisms <A>phiA</A> and <A>phiB</A> are
##  can be composed, i.e. the generators of <A>phiB</A> lie in the image ring of 
##  <A>phiA</A>, and the image relations of <A>phiA</A> are the same as the 
##  source relations of <A>phiB</A>. If the ring homomorphisms are not
##  composable, then an appropriate error is thrown.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_RingHomomorphismsAreComposable,
  function(phiA, phiB)
    local A, B;
    # Check that phiA and phiB are compatible
    if not ForAll(
      SourceGenerators(phiB), i->i in ImagePolynomialRing(phiA)) then
        Error("the generators of the source of <phiB> must be in the image ring of <phiA>");
    fi;
    if Set(ImageRelations(phiA)) <> Set(SourceRelations(phiB)) then
      # If they're not identical, are their GroebnerBasis the same?
      A := ImagePolynomialRing(phiA);
      SetTermOrdering(A, "dp");
      B := SourcePolynomialRing(phiB);
      SetTermOrdering(B, "dp");
      if Set(SingularReducedGroebnerBasis(Ideal(A, ImageRelations(phiA)))) <>
        Set(SingularReducedGroebnerBasis(Ideal(B, SourceRelations(phiB)))) then
          Error("the ideal of the image ring of <phiA> must be the same as the ideal of the source of <phiB>");
      fi;
    fi;
  end
);
#####################################################################
  
