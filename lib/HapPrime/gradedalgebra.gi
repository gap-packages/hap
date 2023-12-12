#############################################################################
##
##  HAPPRIME - gradedalgebrapresentation.gi
##  Functions, Operations and Methods to implement graded algebras
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
##  <#GAPDoc Label="GradedAlgebraPresentationRep_DTmanGradedAlgebraNODOC">
##  <ManSection>
##  <Filt Name="IsGradedAlgebraPresentationRep" Arg="O" Type="Representation"/>
##  <Description>
##  Returns <K>true</K> if the object is in the internal representation used for 
##  a <K>GradedAlgebraPresentation</K>, or <K>false</K> otherwise
##  </Returns>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareRepresentation(
  "IsGradedAlgebraPresentationRep",
  IsComponentObjectRep and IsAttributeStoringRep,
  ["ring", "relations", "degrees"]
);
# Note this also defines the function IsGradedAlgebraPresentationRep
#####################################################################

#####################################################################
# The type for a GradedAlgebraPresentation is a GradedAlgebraPresentation in
# GradedAlgebraPresentationRep representation, in the GradedAlgebraPresentation family
GradedAlgebraPresentationType :=
NewType(NewFamily("GradedAlgebraPresentationFamily"), IsGradedAlgebraPresentation and IsGradedAlgebraPresentationRep);
#####################################################################



#####################################################################
##  <#GAPDoc Label="GradedAlgebraPresentation_DTmanGradedAlgebra_Con">
##  <ManSection Label="GradedAlgebraPresentationConstructors">
##  <Heading>GradedAlgebraPresentation construction functions</Heading>
##  <Oper Name="GradedAlgebraPresentation" Arg="R, I, degs"/>
##  <Oper Name="GradedAlgebraPresentationNC" Arg="R, I, degs"/>
##
##  <Returns>
##  <K>GradedAlgebraPresentation</K>
##  </Returns>
##  <Description>
##  Construct a <K>GradedAlgebraPresentation</K> object representing a presentation
##  of a graded algebra as the quotient of a polynomial ring <A>R</A> by the
##  ideal <A>I</A> (as a list of relations in <A>R</A>) where the indeterminates
##  of <A>R</A> (as returned by <Ref Func="IndeterminatesOfGradedAlgebraPresentation"/> 
##  have degrees <A>degs</A> respectively.
##  <P/>
##  The function <K>GradedAlgebraPresentation</K> checks that the arguments are
##  compatible, while the <C>NC</C> method performs no checks.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(GradedAlgebraPresentation,
  [IsPolynomialRing, IsHomogeneousList, IsHomogeneousList],
  function(ring, relations, degs)
    local i;
    # Check that there is one degree for each ring indeterminant
    if Length(IndeterminatesOfPolynomialRing(ring)) <> Length(degs) then
      Error("the number of degrees in <degs> is not the same as the number of indeterminants in <R>.");
    fi;
    # Check that each element in <relations> is in the ring
    for i in relations do
      if not i in ring then
        Error("the entries in <I> must all be polynomials from <R>.");
      fi;
    od;
    # Check that each entry in <degs> is a positive integer
    #if not ForAll(degs, IsPosInt) then
    if not ForAll(degs, IsInt) then   #Graham changed this October 2023
      Error("the entries in <degs> must all be positive integers");
    fi;
    return GradedAlgebraPresentationNC(ring, relations, degs);
  end
);
#####################################################################
InstallMethod(GradedAlgebraPresentationNC,
  [IsPolynomialRing, IsHomogeneousList, IsHomogeneousList],
  function(ring, relations, degs)
    return Objectify(
      GradedAlgebraPresentationType, 
      rec(ring := ring, relations := relations, degrees := degs)
      );
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="BaseRing_DTmanGradedAlgebra_Dat">
##  <ManSection>
##  <Attr Name="BaseRing" Arg="A"/>
##
##  <Returns>
##    Polynomial ring
##  </Returns>
##  <Description>
##  Returns the base ring of the graded algebra presentation <A>A</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(BaseRing,
  "of graded algebra presentation",
  [IsGradedAlgebraPresentation],
  function(A)
    return A!.ring;
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="CoefficientsRing_DTmanGradedAlgebra_Dat">
##  <ManSection>
##  <Attr Name="CoefficientsRing" Arg="A"/>
##
##  <Returns>
##    Ring
##  </Returns>
##  <Description>
##  Returns the ring of coefficients of the graded algebra presentation 
##  <A>A</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(CoefficientsRing,
  "of graded algebra presentation",
  [IsGradedAlgebraPresentation],
  function(A)
    return CoefficientsRing(A!.ring);
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="IndeterminatesOfGradedAlgebraPresentation_DTmanGradedAlgebra_Dat">
##  <ManSection>
##  <Attr Name="IndeterminatesOfGradedAlgebraPresentation" Arg="A"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  Returns the indeterminates used in the graded algebra presentation <A>A</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(IndeterminatesOfGradedAlgebraPresentation,
  [IsGradedAlgebraPresentation],
  function(A)
    return IndeterminatesOfPolynomialRing(A!.ring);
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="GeneratorsOfPresentationIdeal_DTmanGradedAlgebra_Dat">
##  <ManSection>
##  <Attr Name="GeneratorsOfPresentationIdeal" Arg="A"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  Returns the relations in the ring presentation for the graded algebra
##  <A>A</A>. The relations are returned sorted in order of increasing degree, 
##  and by indeterminate within each degree.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(GeneratorsOfPresentationIdeal,
  "for graded algebra presentation",
  [IsGradedAlgebraPresentation],
  function(A)
    local sortedrels, degs;
    sortedrels := ShallowCopy(A!.relations);
    # sort the relations first of all so that they are in indeterminate order
    Sort(sortedrels);
    # now sort by degree
    degs := List(A!.relations, p->DegreeOfRepresentative(A, p));
    SortParallel(degs, sortedrels);
    return sortedrels;
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="PresentationIdeal_DTmanGradedAlgebra_Dat">
##  <ManSection>
##  <Attr Name="PresentationIdeal" Arg="A"/>
##
##  <Returns>
##    Ideal
##  </Returns>
##  <Description>
##  Returns the ideal in the graded algebra presentation <A>A</A> as a &GAP;
##  ideal <Ref Sect="Ideal" BookName="ref"/>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(PresentationIdeal,
  "of graded algebra presentation",
  [IsGradedAlgebraPresentation],
  function(A)
    return Ideal(A!.ring, A!.relations);
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="IndeterminateDegrees_DTmanGradedAlgebra_Dat">
##  <ManSection>
##  <Attr Name="IndeterminateDegrees" Arg="A"/>
##
##  <Returns>
##  List
##  </Returns>
##  <Description>
##  Returns the degrees of the polynomial ring indeterminates in the graded
##  algebra presentation <A>A</A>. The ordering corresponds to the
##  order of the ring indeterminates returned by
##  <Ref Attr="IndeterminatesOfGradedAlgebraPresentation"/>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(IndeterminateDegrees,
  "of graded algebra presentation",
  [IsGradedAlgebraPresentation],
  function(A)
    return A!.degrees;
  end
);
#####################################################################


  
  
#####################################################################
##  <#GAPDoc Label="ViewObj_DTmanGradedAlgebraNODOC">
##  <ManSection>
##  <Meth Name="ViewObj" Arg="A" Label="for GradedAlgebraPresentation"/>
##
##  <Description>
##  Prints a description of the graded algebra presentation <A>A</A>. This is 
##  the usual description printed by &GAP;.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  ViewObj,
  "for GradedAlgebraPresentation",
  [IsGradedAlgebraPresentation],
  function(A)
    Print("Graded algebra ", CoefficientsRing(A), 
    IndeterminatesOfGradedAlgebraPresentation(A), " / ", 
    GeneratorsOfPresentationIdeal(A), " with indeterminate degrees ", 
    IndeterminateDegrees(A));
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="PrintObj_DTmanGradedAlgebraNODOC">
##  <ManSection>
##  <Meth Name="PrintObj" Arg="A" Label="for GradedAlgebraPresentation"/>
##
##  <Description>
##  Prints a detailed description of the graded algebra <A>A</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  PrintObj,
  "for GradedAlgebraPresentation",
  [IsGradedAlgebraPresentation],
  function(A)
    Print("Graded algebra ", BaseRing(A), " / ", 
    GeneratorsOfPresentationIdeal(A), " with indeterminate degrees ", 
    IndeterminateDegrees(A));
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="Display_DTmanGradedAlgebraNODOC">
##  <ManSection>
##  <Meth Name="Display" Arg="A" Label="for GradedAlgebraPresentation"/>
##
##  <Description>
##  Displays the graded algebra <A>A</A> in a human-readable form.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  Display,
  "for GradedAlgebraPresentation",
  [IsGradedAlgebraPresentation],
  function(A)
    ViewObj(A);
  end
  );
#####################################################################



#####################################################################
##  <#GAPDoc Label="Equals_DTmanGradedAlgebraNODOC">
##  <ManSection>
##  <Oper Name="&#x003D;" Arg="A, B" Label="for GradedAlgebraPresentation"/>
##
##  <Returns>
##  Boolean
##  </Returns>
##  <Description>
##  Returns <K>true</K> if the graded algebra presentations <A>A</A> and 
##  <A>B</A> are the same.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(\=,
  "for GradedAlgebraPresentation",
  IsIdenticalObj,
  [IsGradedAlgebraPresentation, IsGradedAlgebraPresentation],
  function(A, B)
    return BaseRing(A) = BaseRing(B) and 
      GeneratorsOfPresentationIdeal(A) = GeneratorsOfPresentationIdeal(B) and
      IndeterminateDegrees(A) = IndeterminateDegrees(B);
end);
#####################################################################




#####################################################################
##  <#GAPDoc Label="TensorProduct_DTmanGradedAlgebra_Func">
##  <ManSection>
##  <Heading>TensorProductOp</Heading>
##  <Oper Name="TensorProductOp" Arg="A, B" Label="for two algebra presentations"/>
##  <Oper Name="TensorProductOp" Arg="coll" Label="for collection of algebra presentations"/>
##
##  <Returns>
##  GradedAlgebraPresentation
##  </Returns>
##  <Description>
##  Returns a presentation for the graded algebra that is the tensor product of 
##  two graded algebras presented by <A>A</A> and <A>B</A>, or of a list
##  of graded algebras.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallOtherMethod(TensorProductOp,
  "for list of GradedAlgebraPresentations",
  [IsHomogeneousList],
  function(coll)
    local A, i;

    if not IsGradedAlgebraPresentation(coll[1]) then
      TryNextMethod();
    fi;

    A := coll[1];
    for i in [2..Length(coll)] do
      A := TensorProductOp(A, coll[i]);
    od;
    return A;
  end
);
#####################################################################
InstallMethod(TensorProductOp,
  "for GradedAlgebraPresentation",
  IsIdenticalObj,
  [IsGradedAlgebraPresentation, IsGradedAlgebraPresentation],
  function(A, B)
    local Aindets, Bindets, outR, tempR, A2, B2, merged, degrees, indets;
    
    if CoefficientsRing(A) <> CoefficientsRing(A) then
      Error("<A> and <B> must have the same coefficients ring");
    fi;
    # The tensor product will have the combined generators
    # and relations of A and B. We want to make sure that 
    # A, B and the output ring all have distinct indeterminates, so
    # do some shuffling around
    Aindets := IndeterminatesOfGradedAlgebraPresentation(A);
    Bindets := IndeterminatesOfGradedAlgebraPresentation(B);
    outR := PolynomialRing(CoefficientsRing(A), 
      Length(Aindets) + Length(Bindets));
    
    # Express A in new indeterminates that don't clash with outR
    A2 := HAPPRIME_GradedAlgebraPresentationAvoidingIndeterminates(
      A, IndeterminatesOfPolynomialRing(outR));
    Aindets := IndeterminatesOfGradedAlgebraPresentation(A2);

    # Express B in new indeterminates that don't clash with outR or A2
    B2 := HAPPRIME_GradedAlgebraPresentationAvoidingIndeterminates(
      B, Concatenation(IndeterminatesOfPolynomialRing(outR), Aindets));
    Bindets := IndeterminatesOfGradedAlgebraPresentation(B2);
    
    # Merge A2 and B2 into a new graded algebra presentation
    merged := GradedAlgebraPresentationNC(
      PolynomialRing(CoefficientsRing(A), Concatenation(Aindets, Bindets)),
      Concatenation(GeneratorsOfPresentationIdeal(A2), GeneratorsOfPresentationIdeal(B2)),
      Concatenation(IndeterminateDegrees(A2), IndeterminateDegrees(B2)));
      
    # Now convert the indeterminates back to the outR ring
    # sorting the degrees into increasing order
    degrees := List(IndeterminatesOfGradedAlgebraPresentation(merged), 
      i->DegreeOfRepresentative(merged, i));

    indets := ShallowCopy(IndeterminatesOfGradedAlgebraPresentation(merged));
    SortParallel(degrees, indets);
      
    return GradedAlgebraPresentationNC(
      outR, 
      HAPPRIME_SwitchPolynomialIndeterminates(
        BaseRing(merged), outR, GeneratorsOfPresentationIdeal(merged)), 
        degrees);
  end
);
#####################################################################




#####################################################################
##  <#GAPDoc Label="AreIsomorphicGradedAlgebras_DTmanGradedAlgebra_Func">
##  <ManSection>
##  <Oper Name="AreIsomorphicGradedAlgebras" Arg="A, B"/>
##
##  <Description>
##  Returns <K>true</K> if the graded algebras <A>A</A> and <A>B</A> are
##  isomorphic, or <K>false</K> otherwise. This function tries all possible
##  ring isomorphisms, so may take a considerable length of time for 
##  graded algebras with a large number of dimensions in each degree.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(AreIsomorphicGradedAlgebras,
  [IsGradedAlgebraPresentation, IsGradedAlgebraPresentation],
  function(A, B)
    local a, b, Aindets, Bindets, B2, ord, Arels, Aleading, B2rels,
    maxdegree, generatordegrees, generatordegrees2, i, d, u, v, 
    totals1, totals2, m, n, total, degrees, vectorspacebasis, noInducedBasis, 
    switchAB, isomorphism, 
    Therm, Thermometer, DegreeOfPolynomial, FindIsomorphism;

    ##########################################
    Thermometer := function(totals)
      local n, space, cols, divs, spinset, spin, counts, bar, nononecols, 
      Update;

      n := Length(totals);
      # How much space do we have available for thermometer blocks?
      # |===   | |====  | |======| /
      # Each thermometer takes up two columns, plus one space between,
      # plus three spaces at the end for the spinner
      space := SizeScreen()[1] - 3*n - 4;
#      space := SizeScreen()[1] - 3*n - 24;
      if 2*n > space then
        Error("too many things for the thermometer");
      fi;

      if Sum(totals) < space then
        # If we just have small totals, we can probably do this exactly
        cols := totals;
      else
        # Divide this space proportiantly between the thermometers 
        cols := List(totals, i->QuoInt(space * i, Sum(totals)));
        # if any are less than one, make it one and reduce the largest to 
        # make sure that we don't exceed the space
        repeat
          if Minimum(cols) < 1 then
            cols[Position(cols, Minimum(cols))] := 1;
            if Sum(cols) > space then
              cols[Position(cols, Maximum(cols))] := 
              cols[Position(cols, Maximum(cols))] - (Sum(cols) - space);
            fi;
          fi;
        until Minimum(cols) >= 1;
        # We in fact want to make sure that any whose total isn't one is 
        # of size two
        nononecols := Filtered([1..Length(totals)], i->totals[i] > 1);
        repeat
          if Minimum(cols{nononecols}) < 2 then
            cols[Position(cols, Minimum(cols{nononecols}))] := 2;
            if Sum(cols) > space then
              cols[Position(cols, Maximum(cols))] := 
              cols[Position(cols, Maximum(cols))] - (Sum(cols) - space);
            fi;
          fi;
        until Minimum(cols{nononecols}) >= 2;
      fi;

      # What do we divide things by?
      divs := List([1..n], i->totals[i]/cols[i]);

      spinset := ["\\", "/", "-"];
      spin := 0;
      counts := ListWithIdenticalEntries(n, 0);

      if noInducedBasis then
        bar := "=";
      else
        bar := ">";
      fi;

 #     Print("Totals are ", totals, "\n");
#      Print("Cols are ", cols, "\n");
#      Print("Divs are ", divs, "\n");


      ######################################
      ## The function that does an update and displays
      Update := function(i, val)
        local curs, j;

        counts[i] := val;
        if InfoLevel(InfoHAPprime) > 0 then
          spin := (spin + 1) mod 3;

          curs := List([1..n], 
            i->Minimum(Int(counts[i]/divs[i] + 1/2), cols[i]));
          # Make sure that its only zero if counts is zero,
          # and that it is only cols if counts is totals
          for i in [1..Length(curs)] do
            if curs[i] = 0 and counts[i] <> 0 then
              curs[i] := 1;
            fi;
            if curs[i] = cols[i] and counts[i] <> totals[i] then
              curs[i] := curs[i] - 1;
            fi;
          od;
          for i in [1..n] do
            Print("|");
            for j in [1..curs[i]] do
              Print(bar);
            od;
            for j in [(curs[i]+1)..cols[i]] do
              Print(".");
            od;
            Print("| ");
          od;
          Print(spinset[spin+1], "\r");
#          Print(spinset[spin+1], " ", counts, "\n");
        fi;
      end;
      ######################################

      return rec(
        Update := Update, 
        totals := totals,
        cols := cols,
        divs := divs,
        counts := counts);
    end;
    ##########################################

    # Check that the coefficients of the ring are the same
    if CoefficientsRing(A) <> CoefficientsRing(B) then
      return false;
    fi;

    # Check that the number of indeterminates at least agrees
    if Length(IndeterminateDegrees(A)) <> Length(IndeterminateDegrees(B)) then
      return false;
    fi;

    # Check that the degrees agree
    a := ShallowCopy(IndeterminateDegrees(A)); Sort(a);
    b := ShallowCopy(IndeterminateDegrees(B)); Sort(b);
    if a <> b then
      return false;
    fi;

    # See if the relations are the same
    if Set(GeneratorsOfPresentationIdeal(A)) = Set(GeneratorsOfPresentationIdeal(B)) then
      return true;
    fi;

    # Check the Hilbert series
    a := HilbertPoincareSeries(A);
    b := HilbertPoincareSeries(B);
    if a <> b then
      return false;
    fi;

    # OK, so we'll have to search for an isomorphism 
    # Choose whichever of our A or B has the most relations as our 'A'
    switchAB := false;
    if Length(GeneratorsOfPresentationIdeal(B)) > Length(GeneratorsOfPresentationIdeal(A)) then
      switchAB := ShallowCopy(A);
      A := ShallowCopy(B);
      B := switchAB;
    fi;
    

    # Convert B to a new set of indeterminates so that we can then map it 
    # back to something else
    Aindets := IndeterminatesOfGradedAlgebraPresentation(A);
    Bindets := IndeterminatesOfGradedAlgebraPresentation(B);
    B2 := HAPPRIME_GradedAlgebraPresentationAvoidingIndeterminates(
      B, Set(Concatenation(Aindets, Bindets)));

    # Make sure that we have Groebner Bases for A's relations
    ord := MonomialLexOrdering();
    Arels := Set(HAPPRIME_SingularReducedGroebnerBasis(
      GeneratorsOfPresentationIdeal(A), ord));
    # And also get the leading monomials from Arels
    Aleading := List(Arels, i->LeadingMonomialOfPolynomial(i, ord));

    # Make a list of the possible degrees of the indeterminates
    degrees := Set(IndeterminateDegrees(A));
    
    # What is the maximum degree of a relation in B?
    maxdegree := Maximum(List(
      GeneratorsOfPresentationIdeal(B), i->DegreeOfRepresentative(B, i)));
    # We can ignore generators in higher degrees
    degrees := Filtered(degrees, i->i <= maxdegree);
    # Make a list of the generators of each degree
    # in A and in B2
    generatordegrees := List([1..maxdegree], i->[]);
    generatordegrees2 := List([1..maxdegree], i->[]);
    for i in [1..Length(Aindets)] do
      # We can ignore generators with degree higher than any relation
      # since they have no relations!
      if IndeterminateDegrees(A)[i] <= maxdegree then
        Add(generatordegrees[IndeterminateDegrees(A)[i]], Aindets[i]);
      fi;
      if IndeterminateDegrees(B2)[i] <= maxdegree then
        Add(generatordegrees2[IndeterminateDegrees(B2)[i]], 
          IndeterminatesOfGradedAlgebraPresentation(B2)[i]);
      fi;
    od;

    # Form the basis for the whole vector space of the algebra A 
    # up to maxdegree
    vectorspacebasis := List([1..maxdegree], i->[]);
    for d in [1..maxdegree] do
      for i in [1..(d-1)] do
        # Add in the product of every basis element in degree i with every 
        # basis element in degree (d-i)
        for u in vectorspacebasis[i] do
          for v in vectorspacebasis[d-i] do
            AddSet(vectorspacebasis[d], u*v);
          od;
        od;
      od;
      # Reduce everything using the ideal to get it into its simplest form
      vectorspacebasis[d] := List(vectorspacebasis[d],
        i->PolynomialReducedRemainder(i, Arels, ord)); 
      # Remove anything that is zero
      vectorspacebasis[d] := Filtered(vectorspacebasis[d], i->not IsZero(i)); 
      # And that are not monomials
      vectorspacebasis[d] := Filtered(vectorspacebasis[d], IsMonomial); 
      # Remove duplicates
      vectorspacebasis[d] := Set(vectorspacebasis[d]);
      # And add any new generators in this degree
      Append(vectorspacebasis[d], generatordegrees[d]);
    od; 

    # As a sanity check, we'll make sure that this agrees with the Hilbert Series
    if List(vectorspacebasis, Length) <> 
      CoefficientsOfPoincareSeries(A, maxdegree+1){[2..maxdegree+1]} then
        Error("something has gone wrong computing the vector space basis");
    fi;

    ## Calculate the maximum combinations in each degree
    totals1 := [];
    totals2 := [];
    for d in degrees do
      m := Length(Filtered(
        vectorspacebasis[d], v->not v in generatordegrees[d]));
      n := Length(generatordegrees[d]);
      total := Length(EnumeratorByBasis(CanonicalBasis(CoefficientsRing(A)^m)));
      total := total^n * Order(GL(n, Size(CoefficientsRing(A))));
      Add(totals2, total);
      Add(totals1, Order(GL(n, Size(CoefficientsRing(A)))));
    od;

    #############################################
    # Recursive function to find an isomorphism if one exists
    # this enmumerates all isomorphisms in degree given by degrees[<di>] and 
    # calls FindIsomorphism for each of these in turn. <partialmap> contains
    # the images under the proposed isomorphism of all the generators of degree
    # less than <degree>
    FindIsomorphism := function(di, partialmap)
      
      local degree, nextdegree, Arelsdeg, Brelsdeg, B2relsdeg, inducedbasis,
      m, n, V, Venumerators, M, Mhom, Genumerator, mapping, count,
      generatormap, i, newmap, message;
      
      # What is this degree, and what is the next degree?
      degree := degrees[di];
      if di < Length(degrees) then
        nextdegree := degrees[di+1];
      else
        nextdegree := false;
      fi;
      # Filter the ideals Arels and Brels so that we only have relations
      # up to the ones that will be considered by the next call to 
      # FindIsomorphism, i.e. up to degree (next degree - 1)
      if nextdegree <> false then
        Arelsdeg := Filtered(Arels, i->DegreeOfRepresentative(A, i) < nextdegree);
        B2relsdeg := Filtered(GeneratorsOfPresentationIdeal(B2), i->DegreeOfRepresentative(B2, i) < nextdegree);
      else;
        Arelsdeg := Arels;
        B2relsdeg := GeneratorsOfPresentationIdeal(B2);
      fi;

      # Get the vector-space basis elements of this degree that 
      # are not new generators, i.e. the ones that are generated by
      # generators of smaller degree
      if noInducedBasis then
        inducedbasis := [];
      else
        inducedbasis := Filtered(
          vectorspacebasis[degree], v->not v in generatordegrees[degree]);
      fi;

      # How many induced basis elements are there?
      m := Length(inducedbasis);
      # How many new generators do I have?
      n := Length(generatordegrees[degree]);

      # Create enumerators to enumerate possible combinations of the induced 
      # basis. We want one of these for each new generator
      V := CoefficientsRing(A)^m;
      Venumerators := List([1..n], i->EnumeratorByBasis(CanonicalBasis(V)));

      # Also enumerate the linear mappings of my new generators
      M := GL(n, Size(CoefficientsRing(A)));
      # Because GAP can't enumerate GL(n,q) nicely, let's instead convert this 
      # to a permutation group and enumerate over this instead
      Mhom := IsomorphismPermGroup(M);
      M := Image(Mhom, M);
      Genumerator := Enumerator(M);

      # And start off with the first entry in each list
      # The first entry in mapping corresponds to the Genumerator, and
      # the rest to the elements of Venumerators
      mapping := ListWithIdenticalEntries(Length(Venumerators) + 1, 1);
      # How many possible mappings are there?
      count := 1;
      repeat
        Therm.Update(di, count);

        # Map the new generators under this mapping
        # First map the generators under the linear mapping
        generatormap := generatordegrees[degree]*
          PreImage(Mhom, Genumerator[mapping[1]]);

        # Now add on some combination of the induced basis
        if not IsEmpty(inducedbasis) then
          for i in [1..n] do
            generatormap[i] := generatormap[i] + 
              Venumerators[i][mapping[i+1]]*inducedbasis;
          od;
        fi;

        # So we can now augment our partialmap with the new mapping
        newmap := [
          Concatenation(partialmap[1], generatordegrees2[degree]),
          Concatenation(partialmap[2], generatormap)];

        # Convert our Brelsdeg with this mapping and see if its Groebner basis
        # matches Arelsdeg
        Brelsdeg := Set(HAPPRIME_SingularReducedGroebnerBasis(
          HAPPRIME_MapPolynomialIndeterminates(
            newmap[1], newmap[2], B2relsdeg), ord));
        # The Groebner Basis might add some more relations with a degree
        # that we're not considering yet. Filter those out 
        Brelsdeg := Filtered(Brelsdeg, i->DegreeOfRepresentative(B, i) < nextdegree);

        # If these are the same then this mapping works (up to the current 
        # degree)
        if Arelsdeg = Brelsdeg then
          if nextdegree <> false then      
            # See if we can find a mapping that works in a higher degree
            if FindIsomorphism(di+1, newmap) then
              Therm.Update(di, 0);
              return true;
            fi;
          else
            # We have found a mapping that works. 
            # The rings are therefore isomorphic. Hurrah!
            if InfoLevel(InfoHAPprime) >= 1 then
              # Convert the indeterminants in newmap[1] back into
              # the original indeterminants
              newmap[1] := HAPPRIME_SwitchPolynomialIndeterminates(
                BaseRing(B2), BaseRing(B), newmap[1]);
              if switchAB = false then
                Info(InfoHAPprime, 1, "isomorphism (B -> A):");
              else
                Info(InfoHAPprime, 1, "isomorphism (A -> B):");
              fi;
              for i in [1..Length(newmap[1])] do
                Info(InfoHAPprime, 1, "  ", String(newmap[1][i]), " -> ", String(newmap[2][i]));
              od;
            fi;
              Therm.Update(di, 0);
            return true;
          fi;
        fi;

        # Step to the next mapping. We enumerate through all the new generator
        # mappings first, and then each of the vector space combinations
        mapping[1] := mapping[1] + 1;
        if mapping[1] > Length(Genumerator) then
          mapping[1] := 1;
          i := 1;
          repeat
            mapping[i+1] := mapping[i+1] + 1;
            if mapping[i+1] > Length(Venumerators[i]) then
              mapping[i+1] := 1;
              i := i + 1;
              if i > Length(Venumerators) then
                mapping := false;
              fi;
            else
              break; 
            fi;
          until mapping = false;
        fi;

        count := count + 1;
      until mapping = false;

      Therm.Update(di, 0);
      return false;
    end;
    #############################################


    
    # First see if we can find a ismorphism by just mapping the generators
    noInducedBasis := true;
    Therm := Thermometer(totals1);
    if FindIsomorphism(1, [[],[]]) = true then
      isomorphism := true;
    else
      noInducedBasis := false;
      # Now try with a full mapping of the vector space basis
      Therm := Thermometer(totals2);
      if FindIsomorphism(1, [[],[]]) = true then
        isomorphism := true;
      else
        isomorphism := false;
      fi;
    fi;
    
    # Switch A and B back if they were switched
    if switchAB <> false then
      A := ShallowCopy(B);
      A := switchAB;
    fi;
    
    return isomorphism;

  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="IsAssociatedGradedRing_DTmanGradedAlgebra_Func">
##  <ManSection>
##  <Oper Name="IsAssociatedGradedRing" Arg="A, B"/>
##
##  <Description>
##  Returns <K>true</K> if the algebra <A>A</A> is an associated
##  graded ring of the algebra <A>B</A>. This is the case if the additive structure is the
##  same (i.e. the Hilbert-Poincaré series is the same), and the generators
##  for <A>A</A> (and their degrees) are included in the generators for
##  <A>B</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(IsAssociatedGradedRing,
  "for graded algebra presentation",
  [IsGradedAlgebraPresentation, IsGradedAlgebraPresentation],
  function(A, B)

    local Agens, Bgens, d, p, a, b;
    # Check that the generator degrees for A are a subset of those for B
    Agens := ShallowCopy(IndeterminateDegrees(A)); Sort(Agens);
    Bgens := ShallowCopy(IndeterminateDegrees(B)); Sort(Bgens);
    for d in Agens do
      p := Position(Bgens, d);
      if p = fail then
        return false;
      fi;
      Remove(Bgens, p);
    od;

    # Check the Hilbert series
    a := HilbertPoincareSeries(A);
    b := HilbertPoincareSeries(B);
    if a <> b then
      return false;
    fi;

    return true;

  end
);
#####################################################################



#####################################################################
##  <#GAPDoc Label="DegreeOfRepresentative_DTmanGradedAlgebra_Func">
##  <ManSection>
##  <Oper Name="DegreeOfRepresentative" Arg="A, p"/>
##
##  <Returns>
##    Integer
##  </Returns>
##  <Description>
##  Returns the degree of a polynomial representative <A>p</A> from the graded 
##  ring presentation <A>A</A>. 
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(DegreeOfRepresentative,
  "for graded algebra presentation",
  [IsGradedAlgebraPresentation, IsPolynomial],
  function(A, p)

  local degree, terms, t, unimons, deg, i, indetExp, k;

  degree := 0;
  terms := TermsOfPolynomial(p);
  for t in terms do
    unimons := UnivariateMonomialsOfMonomial(t[1]);
    deg := 0;
    for i in unimons do
      indetExp := IndeterminateAndExponentOfUnivariateMonomial(i);
      k := Position(IndeterminatesOfGradedAlgebraPresentation(A), indetExp[1]);
      if k = fail then
        Error("<p> is not a representative in the algebra <A>");
      fi;
      deg := deg + IndeterminateDegrees(A)[k] * indetExp[2];
    od;
    if degree = 0 then
      degree := deg;
    else
      if degree <> deg then
        Error("<p> is not a homogeneous polynomial in the algebra <A>");
      fi;
    fi;
  od;

  return deg;
end);
#########################################


#####################################################################
##  <#GAPDoc Label="MaximumDegreeForPresentation_DTmanGradedAlgebra_Func">
##  <ManSection>
##  <Attr Name="MaximumDegreeForPresentation" Arg="A"/>
##
##  <Returns>
##    Integer
##  </Returns>
##  <Description>
##  Returns the maximum degree in generators or relations that is needed to 
##  generate the graded algebra presentation <A>A</A>. 
##  This is not necessarily the same as the 
##  largest degree in any of the relations and generators - some relations
##  may be redundant (for example due to being a Groebner basis), so 
##  this routine checks for the largest degree of a required generator, and
##  returns the maximum of this and the generator degrees.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(MaximumDegreeForPresentation,
  "for graded algebra presenatation",
  [IsGradedAlgebraPresentation],
  function(A)

    local sortedrels, degs, PS, maxdeg, maxgen;
    
    # If there are no relations then it is just the largest generator degree
    if IsEmpty(GeneratorsOfPresentationIdeal(A)) then 
      return Maximum(IndeterminateDegrees(A));
    fi;
    
    # Sort the relations by degree
    sortedrels := ShallowCopy(GeneratorsOfPresentationIdeal(A));
    degs := List(GeneratorsOfPresentationIdeal(A), p->DegreeOfRepresentative(A, p));
    SortParallel(degs, sortedrels);

    maxgen := Maximum(IndeterminateDegrees(A));
    PS := HilbertPoincareSeries(A);
    repeat
      # Remove the largest degree and see if the Hilbert Series is the same
      maxdeg := degs[Length(degs)];
      repeat
        Remove(degs, Length(degs));
        Remove(sortedrels, Length(sortedrels));
      until IsEmpty(degs) or degs[Length(degs)] <> maxdeg;

      if HilbertPoincareSeries(
        GradedAlgebraPresentation(
          BaseRing(A), sortedrels, IndeterminateDegrees(A))) <> PS then
        # We need the relations of maxdeg. Return the larger of
        # maxdeg and the generator degrees
        return Maximum(maxdeg, maxgen);
      fi;
    until IsEmpty(sortedrels) or maxdeg <= maxgen;

    return maxgen;

  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="SubspaceDimensionDegree_DTmanGradedAlgebra_Func">
##  <ManSection>
##  <Heading>SubspaceDimensionDegree</Heading>
##  <Oper Name="SubspaceDimensionDegree" Arg="A, d" Label="for one degree"/>
##  <Oper Name="SubspaceDimensionDegree" Arg="A, degs" Label="for list of degrees"/>
##
##  <Returns>
##    Integer or list
##  </Returns>
##  <Description>
##  Returns the dimension of degree <A>d</A> of the graded algebra <A>A</A>, 
##  or a list of dimensions corresponding to the list of degrees <A>degs</A>. 
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(SubspaceDimensionDegree,
  "for graded algebra presentation and one degree",
  [IsGradedAlgebraPresentation, IsPosInt],
  function(A, d)
    return SubspaceDimensionDegree(A, [d])[1];
  end
);
#####################################################################
InstallOtherMethod(SubspaceDimensionDegree,
  "for graded algebra presentation and list of degrees",
  [IsGradedAlgebraPresentation, IsHomogeneousList],
  function(A, degs)
    
    return CoefficientsOfPoincareSeries(A, Maximum(degs+1)){degs+1};
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="SubspaceBasisRepsByDegree_DTmanGradedAlgebra_Func">
##  <ManSection>
##  <Heading>SubspaceBasisRepsByDegree</Heading>
##  <Oper Name="SubspaceBasisRepsByDegree" Arg="A, d" Label="for one degree"/>
##  <Oper Name="SubspaceBasisRepsByDegree" Arg="A, degs" Label="for list of degrees"/>
##
##  <Returns>
##    List or list of lists
##  </Returns>
##  <Description>
##  Returns a basis for degree <A>d</A> of the graded algebra <A>A</A>, 
##  or a list of bases for the list of degrees <A>degs</A>. Each basis
##  is returned as a list of representatives.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(SubspaceBasisRepsByDegree,
  "for one degree",
  [IsGradedAlgebraPresentation, IsPosInt],
  function(A, d)
    return SubspaceBasisRepsByDegree(A, [d])[1];
  end
);
#####################################################################
InstallOtherMethod(SubspaceBasisRepsByDegree,
  "for list of degrees",
  [IsGradedAlgebraPresentation, IsHomogeneousList],
  function(A, degs)
    local maxdegree, vectorspacebasis, d, i, u, v, Arels, ord, generatordegrees, 
      Aindets;

    maxdegree := Maximum(degs);
    ord := MonomialLexOrdering();
    Arels := Set(HAPPRIME_SingularReducedGroebnerBasis(
      GeneratorsOfPresentationIdeal(A), ord));
    generatordegrees := List([1..maxdegree], i->[]);
    Aindets := IndeterminatesOfGradedAlgebraPresentation(A);
    for i in [1..Length(Aindets)] do
      # We can ignore generators with degree higher than any relation
      # since they have no relations!
      if IndeterminateDegrees(A)[i] <= maxdegree then
        Add(generatordegrees[IndeterminateDegrees(A)[i]], Aindets[i]);
      fi;
    od;
  
    # Form the basis for the whole vector space of the algebra A 
    # up to maxdegree
    vectorspacebasis := List([1..maxdegree], i->[]);
    for d in [1..maxdegree] do
      for i in [1..(d-1)] do
        # Add in the product of every basis element in degree i with every 
        # basis element in degree (d-i)
        for u in vectorspacebasis[i] do
          for v in vectorspacebasis[d-i] do
            AddSet(vectorspacebasis[d], u*v);
          od;
        od;
      od;
      # Reduce everything using the ideal to get it into its simplest form
      vectorspacebasis[d] := List(vectorspacebasis[d],
        i->PolynomialReducedRemainder(i, Arels, ord)); 
      # Remove anything that is zero
      vectorspacebasis[d] := Filtered(vectorspacebasis[d], i->not IsZero(i)); 
      # And that are not monomials
      vectorspacebasis[d] := Filtered(vectorspacebasis[d], IsMonomial); 
      # Remove duplicates
      vectorspacebasis[d] := Set(vectorspacebasis[d]);
      # And add any new generators in this degree
      Append(vectorspacebasis[d], generatordegrees[d]);
    od; 

    # As a sanity check, we'll make sure that this agrees with the Hilbert Series
    if List(vectorspacebasis, Length) <> 
      CoefficientsOfPoincareSeries(A, maxdegree+1){[2..maxdegree+1]} then
        Error("something has gone wrong computing the vector space basis");
    fi;

    return vectorspacebasis{degs};
  end
);
#####################################################################



#if LoadPackage("singular") = true then
if IsPackageMarkedForLoading("singular","0") then
  #####################################################################
  ##  <#GAPDoc Label="CoefficientsOfPoincareSeries_DTmanGradedAlgebra_Func">
  ##  <ManSection>
  ##  <Oper Name="CoefficientsOfPoincareSeries" Arg="A, n"/>
  ##
  ##  <Returns>
  ##    List
  ##  </Returns>
  ##  <Description>
  ##  Returns the first <A>n</A> coefficients of the Poincaré series for the 
  ##  graded algebra with <A>A</A>. These are equal to the dimensions of degrees 
  ##  0 to <M>n-1</M> of the algebra (a fact that is used in the function
  ##  <Ref Func="SubspaceDimensionDegree" Label="for one degree"/>).
  ##  <P/>
  ##  <E>This function uses the <Package>singular</Package> package.</E>
  ##  </Description>
  ##  </ManSection>
  ##  <#/GAPDoc>
  #####################################################################
  InstallMethod(CoefficientsOfPoincareSeries, 
    [IsGradedAlgebraPresentation, IsPosInt],
    function(A, n)
  
      local series, T, t, poly, den, d;

      # Get the Hilbert-Poincare series
      series := HAPPRIME_HilbertSeries(A);

      # Turn this into a polynomial
      T := PolynomialRing(Integers, 1);
      t := IndeterminatesOfPolynomialRing(T)[1];
      series := series * List([0..Length(series)-1], i->t^i);
    
      # The Hilbert-Poincare Series is Q(t)/(1-t)^n
      # We have Q(t) - this is what we call the Hilbert Series - now we
      # want to generate (1-t)^n. When we have degrees \neq 1, this is 
      # in fact \prod (1-t^d)
      # Calculate this up to degree n using the binomial theorem
      den := One(t);
      for d in IndeterminateDegrees(A) do
        poly := Sum(List([0..n], i->t^(d*i)));
        den := den*poly;
      od;

      # The product of these two is the Hilbert-Poincare Series (accurate to
      # degree n)
      series := series * den;

      # And return the first n terms
      return CoefficientsOfUnivariatePolynomial(series){[1..n]};
    end
  );
  #####################################################################
    
  
  #####################################################################
  ##  <#GAPDoc Label="HilbertPoincareSeries_DTmanGradedAlgebra_Func">
  ##  <ManSection>
  ##  <Attr Name="HilbertPoincareSeries" Arg="A"/>
  ##
  ##  <Returns>
  ##    Rational function
  ##  </Returns>
  ##  <Description>
  ##  Returns the Poincaré series for the graded algebra <A>A</A>. 
  ##  This is a rational function <M>P(t)/Q(t)</M> which is a 
  ##  is a polynomial whose coefficients are the dimensions of each degree 
  ##  of the algebra. 
  ##  <P/>
  ##  <E>This function uses the <Package>singular</Package> package.</E>
  ##  </Description>
  ##  </ManSection>
  ##  <#/GAPDoc>
  #####################################################################
  # EVIL! This should be called just PoincareSeries for consistency, but that
  # name is taken for a HAP global function, so I can't use that.
  InstallMethod(HilbertPoincareSeries,
    [IsGradedAlgebraPresentation],
    function(A)
      local H, T, t, series, den, d, coeffs, fun, gcd;

      series := HAPPRIME_HilbertSeries(A);

      # Turn this into a polynomial
      T := PolynomialRing(Integers, 1);
      t := IndeterminatesOfPolynomialRing(T)[1];
      series := series * List([0..Length(series)-1], i->t^i);

      # The denominator is (1-t)^n, or when we have degrees \neq 1, this is
      # in fact \prod (1-t^d)
      den := One(t);
      for d in IndeterminateDegrees(A) do
        den := den*(1-t^d);
      od;

      # Turn this into a rational function so that we can reduce it
      fun := series/den;
      # This function guarantees to cancel the polynomials
      coeffs := List(CoefficientsOfUnivariateRationalFunction(fun), 
        x->ShallowCopy(x));
      # We shall make sure that the leading coefficients are in their
      # smallest form, and that the leading coefficient in the numerator
      # is positive.
      for d in [1..Length(coeffs)] do                 # Graham added this 
      if IsList(coeffs[d]) then                       # Oct 2023
      if Length(coeffs[d])=0 then coeffs[d]:=[0]; fi; #
      fi;                                             #
      od;
      gcd := Gcd(coeffs[1][1], coeffs[2][1]) * SignInt(coeffs[1][1]);  
      coeffs[1] := coeffs[1] / gcd;
      coeffs[2] := coeffs[2] / gcd;
      # Turn these back into polynomials
      series := coeffs[1]* List([0..Length(coeffs[1])-1], i->t^i) * t^coeffs[3];
      den := coeffs[2]* List([0..Length(coeffs[2])-1], i->t^i);

      return series/den;
    end
  );
  #####################################################################
    
  
  #####################################################################
  ##  <#GAPDoc Label="HAPPRIME_SingularHilbertSeries_DTmanGradedAlgebraInt">
  ##  <ManSection>
  ##  <Attr Name="HAPPRIME_HilbertSeries" Arg="P"/>
  ##
  ##  <Returns>
  ##    List
  ##  </Returns>
  ##  <Description>
  ##  The Hilbert-Poincaré series for a graded ring is a polynomial 
  ##  <M>H_P(t)</M> whose coefficients are the dimensions of each degree of the 
  ##  ring. It can be written as 
  ##  <Alt Only="LaTeX">
  ##    <Display>
  ##      H_P(t) = \frac{Q(t)}{\prod_{i=1}^n (1-t^{d_i})}
  ##    </Display>
  ##  </Alt>
  ##  <Alt Not="LaTeX">
  ##    <Display>
  ##      H_P(t) = Q(t) / (1-t^a)(1-t^b)...(1-t^n)
  ##    </Display>
  ##  </Alt>
  ##  This function returns the list of coefficients for the polynomial
  ##  <M>Q(t)</M> in the Hilbert Series for the graded algebra with 
  ##  presentation <A>P</A>. 
  ##  <P/>
  ##  This function simply calls the Singular function <C>hilb</C> via the 
  ##  &singular; package and returns the result.
  ##  </Description>
  ##  </ManSection>
  ##  <#/GAPDoc>
  #####################################################################
  InstallMethod(HAPPRIME_HilbertSeries, 
    [IsGradedAlgebraPresentation],
    function(A)
  
      local R, I, degs, series, Ig;
    
      R := BaseRing(A);
      I := GeneratorsOfPresentationIdeal(A);
      degs := IndeterminateDegrees(A);
  
      if Length(degs) <> Length(
        IndeterminatesOfPolynomialRing(R)) then
          Error("the list of degrees must be the same length as the number of ring indeterminates");
      fi;
    
      if IsEmpty(I) then 
        # if the ideal is empty then the Hilbert Series is just one.
        return [1];
      else
        # Use Singular to calculate the Hilbert Series
        SetTermOrdering(R, "dp");
        SingularSetBaseRing(R);
        # Make sure that I is a Groebner Basis with respect to the ordering that
        # Singular is using
        Ig := Ideal(R, GroebnerBasis(Ideal(R,I)));
        # Now get the Hilbert Series
        series := SingularInterface("hilb", [Ig, 1, degs], "intvec");
    
        # Remove the last one element of the list since that is not part of
        # the series (see the Singular help)
        Remove(series, Length(series));
    
        return series;
      fi;
    end
  );
  #####################################################################
else
  InstallMethod(CoefficientsOfPoincareSeries, 
    [IsGradedAlgebraPresentation, IsPosInt],
    function(A, n)
      Error("The package 'singular' cannot be loaded, so this 'HAPprime' function is not available");
    end
  );
  InstallMethod(HilbertPoincareSeries, 
    "for GradedAlgebraPresentation",
    [IsGradedAlgebraPresentation],
    function(A)
      Error("The package 'singular' cannot be loaded, so this 'HAPprime' function is not available");
    end
  );
  InstallGlobalFunction(HAPPRIME_HilbertSeries, 
    function(A)
      Error("The package 'singular' cannot be loaded, so this 'HAPprime' function is not available");
    end
  );
fi;

#####################################################################
##  <#GAPDoc Label="HAPPRIME_SwitchGradedAlgebraRing_DTmanGradedAlgebraInt">
##  <ManSection>
##  <Attr Name="HAPPRIME_SwitchGradedAlgebraRing" Arg="A, R"/>
##
##  <Returns>
##    GradedAlgebraPresentation
##  </Returns>
##  <Description>
##  Returns a new presentation for the graded algebra <A>A</A> which uses the 
##  ring <A>R</A> instead of the one in <A>A</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_SwitchGradedAlgebraRing,
  function(A, R)
    
    if CoefficientsRing(R) <> CoefficientsRing(A) then
      Error("<A> and <R> must have the same coefficients ring");
    fi;
    if Length(IndeterminatesOfPolynomialRing(R)) <>
      Length(IndeterminatesOfGradedAlgebraPresentation(A)) then
      Error("<A> and <R> must have the same number of indeterminates ring");
    fi;

    return GradedAlgebraPresentationNC(
      R, 
      HAPPRIME_SwitchPolynomialIndeterminates(
        BaseRing(A), R, GeneratorsOfPresentationIdeal(A)), 
      IndeterminateDegrees(A));
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="HAPPRIME_GradedAlgebraPresentationAvoidingIndeterminates_DTmanGradedAlgebraInt">
##  <ManSection>
##  <Attr Name="HAPPRIME_GradedAlgebraPresentationAvoidingIndeterminates" Arg="A, avoid"/>
##
##  <Returns>
##    GradedAlgebraPresentation
##  </Returns>
##  <Description>
##  Returns a new presentation for the graded algebra <A>A</A> which avoids the
##  indeterminates listed in <A>avoid</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_GradedAlgebraPresentationAvoidingIndeterminates,
  function(A, avoid)
    local R;
    R := PolynomialRing(
      CoefficientsRing(A), 
      Length(IndeterminatesOfGradedAlgebraPresentation(A)),
      avoid);
    return HAPPRIME_SwitchGradedAlgebraRing(A, R);
  end
);
#####################################################################


