#############################################################################
##
##  HAPPRIME - polynomials.gi
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


#####################################################################
##  <#GAPDoc Label="TermsOfPolynomial_manPolynomial">
##  <ManSection>
##  <Attr Name="TermsOfPolynomial" Arg="poly"/>
##
##  <Returns>
##    List of pairs
##  </Returns>
##  <Description>
##  Returns a list of the terms in the polynomial.
##  This list is a list of pairs of the form <C>[mon, coeff]</C> where
##  <C>mon</C> is a monomial and <C>coeff</C> is the coefficient of that
##  monomial in the polynomial. The monomials are sorted according to
##  the total degree/lexicographic order (the same as the in 
##  <Ref Oper="MonomialGrLexOrdering" BookName="ref"/>).
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(TermsOfPolynomial,
  [IsPolynomial],
  function(poly)
  
    local terms, extrep, fam, one, i, term, mon, coeff;
    
    terms := [];
    extrep := ExtRepPolynomialRatFun(poly);
    if IsEmpty(extrep) then
      return [];
    fi;
    
    fam := FamilyObj(poly);
    i := 1;
    repeat
      mon := extrep[i];
      coeff := extrep[i+1];
      
      term := [PolynomialByExtRep(fam, [mon, One(coeff)]), coeff];
      Add(terms, term);
      i := i + 2;
    until i > Length(extrep);
    return terms;
  end
);        
#####################################################################


#####################################################################
##  <#GAPDoc Label="IsMonomial_manPolynomial">
##  <ManSection>
##  <Attr Name="IsMonomial" Arg="poly" Label="for polynomial"/>
##
##  <Returns>
##    Boolean
##  </Returns>
##  <Description>
##  Returns <K>true</K> if <A>poly</A> is a monomial, i.e. the polynomial 
##  contains only one term.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(IsMonomial,
  "for polynomial",
  [IsPolynomial],
  function(poly)
  
    return Length(ExtRepPolynomialRatFun(poly)) = 2;
  end
);        
#####################################################################


#####################################################################
##  <#GAPDoc Label="UnivariateMonomialsOfMonomial_manPolynomial">
##  <ManSection>
##  <Attr Name="UnivariateMonomialsOfMonomial" Arg="mon"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  Returns a list of the univariate monomials of the largest order 
##  whose product equals <A>mon</A>.
##  The univariate monomials are sorted according to
##  their indeterminate number.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(UnivariateMonomialsOfMonomial,
  [IsPolynomial],
  function(poly)

    local extrep, fam, coeffring, one, mon, i, unimons;

    if IsOne(poly) or IsZero(poly) then
      return [poly];
    fi;

    extrep := ExtRepPolynomialRatFun(poly);
    fam := FamilyObj(poly);
    coeffring := Ring(extrep[2]);
    one := One(coeffring);
    
    if Length(extrep) > 2 or extrep[2] <> one then
      Error("<mon> must be a monomial");
    fi;

    mon := extrep[1];

    i := 1;
    unimons := [];
    repeat
      Add(unimons, PolynomialByExtRep(fam, [[mon[i], mon[i+1]], one]));
      i := i + 2;
    until i > Length(mon);
    
    return unimons;
  end
);        
#####################################################################
  

#####################################################################
##  <#GAPDoc Label="IndeterminateAndExponentOfUnivariateMonomial_manPolynomial">
##  <ManSection>
##  <Attr Name="IndeterminateAndExponentOfUnivariateMonomial" Arg="mon"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  Returns a list <C>[indet, exp]</C> where <C>indet</C> is the indeterminate 
##  of the univariate monomial <A>mon</A> and <C>exp</C> is the exponent
##  of that indeterminate in the monomial. If <A>mon</A> is an element in the 
##  coefficient ring (i.e. the monomial contains no indeterminates) then the 
##  first element will be <A>mon</A> with an exponent of zero.
##  If <A>mon</A> is not a univariate monomial, then <K>fail</K> is returned.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(IndeterminateAndExponentOfUnivariateMonomial,
  [IsPolynomial],
  function(poly)

    local extrep, one, fam;

    extrep := ExtRepPolynomialRatFun(poly);
    if IsEmpty(extrep[1]) then
      return [extrep[2], 0];
    fi;


    if Length(extrep) <> 2 then
      return fail;
    fi;
    one := One(extrep[2]);
    fam := FamilyObj(poly);

    return [PolynomialByExtRep(fam, [[extrep[1][1], 1], one]), extrep[1][2]];
  end
);        
#####################################################################
  
  


#####################################################################
##  <#GAPDoc Label="IndeterminatesOfPolynomial_manPolynomial">
##  <ManSection>
##  <Attr Name="IndeterminatesOfPolynomial" Arg="poly"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  Returns a list of the indeterminates used in the polynomial <A>poly</A>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(IndeterminatesOfPolynomial,
  [IsPolynomial],
  function(poly)

    # Is it a constant?
    if IsEmpty(ExtRepPolynomialRatFun(poly)) or 
      IsEmpty(ExtRepPolynomialRatFun(poly)[1]) then
        return [];
    else
      return IndeterminatesOfPolynomialRing(DefaultRing(poly));
    fi;
  end
);        
#####################################################################
  
  


#####################################################################
##  <#GAPDoc Label="ReduceIdeal_manPolynomial">
##  <ManSection>
##  <Heading>ReduceIdeal</Heading>
##  <Oper Name="ReduceIdeal" Arg="I, O" Label="for Ideal"/>
##  <Oper Name="ReduceIdeal" Arg="rels, O" Label="for list of relations"/>
##
##  <Returns>
##    Ideal or list 
##  </Returns>
##  <Description>
##  For an ideal <A>I</A> returns an ideal containing a reduced generating set 
##  for the ideal, i.e. one in which no monomial in a relation in <A>I</A> is 
##  divisible by the leading term of another polynomial in <A>I</A>. 
##  The monomial ordering to be used is
##  specified by <A>O</A> (see <Ref Sect="Monomial Orderings" BookName="ref"/>).
##  The ideal can instead be specified by a list of relations <A>rels</A>, 
##  in which case a reduced list of relations is returned.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallOtherMethod(ReduceIdeal, 
  "for empty ideal",
  [IsEmpty, IsMonomialOrdering],
  function(ideal, order)
    return [];
  end
);
#####################################################################
InstallOtherMethod(ReduceIdeal, 
  "for Ideal",
  [IsPolynomialRingIdeal, IsMonomialOrdering],
  function(ideal, order)
    return Ideal(
      LeftActingRingOfIdeal(ideal), 
      ReduceIdeal(GeneratorsOfIdeal(ideal), order));
  end
);
#####################################################################
InstallMethod(ReduceIdeal, 
  "for list of relations",
  [IsHomogeneousList and IsRationalFunctionCollection, IsMonomialOrdering],
  function(I, order)

    local ideal, i, len, changed, poly;

    ideal := ShallowCopy(I);

    len := Length(ideal);
    repeat
      i := 1;
      changed := false;
      while i <= len do
        poly := PolynomialReducedRemainder(
          ideal[i], ideal{Difference([1..len],[i])}, order);
        if poly <> ideal[i] then 
          changed := true;
        fi;
        if IsZero(poly) then
          ideal[i] := ideal[len];
          Unbind(ideal[len]);
          len := len - 1;
        else 
          ideal[i] := poly;
          i := i + 1;
        fi;
      od;
    until not changed;
    for i in [1..Length(ideal)] do
      ideal[i] := ideal[i] / LeadingCoefficientOfPolynomial(ideal[i], order);
    od;

    return ideal;

  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="ReducedPolynomialRingPresentation_manPolynomial">
##  <ManSection>
##  <Oper Name="ReducedPolynomialRingPresentation" Arg="R, I[, avoid]"/>
##  <Oper Name="ReducedPolynomialRingPresentationMap" Arg="R, I[, avoid]"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  For a polynomial ring <A>R</A> and a list of relations <A>I</A> in that 
##  ring, returns a list <C>[S, J]</C> representing a polynomial quotient ring 
##  <M>S/J</M> which is isomorphic to the ring <M>R/I</M>, but which involves 
##  the minimal number of ring indeterminates. The indeterminates in <C>S</C>
##  will be distinct from thise in <A>R</A>, and an optional argument 
##  <A>avoid</A> can be used to give a list of further indeterminates to avoid 
##  when creating the ring <C>S</C>.
##  <P/>
##  The extended version of this function,
##  <Ref Oper="ReducedPolynomialRingPresentationMap"/>, returns an additional
##  third element to the list, which contains two lists giving the mapping 
##  between the new ring indeterminates and the old ring indeterminates. The
##  first list is of polynomials in the original indeterminates, the 
##  second the equivalent polynomials in the new ring indeterminates.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(ReducedPolynomialRingPresentation, 
  "with nothing to avoid",
  [IsPolynomialRing, IsHomogeneousList and IsRationalFunctionCollection, IsHomogeneousList],
  function(ring, I, avoid)
    return ReducedPolynomialRingPresentationMap(ring, I, avoid){[1..2]};
  end
  );
#####################################################################
InstallOtherMethod(ReducedPolynomialRingPresentation, 
  "with nothing to avoid",
  [IsPolynomialRing, IsHomogeneousList and IsRationalFunctionCollection],
  function(ring, I)
    return ReducedPolynomialRingPresentationMap(ring, I, []){[1..2]};
  end
  );
#####################################################################
InstallOtherMethod(ReducedPolynomialRingPresentation, 
  "for empty relations",
  [IsPolynomialRing, IsEmpty, IsHomogeneousList],
  function(ring, I, avoid)
    local newring, newrelations, map;
    # Just copy the ring
    newring := PolynomialRing(
      CoefficientsRing(ring), Length(IndeterminatesOfPolynomialRing(ring)),
      Concatenation(avoid, IndeterminatesOfPolynomialRing(ring)));

    newrelations := HAPPRIME_SwitchPolynomialIndeterminates(
      ring, newring, I);

      return [newring, newrelations ];
  end
  );
#####################################################################
InstallOtherMethod(ReducedPolynomialRingPresentationMap, 
  "with nothing to avoid",
  [IsPolynomialRing, IsHomogeneousList and IsRationalFunctionCollection],
  function(ring, I)
    return ReducedPolynomialRingPresentationMap(ring, I, []);
  end
  );
#####################################################################
InstallOtherMethod(ReducedPolynomialRingPresentationMap, 
  "for empty relations and nothing to avoid",
  [IsPolynomialRing, IsEmpty],
  function(ring, I)
    return ReducedPolynomialRingPresentationMap(ring, I, []);
  end
  );
#####################################################################
InstallOtherMethod(ReducedPolynomialRingPresentationMap, 
  "for empty relations",
  [IsPolynomialRing, IsEmpty, IsHomogeneousList],
  function(ring, I, avoid)
    local newring, newrelations, map;
    # Just copy the ring
    newring := PolynomialRing(
      CoefficientsRing(ring), Length(IndeterminatesOfPolynomialRing(ring)),
      Concatenation(avoid, IndeterminatesOfPolynomialRing(ring)));

    newrelations := HAPPRIME_SwitchPolynomialIndeterminates(
      ring, newring, I);

      return [newring, newrelations, 
        [IndeterminatesOfPolynomialRing(ring), 
          IndeterminatesOfPolynomialRing(newring)] ];
  end
  );
#####################################################################
InstallMethod(ReducedPolynomialRingPresentationMap, 
  [IsPolynomialRing, IsHomogeneousList and IsRationalFunctionCollection, IsHomogeneousList],
  function(ring, I, avoid)

    local ideal, indetorder, order, indet, i, t, t2, indets, allindets, ord,
    removedindets, removedrelations, removedring, unimons, relation, poly, len, 
    nomod, newring, newrelations, map;

    # Check for unit in the ideal - if so, we return a trivial ring
    if Length(I) = 1 and IsOne(I[1]) then
      return [Ring(Zero(I[1])), [ ] ];
    fi;

    ideal := ShallowCopy(I);
    # indeterminates will be removed from indets and added to removedindets
    # as we do our reduction 
    indets := ShallowCopy(IndeterminatesOfPolynomialRing(ring));
    removedindets := [];
    removedrelations := [];

    # Are any of the leading terms of the ideal single indeterminates?
    # If so then we can (after reduction) remove those indeterminates and those 
    # terms of the ideal
    repeat
      indet := false;
      # Find a relation in the ideal that involves a solitary indeterminate
      for i in [1..Length(ideal)] do
        for t in TermsOfPolynomial(ideal[i]) do
          unimons := UnivariateMonomialsOfMonomial(t[1]);
          if Length(unimons) = 1 and unimons[1] = IndeterminateOfUnivariateRationalFunction(unimons[1]) then
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
          Add(removedrelations, relation);
          break;
        fi;
      od;

      if indet <> false then
        # Remove this indeterminate from the indets list and put it
        # onto the removedindets list instead
        Remove(indets, Position(indets, indet));
        Add(removedindets, indet);

        # And create a new (temporary) order which has indet as the most important
        order := MonomialLexOrdering(Concatenation([indet], indets));
        # And make sure that our relation has a unit coefficient
        relation := relation / LeadingCoefficientOfPolynomial(relation, order);
  
        # Now reduce all the other relations in the ideal with this one, which
        # will get rid of this indeterminate
        i := 1;
        len := Length(ideal);
        while i <= len do
          poly := PolynomialReducedRemainder(ideal[i], [relation], order);
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
  
    # Tidy up the ideal
    ideal := ReduceIdeal(ideal, MonomialLexOrdering());

    # Now make the new ring and convert the indeterminates
    newring := PolynomialRing(CoefficientsRing(ring), Length(indets),
      Concatenation(avoid, IndeterminatesOfPolynomialRing(ring)));
    newrelations := HAPPRIME_SwitchPolynomialIndeterminates(
      PolynomialRing(CoefficientsRing(ring), indets), newring, ideal);
    # And remember the mapping
    map := [indets, ShallowCopy(IndeterminatesOfPolynomialRing(newring))];

    if not IsEmpty(removedrelations) then
      # Finally, sort out the map between the old and new indeterminates
      # for the removed ones
      # Start off by adding relations that tell us what the new indeterminates are
      ideal := IndeterminatesOfPolynomialRing(newring) - indets;
      # Now add the relations that we have removed
      Append(ideal, removedrelations);
      # Create an ordering that puts the removed relations largest so that 
      # they will be removed as much as possible, and the new relations smallest
      # so they will be kept, then reduce this set of relations
      ord := MonomialLexOrdering(Concatenation(removedindets, indets, 
        IndeterminatesOfPolynomialRing(newring)));
      ideal := ReduceIdeal(ideal, ord);

      # Now add to the map. The polynomials in the ideal that have
      # leading monomials which still involve the removedindets are the ones
      # we want
      removedring := PolynomialRing(CoefficientsRing(ring), removedindets);
      for i in ideal do
        if LeadingMonomialOfPolynomial(i, ord) in removedring then
          # is this a single indeterminate (i.e. a relation of the form x_i)?
          if Length(TermsOfPolynomial(i)) = 1 then
            poly := [[i], [Zero(i)]];
          else
            poly := [[], []];
            for t in TermsOfPolynomial(i) do
              if t[1] in ring then 
                Add(poly[1], t[1]*t[2]);
              elif t[1] in newring then
                Add(poly[2], t[1]*t[2]);
              else
                Error("Relation is not seperable when creating map. Please consult the package maintainer.");
              fi;
            od;
            if IsEmpty(poly[2]) then
              Error("Unexpected relation with no new indeterminates. Please consult the package maintainer.");
            fi;
          fi;
          Add(map[1], Sum(poly[1]));
          Add(map[2], Sum(poly[2]));
        fi;
      od;
    fi;

    return [newring, newrelations, map];
  
  end
);
#####################################################################


#################################
## This functions new to GAP 4.4.10
## so needs defining if using an earlier version
if not IsBound(EmptyPlist) then
  EmptyPlist := function(n)
    return [];
  end;
fi;
#################################


#####################################################################
##  <#GAPDoc Label="HAPPRIME_SwitchPolynomialIndeterminates_manDTPolynomialInt">
##  <ManSection>
##  <Func Name="HAPPRIME_SwitchPolynomialIndeterminates" Arg="R, S, poly"/>
##
##  <Returns>
##    Polynomial or List of Polynomials
##  </Returns>
##  <Description>
##  Changes the indeterminates in <A>poly</A>, which should be a polynomial or 
##  a list of polynomials, substituting the indeterminates of the polynomial 
##  ring <A>S</A> one-for-one for those in <A>R</A> (from which all polynomials 
##  in <A>poly</A> must come). The returned object is either a polynomial or a 
##  list of polynomials in the new indeterminates, depending on the input object.
##  <P/>
##  See <Ref Func="HAPPRIME_MapPolynomialIndeterminates"/> for a function that
##  can work with more general indeterminate maps.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_SwitchPolynomialIndeterminates,
  function(R, S, polys)
    local onepoly, Rindetnums, Sindetnums, newpolys, fam, p, extrep, extrepnew,
      i, j, e;
    
    if CoefficientsRing(R) <> CoefficientsRing(S) then
      Error("<R> and <S> must have the same coefficient ring");
    fi;

    onepoly := false;
    if Length(IndeterminatesOfPolynomialRing(R)) <> 
      Length(IndeterminatesOfPolynomialRing(S)) then 
      Error("<R> and <S> must have the same number of indeterminates");
    fi;
    
    if IsEmpty(polys) then
      return polys;
    fi;
    if IsPolynomial(polys) then
      polys := [polys];
      onepoly := true;
    elif not IsHomogeneousList(polys) then
      Error("<polys> must be a polynomial or a list of polynomials in <R>");
    fi;

    Rindetnums := List(IndeterminatesOfPolynomialRing(R), 
      IndeterminateNumberOfUnivariateRationalFunction);
    Sindetnums := List(IndeterminatesOfPolynomialRing(S), 
      IndeterminateNumberOfUnivariateRationalFunction);

    newpolys := EmptyPlist(Length(polys));
    fam := FamilyObj(polys[1]);
    for p in polys do
      if not p in R then
        Error("<polys> must be a polynomial or a list of polynomials in <R>");
      fi;

      extrep := ExtRepPolynomialRatFun(p);
      extrepnew := EmptyPlist(Length(extrep));
      i := 1;
      repeat
        e := ShallowCopy(extrep[i]);
        if not IsEmpty(e) then
          j := 1;
          repeat
            # Swap the indeterminate number
            e[j] := Sindetnums[Position(Rindetnums, e[j])];
            j := j+2;
          until j > Length(e);
        fi;
        Add(extrepnew, e);
        Add(extrepnew, extrep[i+1]);
        i := i+2;
      until i > Length(extrep);
      Add(newpolys, PolynomialByExtRep(fam, extrepnew));
    od;


    # And sort out the return type
    if onepoly then
      return newpolys[1];
    else
      return newpolys;
    fi;
end);
######################################################


#####################################################################
##  <#GAPDoc Label="HAPPRIME_MapPolynomialIndeterminates_manDTPolynomialInt">
##  <ManSection>
##  <Func Name="HAPPRIME_MapPolynomialIndeterminates" Arg="old, new, poly"/>
##
##  <Returns>
##    Polynomial or List of Polynomials
##  </Returns>
##  <Description>
##  Changes the indeterminates in <A>poly</A>, which can be a polynomial or a 
##  list of polynomials, substituting the polynomials in <A>old</A> for those
##  in <A>new</A>. The returned object is either a polynomial or a list of 
##  polynomials in the new indeterminates, depending on the input object.
##  The change of variable arguments, <A>old</A> and <A>new</A>, do not 
##  have to be simply indeterminates: they can be can be lists of polynomials 
##  which are equivalent in the two different sets of indeterminates.
##  If a polynomial cannot be converted (i.e. if it cannot be generated from the
##  polynomials in <A>old</A>) then <K>fail</K> is returned for that polynomial.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_MapPolynomialIndeterminates,
  function(old, new, polys)
  
    local old2, new2, oldindets, newindets, ord, newpolys, p, newp, onepoly, I, 
      lead, newI, ord2;
    
    onepoly := false;
    if not IsHomogeneousList(old) then
      Error("<old> must a list of polynomials");
    fi;
    if not IsHomogeneousList(new) then
      Error("<new> must a list of polynomials");
    fi;
    if Length(old) <> Length(new) then
      Error("<old> and <new> must both be the same length");
    fi;
    if IsPolynomial(polys) then
      polys := [polys];
      onepoly := true;
    elif not IsHomogeneousList(polys) then
      Error("<polys> must be a polynomial or a list of polynomials in the old indeterminates");
    fi;

    oldindets := [];
    for p in old do
      UniteSet(oldindets, IndeterminatesOfPolynomial(p));
    od;
    newindets := [];
    for p in new do
      UniteSet(newindets, IndeterminatesOfPolynomial(p));
    od;
    if not IsEmpty(Intersection(oldindets, newindets)) then
      Error("the indeterminates in <old> and <new> must be distinct sets.");
    fi;
  
    # Order the indeterminates so that the old ones are larger, so will
    # be replaced
    ord := MonomialLexOrdering(Concatenation(oldindets, newindets));

    # Are there any zeros in the old list? If so, remove them and the 
    # corresponding element of the new list
    old2 := ShallowCopy(old);
    new2 := ShallowCopy(new);
    repeat
      p := Position(old2, Zero(old2[1]));
      if p <> fail then
        Remove(old2, p);
        Remove(new2, p);
      fi;
    until p = fail;
    # The conversion is relations in an ideal. Make sure it is a GroebnerBasis
    I := HAPPRIME_SingularGroebnerBasis(old2 - new2, ord);
    
    # Extract out the relations in the new ring
    newI := Filtered(I, i->IsSubset(newindets, IndeterminatesOfPolynomial(i)));
    # and turn this into a Groebner Basis with grevlex ordering
    ## TODO I would like to make it MonomialGrevlexOrdering here,
    ## but a bug in that function (reported by me on 27/6/08) means
    ## that that ordering doesn't work!
    #ord2 := MonomialGrevlexOrdering(); 
    ord2 := MonomialGrlexOrdering();
    newI := HAPPRIME_SingularGroebnerBasis(newI, ord2);

    newpolys := [];
    for p in polys do
      # is this just a constant?
      lead := LeadingMonomial(p);
      if IsEmpty(lead) or lead[2] = infinity or lead[2] = 0 then
        Add(newpolys, p);
      else
        # Reduce this. Given the ordering, this will substitute the 
        # new in place of the old
        # We further reduce this using the grevlex ordering to get the 
        # simplest (smallest degree) form in the new indeterminates
        newp := PolynomialReducedRemainder(
          PolynomialReducedRemainder(p, I, ord), newI, ord2);
        # Check that this really is in the new indeterminates - it might not
        # be if the map is not complete
        if not IsSubset(newindets, IndeterminatesOfPolynomial(newp)) then
          Add(newpolys, fail);
        else
          Add(newpolys, newp);
        fi;
      fi;
    od;

    # And sort out the return type
    if onepoly then
      return newpolys[1];
    else
      return newpolys;
    fi;
end);
######################################################


#####################################################################
##  <#GAPDoc Label="HAPPRIME_CombineIndeterminateMaps_manDTPolynomialInt">
##  <ManSection>
##  <Func Name="HAPPRIME_CombineIndeterminateMaps" Arg="coeff, M, N"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  Returns the indeterminate map that results from applying map <A>M</A>
##  followed by map <A>N</A>. An indeterminate map is a list containing two
##  lists, the first of which is a list of polynomials in the original indeterminates, 
##  the second the equivalent polynomials in the new ring indeterminates.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_CombineIndeterminateMaps,
  function(M, N)
      
    local indetsA, indetsB, indetsC, i, ord, I, J, relations, map, poly, t;
    
    # Recreate the three rings
    indetsA := [];
    indetsB := [];
    indetsC := [];
        
    for i in M[1] do
      Append(indetsA, IndeterminatesOfPolynomial(i));
    od;
    indetsA := AsSet(indetsA);
    
    for i in M[2] do
      Append(indetsB, IndeterminatesOfPolynomial(i));
    od;
    indetsB := AsSet(indetsB);
    
    for i in N[1] do
      Append(indetsC, IndeterminatesOfPolynomial(i));
    od;
    indetsC := AsSet(indetsC);
    
    # Check that the second ring in M and the first in N are compatible
    if not IsSubset(indetsC, indetsB) and not IsSubset(indetsB, indetsC) then
      Error("the maps <M> and <N> are incompatible");
    fi;
        
    indetsC := [];
    for i in N[2] do
      Append(indetsC, IndeterminatesOfPolynomial(i));
    od;
    indetsC := AsSet(indetsC);
    
    # And check that the three rings are distinct
    if not IsEmpty(Intersection(indetsA, indetsB)) then 
      Error("the source and target rings of <M> are not different");
    fi;
    if not IsEmpty(Intersection(indetsA, indetsC)) then 
      Error("the source ring of <M> and target ring of <N> are not different");
    fi;
    if not IsEmpty(Intersection(indetsB, indetsC)) then 
      Error("the source and target rings of <N> are not different");
    fi;
        
    # Order the middle indeterminates first so that they are removed
    ord := MonomialLexOrdering(Concatenation(indetsB, indetsA, indetsC));
    # And get the two set of relations
    I := M[1] - M[2];
    J := N[1] - N[2];
    # Now reduce each with the other
    relations := [];
    for i in I do
      Add(relations, PolynomialReducedRemainder(i, J, ord));
    od;
    for i in J do
      Add(relations, PolynomialReducedRemainder(i, I, ord));
    od;
    relations := AsSet(relations);

    # Now create the new map. The polynomials in the ideal have
    # leading monomials which involve the source ring
    map := [[], []];
    for i in relations do
      poly := [[Zero(i)], [Zero(i)]];
      for t in TermsOfPolynomial(i) do
        if IsSubset(indetsA, IndeterminatesOfPolynomial(t[1])) then 
          Add(poly[1], t[1]*t[2]);
        elif IsSubset(indetsC, IndeterminatesOfPolynomial(t[1])) then
          Add(poly[2], t[1]*t[2]);
        else
          # This relation still involves an old indeterminate, so we
          # can't do anything with it - ignore it
          poly := "ignore";
          break;
        fi;
      od;
      if poly <> "ignore" then
        Add(map[1], Sum(poly[1]));
        Add(map[2], Sum(poly[2]));
      fi;
    od;

    return map;
end);
######################################################

    
#####################################################################
##  <#GAPDoc Label="HAPPRIME_SingularGroebnerBasis_manDTGroebnerInt">
##  <ManSection>
##  <Func Name="HAPPRIME_SingularGroebnerBasis" Arg="pols, O"/>
##  <Func Name="HAPPRIME_SingularReducedGroebnerBasis" Arg="pols, O"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  Returns the Gröbner basis (or reduced Gröbner basis) with respect to the 
##  ordering <A>O</A> for the ideal generated by the polynomials <A>pols</A>. 
##  This function uses the Gröbner basis implementation 
##  from &singular;, for preference, if available 
##  (<Ref Func="GroebnerBasis" BookName="singular"/>), and if so it also 
##  manuipulates the result to fix a bug in &singular; where the returned 
##  polynomials are not necessarily returned with a value external 
##  representation (see 
##  <Ref Sect="The Defining Attributes of Rational Functions" BookName="ref"/>).
##  <P/>
##  If the option <C>obeyGBASIS</C> is <K>true</K>, then this function will use 
##  whichever algorithm is specified by the <K>GBASIS</K> global variable
##  (see <Ref Var="SINGULARGBASIS" BookName="singular"/>).
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
#if LoadPackage("singular") = true then
if IsPackageMarkedForLoading("singular","0") then
  InstallGlobalFunction(HAPPRIME_SingularGroebnerBasis,
    function(pols, O)
      local gbasis, ideal, fam;

      # If empty, do nothing
      if IsEmpty(pols) then
        return pols;
      fi;
      # Make sure we use Singular
      gbasis := GBASIS;
      if not ValueOption("obeyGBASIS") = true then
        GBASIS := SINGULARGBASIS;
      fi;
      ideal := GroebnerBasis(pols, O);
      GBASIS := gbasis;

      # Hack to make sure that the polynomials returned by Singular are correct
      fam := FamilyObj(pols[1]);
      return List(ideal, x->PolynomialByExtRep(fam, ExtRepPolynomialRatFun(x)));;
    end
  ); 
  #####################################################################
  InstallGlobalFunction(HAPPRIME_SingularReducedGroebnerBasis,
    function(pols , O)
      local  ipr, mcf, R, I, input, out, fam;
 
    # If empty, do nothing
    if IsEmpty(pols) then
      return pols;
    fi;
    
    if ValueOption("obeyGBASIS") = true and GBASIS = GAPGBASIS then
      return ReducedGroebnerBasis(pols, O);
    fi;
    
    if IsPolynomialRingIdeal(pols)  then
      R := LeftActingRingOfIdeal(pols);
      pols := GeneratorsOfTwoSidedIdeal(pols);
    else
      R := DefaultRing(pols);
    fi;
    
    if IsMonomialOrdering(O) then
      ipr := ShallowCopy(IndeterminatesOfPolynomialRing(R));
      mcf := MonomialComparisonFunction(O);
      Sort(ipr, mcf);
      ipr := Reversed(ipr);
      R := PolynomialRing(LeftActingDomain(R), ipr);
    fi;
    
    if not(HasTermOrdering(R) and 
      IsIdenticalObj(TermOrdering(R), O)) then
        SetTermOrdering(R, O);
        SingularSetBaseRing(R);
    fi;
    
    I := Ideal(R, pols);

    # Now use Singular to find the Groebner Basis
    Info( InfoSingular, 2, "running GroebnerBasis..." );
  
    SingularCommand("", "option(redSB)");
    # preparing the input for Singular
    input := "";
  
    Append( input, "ideal GAP_groebner = simplify( groebner( " );
    Append( input, ParseGapIdealToSingIdeal( I ) );
    Append( input, " ), 1 );\n" );
  
    out := SingularCommand( input, "string (GAP_groebner)" );
  
    SingularCommand("", "option(noredSB)");
    Info( InfoSingular, 2, "done GroebnerBasis." );
  
    I := List( SplitString( out, ',' ), ParseSingPolyToGapPoly );
    
    # Hack to make sure that the polynomials returned by Singular are correct
    fam := FamilyObj(pols[1]);
    return List(I, x->PolynomialByExtRep(fam, ExtRepPolynomialRatFun(x)));;
  end );
else
  InstallGlobalFunction(HAPPRIME_SingularGroebnerBasis,
    function(I, O)
      return GroebnerBasis(I, O);
    end
  ); 
  InstallGlobalFunction(HAPPRIME_SingularReducedGroebnerBasis,
    function(I, O)
      return ReducedGroebnerBasis(I, O);
    end
  ); 
fi;
#####################################################################


