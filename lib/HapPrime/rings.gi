#############################################################################
##
##  HAPPRIME - rings.gi
##  Functions, Operations and Methods for dealing with cohomology rings
##  Paul Smith
##
##  Copyright (C)  2007-2008
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
##  <#GAPDoc Label="PresentationOfGradedStructureConstantAlgebra_DTmanRings_Other">
##  <ManSection>
##  <Attr Name="PresentationOfGradedStructureConstantAlgebra" Arg="A"/>
##
##  <Returns>
##    <K>GradedAlgebraPresentation</K>
##  </Returns>
##  <Description>
##  Returns a ring presentation for the graded algebra <A>A</A>.
##  The ring <A>A</A> must be a structure constant algebra with embedded
##  degrees, such as is returned by 
##  <Ref Func="ModPCohomologyRing" BookName="HAP"/>.
##  The generators of the <K>GradedAlgebraPresentation</K> (as returned by
##  <Ref Func="IndeterminatesOfGradedAlgebraPresentation"/> are in
##  one-to-one correspondance with the generators of <A>A</A> as returned by
##  <Ref Func="ModPRingGenerators" BookName="HAP"/> (ignoring the first 
##  generator, which is in degree zero).
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(PresentationOfGradedStructureConstantAlgebra,
  [IsAlgebra],
  function(A)
    local field, zero, one, basis, dim, 
    degrees, maxdegree, degreesize, degreeindex, degreeindices,
    gens, elm, polyring, indets, i,
    newbasis, order, orderfunc, leading_monomials, relations, 
    d, done_xy, degree_xy, x, y, xy, prod,
    NS, coeffs, rels, r,
    ElmDegree;
  
#################################################
##Added by Graham Ellis
if not "degree" in NamesOfComponents(A)  then
Print("This function only applies to commutative graded algebras\n");
return fail;
fi;

if not IsCommutative(A) then 
Print("This function only applies to commutative graded algebras\n");
return fail;
fi;
##
#################################################
    # What is the field?
    field := LeftActingDomain(A);
    zero := Zero(field);
    one := One(field);
  
    # Get the basis for the vector space
    basis := Basis(A);
    dim := Size(basis);

    # Make a list of the degrees of the elements
    # This is evil - degree should not be accessed as a direct member of A
    # Check that degree can be accessed
    if not IsFunction(A!.degree) then
      Error("the algebra <A> must have an A!.degree() function");
    fi;
    # And create a nicer function to access the degree  
    ElmDegree := function(A, e)
      return A!.degree(e);
    end;
    # What's the largest degree?
    degrees := List([1..dim], i->ElmDegree(A, basis[i]));
    maxdegree := Maximum(degrees);
    # We want to index the basis elements of each degree seperately, i.e.
    # the first element of degree n, the second element of degree n and so on
    # Make a list mapping total basis number to number in that degree
    # How many elements of each degree are there in degree n-1? (These are the 
    # Betti numbers)
    degreesize := ListWithIdenticalEntries(maxdegree+1, 0);
    degreeindex := ShallowCopy(degrees);
    for i in [1..dim] do
      d := ElmDegree(A, basis[i]);
      degreesize[d+1] := degreesize[d+1] + 1;
      degreeindex[i] := degreesize[d+1];
    od;
    degreeindices := List([1..maxdegree+1], i->[]);
    for i in [1..dim] do
      Add(degreeindices[ElmDegree(A, basis[i])+1], i);
    od;

    # And get the generators
    gens := ShallowCopy(ModPRingGenerators(A));
    # Remove the first of these, which will be in degree zero
    elm := Remove(gens, 1);
    if ElmDegree(A, elm) <> 0 then
      Error("the degree of the first generator returned by ModPRingGenerators() is not zero");
    fi;
  
    # Create the polynomial ring
    polyring := PolynomialRing(field, Length(gens));
    indets := IndeterminatesOfPolynomialRing(polyring);

    # Convert the generators list into triplets of
    # [degree, indet, algebra element]
    for i in [1..Length(gens)] do
      gens[i] := [ElmDegree(A, gens[i]), indets[i], gens[i]];
    od; 

    # Our new basis is going to be built up solely of our basis elements
    # Start with the basis elements of order 1
    newbasis := [];
    for i in gens do
      if i[1] = 1 then
        Add(newbasis, i);
      fi;
    od;
  
    # We want to make relations in each degree in turn using our basis elements
    order := MonomialLexOrdering();
    orderfunc := MonomialComparisonFunction(order);
    leading_monomials := [];
    relations := [];
    for d in [2..maxdegree] do
      done_xy := [];   # Products that have already been considered
      degree_xy := []; # Products in this degree to do linear algebra on later
      for x in newbasis do   # Step through all of our basis elements
        for y in gens do  # And consider multiplying by each generator
          # We only want multiples which give us the current degree
          if x[1] + y[1] = d then  

            # Multiply the indeterminants
            xy := x[2]*y[2];
            # Have we already considered this one?
            if not xy in done_xy then
              Add(done_xy, xy);
              # Is this one a multiple of something we've already considered?
              if not IsZero(PolynomialReducedRemainder(xy, leading_monomials, order)) then
                # OK - it seems to be original. What is the product according 
                # to the algebra?

                prod := x[3]*y[3];
                if IsZero(prod) then
                  # The product is zero. This is immediately one of our relations
                  Add(relations, xy);
                  Add(leading_monomials, xy);
                else
                  # Otherwise remember this for a little later, as it
                  # might form part of our basis for this degree

                  # What is the product in terms of coefficients of our basis elements?
                  coeffs := Coefficients(basis, prod);
                  # This will only be made up of basis elements of this degree, so
                  # ignore everything else
                  coeffs := coeffs{degreeindices[d+1]};

                  # Remember useful information for later
                  Add(degree_xy, [xy, coeffs, prod]); 
                fi;
              fi;
            fi; 
          fi;
        od; # for y
      od; # for x
    
      if not IsEmpty(degree_xy) then
        # Sort the degree_xy list into decreasing order of leading monomial
        # This will guarantee that when I take my TriangulizedNullspace, each 
        # relation will have a different leading monomial
      Sort(degree_xy, function(u, v) return u[1] > v[1]; end);

      # See whether any combinations of our products in this degree
      # can be combined to make zero
      NS := TriangulizedNullspaceMat(List(degree_xy, x->x[2]));
      # What are these in terms of relations?
      rels := List(NS, ns->Sum([1..Size(degree_xy)], i->ns[i]*degree_xy[i][1]));
      for r in rels do
        Add(relations, r);
        Add(leading_monomials, LeadingMonomialOfPolynomial(r, order));
      od;

      # Any of my entries in degree_xy that have led to a relation are
      # clearly not independent, and can be ignored. The other ones
      # will form a (possibly partial) basis for this degree, so remember 
      # them in my new basis
      for elm in degree_xy do
        if not elm[1] in leading_monomials then
          # newbasis is triplets of
          # [degree, indet, algebra element]
          Add(newbasis, [d, elm[1], elm[3]]);
        fi;
      od;
    fi;

    # The rest of the basis for this degree cannot come from products of 
    # current elements, so must come (by definition) from the generators
    # in this degree. Add these to our basis as well.
    Append(newbasis, Filtered(gens, i->i[1] = d));

  od; # for d

  # Sort the relations to make them pretty
  Sort(relations);
  # And get the degrees of the generators
  degrees := List(gens, i->i[1]);

  return GradedAlgebraPresentation(polyring, relations, degrees);
end
); 
#####################################################################

  
#####################################################################
##  <#GAPDoc Label="Mod2CohomologyRingPresentation_manRings">
##  <ManSection>
##  <Attr Name="Mod2CohomologyRingPresentation" Arg="G" Label="for group"/>
##  <Oper Name="Mod2CohomologyRingPresentation" Arg="G, n" Label="for group and degree"/>
##  <Oper Name="Mod2CohomologyRingPresentation" Arg="R" Label="for resolution"/>
##  <Oper Name="Mod2CohomologyRingPresentation" Arg="A" Label="for algebra"/>
##
##  <Returns>
##    <K>GradedAlgebraPresentation</K>
##  </Returns>
##  <Description>
##  Calculates and returns a cohomology ring presentation for the
##  group <M>G</M>. See 
##  <Ref Chap="Presentations of graded algebras" BookName="HAPprime Datatypes"/>
##  in the datatypes reference manual for details of the 
##  <K>GradedAlgebraPresentation</K> type.
##  <P/>
##  If the only argument is a <M>p</M>-group <A>G</A> then this function 
##  computes and returns the provably-correct cohomology ring presentation. This 
##  version first computes the Lyndon-Hoschild-Serre Spectral Sequence until 
##  convergence to find the additive structure of the cohomology ring, and then 
##  computes the cohomology ring up to and including the maximum necessary 
##  generator or relation, using the <C>(G, n)</C> method described below.
##  For certain groups, the cohomology ring is returned without computation:
##  the known mod-&p; cohomology ring presentation 
##  for cyclic groups is returned without calculation, and for 
##  groups which can be expressed as a direct product, the cohomology ring
##  is computed as a tensor product of its direct factors (thus the 
##  cohomology ring of all Abelian groups are also returned with minimal
##  computation.)
##  <P/> 
##  When given a <M>p</M>-group <A>G</A> and integer <A>n</A>, this function 
##  computes the presentation modulo all elements of degree greater <A>n</A>.
##  Alternatively, a minimal resolution <A>R</A> (with <A>n</A> terms) can be 
##  input, or a structure constant algebra <A>A</A> with embedded degrees
##  (from <Ref Func="ModPCohomologyRing" BookName="HAP"/>).
##  <P/>
##  See Section <Ref Sect="ExPresentation"/> and <Ref Sect="ExPresentationLHS"/>
##  for examples and more description. See also 
##  <Ref Func="LHSSpectralSequence" BookName="HAPprime Datatypes"/> for
##  details of options that can be used to guide the spectral sequence
##  computation.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(Mod2CohomologyRingPresentation,
  "for group",
  [IsGroup],
  function(G)
    local A, P, n;
    
####
#### Added by Graham Ellis (22/12/2008)
    if SSortedList(Factors(Order(G)))<>[2] then 
    Print("This function is currently implemented only for 2-groups.\n");
    return fail;
    fi;
Print("Alpha version of completion test code will be used. This needs further work.\n");  
####
####

    # See if there is an easy answer
    A := HAPPRIME_CohomologyRingWithoutResolution(G);
    if A <> fail then
      return A;
    fi;

    A := LHSSpectralSequenceLastSheet(G);
    if A = fail then 
      return fail;
    fi;
    n := MaximumDegreeForPresentation(A);
    P := Mod2CohomologyRingPresentation(G, n);
    
    # Sanity check: make sure that the Poincare series is the same
    if HilbertPoincareSeries(A) <> HilbertPoincareSeries(P) then
      Error("Panic: Poincare series of the last sheet of the LHS is not the same as that of the cohomology ring");
    fi;
    
    return P;
  end
);
#####################################################################
InstallMethod(Mod2CohomologyRingPresentation,
  "for group and degree",
  [IsGroup, IsPosInt],
  function(G, n)

####
#### Added by Graham Ellis (22/12/2008)
    if SSortedList(Factors(Order(G)))<>[2] then
    Print("This function is currently implemented only for 2-groups.\n");
    return fail;
    fi;
####
####

    return PresentationOfGradedStructureConstantAlgebra(ModPCohomologyRing(G, n));
  end
);
#####################################################################
InstallMethod(Mod2CohomologyRingPresentation,
  "for resolution",
  [IsHapResolution],
  function(R)
####
#### Added by Graham Ellis (22/12/2008)
    if SSortedList(Factors(Order(GroupOfResolution(R))))<>[2] then
    Print("This function is currently implemented only for 2-groups.\n");
    return fail;
    fi;
####
####

    return PresentationOfGradedStructureConstantAlgebra(ModPCohomologyRing(R));
  end
);
#####################################################################
InstallMethod(Mod2CohomologyRingPresentation,
  "for algebra",
  [IsAlgebra],
  function(A)
####
#### Added by Graham Ellis (22/12/2008)
    if Characteristic(LeftActingDomain(A))<>2 then
    Print("This function is currently implemented only for 2-groups.\n");
    return fail;
    fi;
####
####

    return PresentationOfGradedStructureConstantAlgebra(A);
  end
);
#####################################################################

#####################################################################
##  <#GAPDoc Label="LHSSpectralSequence_DTmanGradedAlgebra_Func">
##  <ManSection>
##  <Oper Name="LHSSpectralSequence" Arg="G [, N], n"/>
##  <Oper Name="LHSSpectralSequenceLastSheet" Arg="G [, N]"/>
##
##  <Returns>
##    <K>GradedAlgebraPresentation</K> or list
##  </Returns>
##  <Description>
##  Computes the Lyndon-Hoschild-Serre spectral sequence for the group 
##  extension <M>N \to G \to G/N</M>. 
##  If a normal suggroup <A>N</A> is not provided, then the largest central 
##  subgroup of <M>G</M> is used, or (if the order of the centre is larger
##  than <M>\sqrt{|G|}</M>) then the
##  central subgroup that leads to the smallest initial sheet size is chosen.
##  <P/>
##  The function <Ref Oper="LHSSpectralSequence"/> returns the first <A>n</A> 
##  sheets of the spectral sequence, or all of the sequence up to convergence,
##  if that occurs before the <M>(n+1)</M>th sheet. The Lyndon-Hoschild-Serre
##  spectral sequence starts at the <M>E_2</M> sheet, so the first element in
##  returned list will always be empty. If <A>n</A> is set to 
##  <K>infinity</K> then the length of the returned list equals the number of 
##  sheets for convergence, and the last sheet in the list is the limiting 
##  sheet.
##  <P/>
##  The function <Ref Oper="LHSSpectralSequenceLastSheet"/> returns only the 
##  limiting sheet of the spectral sequence. This ring is
##  an associated graded algebra of the mod-<M>p</M> cohomology ring of <M>G</M>,
##  with the same additive structure while not necessarily being isomorphic to 
##  it.
##  <P/>
##  There are four options <Ref Chap="Options Stack" BookName="ref"/> which can 
##  be used to guide this algorithm:
##  <List>
##    <Item><C>InitialLHSBicomplexSize</C> can be used to specify the
##      initial size of the bicomplex (the default is 5). If, in the process of 
##      computing the spectral sequence, this is found to be too small then the 
##      algorithm restarts with a larger value. Specifying a larger initial
##      value in these cases can save time.</Item>
##    <Item><C>LargerLHSBicomplexBreak</C> if set to <K>true</K> will 
##      force the calculation to enter a break loop before restarting with a 
##      larger bicomplex, should the bicomplex be found to be too small. The user
##      user is prompted to type <C>return;</C> before continuing. The default
##      behaviour is <K>false</K>, i.e. no prompt.</Item>
##    <Item><C>LargerLHSBicomplexFail</C> if set to <K>true</K> will 
##      return <K>fail</K> should the bicomplex be found to be too small. The 
##      default behaviour is <K>false</K>, i.e. to either restart or prompt,
##      depending on the setting of the previous option.</Item>
##    <Item><C>NoInductiveProof</C> if set to <K>true</K> will not check that 
##      the cohomology rings for <M>N</M> and <M>G/N</M> are correct. Instead,
##      it will compute the cohomology rings only up to the degree needed for
##      the bicomplex size (5 by default, or specified by the 
##      <C>InitialLHSBicomplexSize</C> option).</Item>
##  </List>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallOtherMethod(LHSSpectralSequence,
  "method without subgroup",
  [IsGroup, IsPosInt],
  function(G, n)
    local N;
####
#### Added by Graham Ellis (22/12/2008)
    if SSortedList(Factors(Order(G)))<>[2] then
    Print("This function is currently implemented only for 2-groups.\n");
    return fail;
    fi;
####
####

 N := PCentre(G);
    N:=Filtered(GeneratorsOfGroup(N),x->Order(x)>1);
    N:=Group(N[1]);
    # If the centre is larger than sqrt(order(G)) then we might be able
    # to find a better split
#    if Order(N)*Order(N) > Order(G) then 
    # TODO: Small subgroups of the centre don't seem to work. Work out why
    # For the moment, only do this if N = G (in which case we shouldn't be
    # using the LHS for the cohomology anyway, since the answer is seperable
#Graham    if N = G then
#Graham      N := BestCentralSubgroupForResolutionFiniteExtension(G);
#Graham    fi;
    return LHSSpectralSequence(G, N, n);
  end
);
#####################################################################
InstallMethod(LHSSpectralSequence,
  "method for group and normal subgroup",
  [IsGroup, IsGroup, IsPosInt],
  function(G, N, n)

####
#### Added by Graham Ellis (22/12/2008)
    if SSortedList(Factors(Order(G)))<>[2] then
    Print("This function is currently implemented only for 2-groups.\n");
    return fail;
    fi;
####
####

    return HAPPRIME_LHSSpectralSequence(G, N, n);
  end
);
#####################################################################
InstallOtherMethod(LHSSpectralSequence,
  "method without subgroup",
  [IsGroup, IsInfinity],
  function(G, n)
    local N;

####
#### Added by Graham Ellis (22/12/2008)
    if SSortedList(Factors(Order(G)))<>[2] then
    Print("This function is currently implemented only for 2-groups.\n");
    return fail;
    fi;
####
####
    N := PCentre(G);
    N:=Filtered(GeneratorsOfGroup(N),x->Order(x)>1);
    N:=Group(N[1]); 
    # If the centre is larger than sqrt(order(G)) then we might be able
    # to find a better split
#    if Order(N)*Order(N) > Order(G) then 
    # TODO: Small subgroups of the centre don't seem to work. Work out why
    # For the moment, only do this if N = G (in which case we shouldn't be
    # using the LHS for the cohomology anyway, since the answer is seperable
#Graham    if N = G then
#Graham      N := BestCentralSubgroupForResolutionFiniteExtension(G);
#Graham    fi;
    return LHSSpectralSequence(G, N, n);
  end
);
#####################################################################
InstallOtherMethod(LHSSpectralSequence,
  "method for group and normal subgroup",
  [IsGroup, IsGroup, IsInfinity],
  function(G, N, n)
####
#### Added by Graham Ellis (22/12/2008)
    if SSortedList(Factors(Order(G)))<>[2] then
    Print("This function is currently implemented only for 2-groups.\n");
    return fail;
    fi;
####
####

    return HAPPRIME_LHSSpectralSequence(G, N, n);
  end
);
#####################################################################
InstallOtherMethod(LHSSpectralSequenceLastSheet,
  "method without subgroup",
  [IsGroup],
  function(G)
    local N;
####
#### Added by Graham Ellis (22/12/2008)
    if SSortedList(Factors(Order(G)))<>[2] then
    Print("This function is currently implemented only for 2-groups.\n");
    return fail;
    fi;
####
####

 N := PCentre(G);
    N:=Filtered(GeneratorsOfGroup(N),x->Order(x)>1);
    N:=Group(N[1]);
    # If the centre is larger than sqrt(order(G)) then we might be able
    # to find a better split
#    if Order(N)*Order(N) > Order(G) then 
    # TODO: Small subgroups of the centre don't seem to work. Work out why
    # For the moment, only do this if N = G (in which case we shouldn't be
    # using the LHS for the cohomology anyway, since the answer is seperable
#Graham    if N = G then
#Graham      N := BestCentralSubgroupForResolutionFiniteExtension(G);
#Graham    fi;
    return LHSSpectralSequenceLastSheet(G, N);
  end
);
#####################################################################
InstallMethod(LHSSpectralSequenceLastSheet,
  "method for group and normal subgroup",
  [IsGroup, IsGroup],
  function(G, N)
####
#### Added by Graham Ellis (22/12/2008)
    if SSortedList(Factors(Order(G)))<>[2] then
    Print("This function is currently implemented only for 2-groups.\n");
    return fail;
    fi;
####
####

    return HAPPRIME_LHSSpectralSequence(G, N, "Einf");
  end
);
#####################################################################



#####################################################################
##  <#GAPDoc Label="ModPRingGeneratorDegrees_DTmanRings_Dat">
##  <ManSection>
##  <Attr Name="ModPRingGeneratorDegrees" Arg="A"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  Returns a list containing the degree of each generator in a minimal 
##  generating set for the mod-&p; cohomology ring <A>A</A>.
##  The <M>i</M>th degree in the list corresponds to the <M>i</M>th 
##  generator returned by <Ref Func="ModPRingGenerators" BookName="HAP"/>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
 InstallMethod(ModPRingGeneratorDegrees,
   [IsAlgebra],
   function(A)
     return List(ModPRingGenerators(A), i->A!.degree(i));
   end
 );
#####################################################################

#####################################################################
##  <#GAPDoc Label="ModPRingNiceBasis_DTmanRings_Dat">
##  <ManSection>
##  <Attr Name="ModPRingNiceBasis" Arg="A"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  Returns the information needed to convert the basis for <A>A</A>
##  (see <Ref Func="Basis" BookName="ref"/>)
##  into a nicer basis consisting of only products of ring generators.
##  The function returns a pair of lists <C>[Coeff, Bas]</C>. The list <C>Coeff</C> is
##  a change-of-basis matrix, where the <M>i</M>th row gives the standard
##  basis element <M>i</M> in terms of the nice basis. The list 
##  <C>Bas</C> can be used to form the new basis and is a list of integers
##  where the <M>i</M>th 'nice basis' element is given by
##  <C>Product(List(Bas[i], x->Basis(A)[x]))</C>. 
##  <P/>
##  This attribute returns exactly the same list as is provided by 
##  component <C>A!.niceBasis</C> (see <Ref Func="ModPCohomologyRing" BookName="HAP"/>,
##  but automatically constructs this list if it is not available.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(ModPRingNiceBasis,
  [IsAlgebra],
  function(A)

  if IsBound(A!.niceBasis) then
    return A!.niceBasis;
  else
    return ModPCohomologyRing_part_2(A)!.niceBasis;
  fi;
end);
#####################################################################


#####################################################################
##  <#GAPDoc Label="ModPRingNiceBasisAsPolynomials_DTmanRings_Dat">
##  <ManSection>
##  <Attr Name="ModPRingNiceBasisAsPolynomials" Arg="A"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  A list which gives the 'nice basis' of the algebra <A>A</A> 
##  (as returned by the second element of <Ref Attr="ModPRingNiceBasis"/>)
##  in terms of products of the indeterminates in the ring presentation 
##  (as given by 
##  <Ref Attr="PresentationOfGradedStructureConstantAlgebra"/>). The
##  <M>i</M>th entry in the list corresponds to the <M>i</M>th basis element
##  returned by <Ref Attr="ModPRingNiceBasis"/>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(ModPRingNiceBasisAsPolynomials,
  [IsAlgebra],
  function(A)

    local gens, ring, polygens, polybasis, Coeff, Bas, i;
    # Which elements are the ring generators?
    gens := ModPRingGenerators(A);
    # Which elements from the basis are these?
    gens := List(gens, i->Position(Basis(A), i));
  
    # What are the generators in terms of polynomial ring elements?
    ring := BaseRing(PresentationOfGradedStructureConstantAlgebra(A));
    polygens := Concatenation([One(ring)], IndeterminatesOfPolynomialRing(ring));
  
    # Put the polynomial ring generators in the correct positions in a list
    polybasis := [];
    for i in [1..Length(polygens)] do
      polybasis[gens[i]] := polygens[i];
    od;
  
    # what is the 'nice basis'?
    Bas := ModPRingNiceBasis(A)[2];
    
    return List(Bas, b->Product(List(b, i->polybasis[i])));
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="ModPRingBasisAsPolynomials_DTmanRings_Dat">
##  <ManSection>
##  <Attr Name="ModPRingBasisAsPolynomials" Arg="A"/>
##
##  <Returns>
##    List
##  </Returns>
##  <Description>
##  A list which gives the basis of the algebra <A>A</A>
##  (as returned by <C>Basis(A)</C>) in terms of sums of products of the 
##  indeterminates in the ring presentation (as given by 
##  <Ref Attr="PresentationOfGradedStructureConstantAlgebra"/>). The
##  <M>i</M>th entry in the list corresponds to the <M>i</M>th basis element
##  returned by <Ref Attr="Basis" BookName="ref"/>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(ModPRingBasisAsPolynomials,
  [IsAlgebra],
  function(A)

    local nicebasis, Coeff;
    
    nicebasis := ModPRingNiceBasisAsPolynomials(A);
    # And what is the standard basis in terms of these?
    Coeff := ModPRingNiceBasis(A)[1];
    
    return List(Coeff, c->c*nicebasis);
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="HAPPRIME_LastLHSBicomplexSize_DTmanGradedAlgebraInt">
##  <ManSection>
##  <Var Name="HAPPRIME_LastLHSBicomplexSize"/>
##  <Description>
##  Stores the last bicomplex size last used by 
##  <Ref Func="LHSSpectralSequenceLastSheet"/>. Used for testing and generating
##  results for paper.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
HAPPRIME_LastLHSBicomplexSize := 0;

#####################################################################
##  <#GAPDoc Label="HAPPRIME_LHSSpectralSequence_DTmanGradedAlgebraInt">
##  <ManSection>
##  <Oper Name="HAPPRIME_LHSSpectralSequence" Arg="G, N, n"/>
##
##  <Returns>
##    <K>GradedAlgebraPresentation</K> or list
##  </Returns>
##  <Description>
##  This function is called by both <Ref Func="LHSSpectralSequence"/> and
##  <Ref Func="LHSSpectralSequenceLastSheet"/> and does the actual work 
##  of computing the Lyndon-Hoschild-Serre spectral sequence for the group 
##  extension <M>N \to G \to G/N</M>. See the documentation for those functions 
##  for the main details, including options that are passed through to this
##  function.
##  <P/>
##  If <A>n</A> is the string <C>"Einf"</C> then only the limiting sheet is 
##  returned (and the other sheets are discarded during computation, which can 
##  save time and memory). Otherwise, <A>n</A> sheets are returned, or
##  enough until convergence is proved, if that is smaller.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_LHSSpectralSequence,
  function(G, N, nsheets)
    local n, GhomQ, Q, RQ, AQ, ringQ, RN, AN, ringN, zeroR, gensG, gensQ, RG, 
    CCG, k, field, message, abelian, d, D, derivations, converged, originalgens,
    EkOriginal, EkCurrent, EkOriginalToCurrentMap, i, H, indetExp, 
    neededlength, temp, sheets, storeallsheets, sheetstodo,
    GroupDescription, GetResolutionAlgebraAndRing, MakeNiceRingPresentation,
    PQDegreeOfMonomial, DegreeOfPolynomial, DerivationImagesOfMonomial,
    DerivationImagesOfPolynomial, PolynomialModIdeal,tmplst,xx;

    if not IsPGroup(G) or not IsFinite(G) then
      Error("<G> must be a p-group");
    fi;
    if not IsFinite(N) or not IsPGroup(N) then
      Error("<N> must be a p-group");
    fi;
    if not IsSubgroup(G, N) then
      Error("<N> must be in the centre of <G>");
    fi;
    if N = G or Order(N) = 1 then
      Error("<N> must be a non-trivial central subgroup of <G>");
    fi;

    # See if the size of the bicomplex is specified. If not, pick a 
    # sensible initial value
    n := ValueOption("InitialLHSBicomplexSize");
    if n = fail then
      n := 5;  # TODO: Have a cleverer way of picking this
    elif not IsPosInt(n) then
      Error("option <InitialLHSBicomplexSize> must be a positive integer");
    fi;
    HAPPRIME_LastLHSBicomplexSize := n;

    # Print a description of the group
    ################################
    GroupDescription := function(G)
      local string;
      string := "";
      if Order(G) < 512 then
        string := Concatenation("SmallGroup(", String(IdSmallGroup(G)[1]), ",", 
            String(IdSmallGroup(G)[2]), "): ", StructureDescription(G));
      else
        string := Concatenation("a group of order ", String(Order(G)));
      fi;
      return string;
    end;
    ################################
    
    # Find the quotient group Q = G/N, and the resolution and
    # ring presentation of this
    GhomQ := NaturalHomomorphismByNormalSubgroup(G, N);;
    Q := Image(GhomQ);;
    Info(InfoHAPprime, 2, "Calculating mod-p cohomology ring presentation for");
    Info(InfoHAPprime, 2, "G = ", GroupDescription(G));
    Info(InfoHAPprime, 2, "Using extension N >--> G -->> Q");
    Info(InfoHAPprime, 2, "N is ", GroupDescription(N));
    Info(InfoHAPprime, 2, "Q is ", GroupDescription(Q));
    Info(InfoHAPprime, 2, "Computing resolutions and presentations...");
    
    #########################################
    # Returns the resolution, cohomology ring and algebra for the group Q,
    # up to at least minn, or more if there are generators in higher 
    # degrees. Also returns the value of n that was used
    GetResolutionAlgebraAndRing := function(Q, minn)
      local ring, n, R, A;
      if ValueOption("NoInductiveProof") = false then
        # Can we get the cohomology ring for Q easily?
        ring := HAPPRIME_CohomologyRingWithoutResolution(Q);
        if ring = fail then
          # Otherwise use the LHS spectral sequence to find the best n
          ring := LHSSpectralSequenceLastSheet(Q);
        fi;
        n := Maximum(Concatenation(IndeterminateDegrees(ring) + 1, [minn]));
#n:=10;#Graham fiddled here
      else
        n := minn;
#n:=10;#Graham fiddled here
      fi;

      R := ResolutionPrimePowerGroup(Q, n);;
      A := ModPCohomologyRing(R);;
      # And now calculate the ring presentation from the algebra
      # (even if we already have it) to guarantee that the generators match up
      ring := Mod2CohomologyRingPresentation(A);;
      return [R, A, ring, n];
    end;
    #########################################

    # Calculate a minimal resolution of Q,
    # and its ring presentation
    RQ := GetResolutionAlgebraAndRing(Q, n);
    Info(InfoHAPprime, 2, "Resolution and presentation of Q found");
    if RQ[4] > n then
      Info(InfoHAPprime, 2, "needed a resolution of length ", RQ[4], " to capture all generators");
    fi;

    # Calculate a minimal resolution of N,
    # and its ring presentation
    RN := GetResolutionAlgebraAndRing(N, RQ[4]);
    Info(InfoHAPprime, 2, "Resolution and presentation of N found");
    
    # Does N need a bigger resolution than Q? If so, we need to
    # redo Q with a bigger resolution. This is unlikely, since Q is
    # typically more complex than N
    if RN[4] > RQ[4] then
      Info(InfoHAPprime, 2, 
        "Redoing resolution and presentation for Q with length of ", RN[4]);
      RQ := GetResolutionAlgebraAndRing(Q, RN[4]);
      Info(InfoHAPprime, 2, "Resolution and presentation of Q found");
    fi;

    # Now extract all the information
    ringQ := RQ[3];
    AQ := RQ[2];
    # Pre-emptively call this attribute here to make sure that One(AQ) is set
    # TODO this won't be needed once HAP is fixed to Set One(A) in 
    # ModPCohomologyRing_partI
    ModPRingBasisAsPolynomials(AQ);
    RQ := RQ[1];
    n := RN[4];
    ringN := RN[3];
    zeroR := Zero(BaseRing(ringN));
    AN := RN[2];
    # Pre-emptively call this attribute here to make sure that One(AN) is set
    # TODO this won't be needed once HAP is fixed
    ModPRingBasisAsPolynomials(AN);
    RN := RN[1];
    # Convert the ring for N to avoid the indeterminants in Q
    ringN := HAPPRIME_GradedAlgebraPresentationAvoidingIndeterminates(
      ringN, IndeterminatesOfGradedAlgebraPresentation(ringQ));
      Info(InfoHAPprime, 2, "Found mapping between presentation and algebra basis");

    # Now calculate the resolution of the extension
    Info(InfoHAPprime, 2, "Building extension bicomplex...");
    gensG := GeneratorsOfGroup(G);;
    gensQ := List(gensG, x->Image(GhomQ,x));;
    RG := ResolutionFiniteExtension(gensG, gensQ, RQ, n, false, RN);;
    # And calculate the cochain of this resolution. This will be used
    # to compute the derivations
    CCG := HomToIntegersModP(RG, 2);
    Info(InfoHAPprime, 2, "Resolution of extension found");

    #########################################
    # Returns the (p,q) location of a monomial in the bicomplex
    PQDegreeOfMonomial := function(mon)

      local unimons, pdeg, qdeg, i, indetExp, k;

      unimons := UnivariateMonomialsOfMonomial(mon);
      pdeg := 0;
      qdeg := 0;
      for i in unimons do
        indetExp := IndeterminateAndExponentOfUnivariateMonomial(i);
        k := Position(IndeterminatesOfGradedAlgebraPresentation(ringQ), indetExp[1]);
        if k <> fail then
          pdeg := pdeg + IndeterminateDegrees(ringQ)[k] * indetExp[2];
        else
          k := Position(IndeterminatesOfGradedAlgebraPresentation(ringN), indetExp[1]);
          qdeg := qdeg + IndeterminateDegrees(ringN)[k] * indetExp[2];
        fi;
      od;

      return [pdeg, qdeg];
    end;
    #########################################

    #########################################
    # Returns the degree of a polynomial in the bicomplex
    DegreeOfPolynomial := function(poly)

      local degree, terms, t, deg;

      degree := [];
      terms := TermsOfPolynomial(poly);
      for t in terms do
        deg := PQDegreeOfMonomial(t[1]);
        if degree = [] then
          degree := deg;
        else
          if degree <> deg then
            Error("terms in polynomial are not all in the same pq-degree");
          fi;
        fi;
      od;

      return Sum(degree);
    end;
    #########################################

    #########################################
    # Returns the image of a monomial (in the presentation of
    # N and Q) under the boundary map as a list of images under
    # the various derivations, in the form [d0, d1, d2, ...]
    # The list is as long as the number of derivations that stay on the 
    # bicomplex - any higher derivation is zero
    DerivationImagesOfMonomial := function(mon)

      local unimons, Nmon, Qmon, Ncoeffs, Qcoeffs, n, Nc, q, NQc, Npair, Qpair, deg,
      gen, dnum, image, p, i, imc, j, vec;
    
      # Split this monomial into Q and N monomials
      unimons := UnivariateMonomialsOfMonomial(mon);
      Nmon := Product(Filtered(unimons, i->i in BaseRing(ringN)));
      Qmon := Product(Filtered(unimons, i->i in BaseRing(ringQ)));
      # Product gives 1 if the list is empty, but its not the 1 that's in the
      # correct ring, so we need to convert to the correct one if that happens
      if IsOne(Nmon) then
        Nmon := One(BaseRing(ringN));
      fi;
      if IsOne(Qmon) then
        Qmon := One(BaseRing(ringQ));
      fi;
      
      # What is this polynomial as multiples of the various generators?
      Ncoeffs := Coefficients(Basis(AN), HAPPRIME_Polynomial2Algebra(AN, ringN, Nmon));
      Qcoeffs := Coefficients(Basis(AQ), HAPPRIME_Polynomial2Algebra(AQ, ringQ, Qmon));
    
      image := [];
      for n in [1..Length(Ncoeffs)] do
        Nc := Ncoeffs[n];
        if IsZero(Nc) then
          continue;
        fi;
        for q in [1..Length(Qcoeffs)] do
          NQc := Nc * Qcoeffs[q];
          if IsZero(NQc) then
            continue;
          fi;
    
          # Find out which [deg, num] this coefficient represents
          Npair := AN!.intToPairModified(n);
          Qpair := AQ!.intToPairModified(q);
          # The total degree is the sum of the two degrees (graded algebra)
          deg := Npair[1] + Qpair[1];
          # Now find the image of it according to the cochain for G
          # First find which generator of G this is
          # vectorToInt takes (p,q,i,j) where p is the Q-degree, q is the N-degree, 
          # and i and j are the generator numbers in the Q and N degrees respectively.
          gen := RG!.vectorToInt(Qpair[1], Npair[1], Qpair[2], Npair[2]);
          # we can now find the boundary of this (in the cochain)
          # The generator number in the cochain matches up with the generator 
          # number in RG, and the degree is the same as the total degree
          imc := CCG!.boundary(deg, gen);
          # im is a vector in deg+1. Divide this up into parts which come from
          # d0 (up), d1(left), d2(DLL), d3(DDLLL) and so on
          for j in [1..Length(imc)] do
            if IsZero(imc[j]) then
              continue;
            fi;
            # Where in the bicomplex is this generator?
            vec := RG!.intToVector(deg+1, j);
            # What is this element as a polynomial?
            p := ModPRingBasisAsPolynomials(AQ)[AQ!.pairToIntModified([vec[1], vec[3]])] *
              ModPRingBasisAsPolynomials(AN)[AN!.pairToIntModified([vec[2], vec[4]])];
            # Look at the Q-degree to see which boundary map dn this corresponds to
            # We have gone from location Qpair[1] to location vec[1]
            dnum := vec[1] - Qpair[1];
            # And fill in the image. If the image hasn't been prepared, do this
            if Length(image) < (dnum + 1) then
              for i in [(Length(image)+1)..(dnum+1)] do
                image[i] := zeroR;
              od;
            fi;
            image[dnum+1] := image[dnum+1] + NQc * imc[j]*p;
          od;
        od;
      od;
      return image;
    end;
        
    #########################################
    # Returns the image of a monomial (in the presentation of
    # N and Q) under the boundary map as a list of images under
    # the various derivations, in the form [d0, d1, d2, ...]
    # The list is as long as the number of derivations that stay on the 
    # bicomplex - any higher derivation is zero
    DerivationImagesOfPolynomial := function(poly)
      local terms;

      # If poly is all in ringQ then all derivations are always zero
      if poly in BaseRing(ringQ)  then
        return [];
      fi;
      
      # Check that the resolution is long enough for this degree
#      if(DegreeOfPolynomial(poly)+1) > ResolutionLength(RG) then
if(DegreeOfPolynomial(poly)+1) > Length(RG) then
##Modified by Graham Ellis

        return fail;
      fi;
      # Calculate the image from the images of the terms
      terms := TermsOfPolynomial(poly);
      return Sum(List(terms, t->DerivationImagesOfMonomial(t[1])*t[2]));
    end;
    #########################################

    # Now we want to build up the spectral sequence rings E_k
    # and work out the derivations of the generators.

    k := 2;
    # Store the field of the ring seperately
    field := CoefficientsRing(ringQ);

    # Create our current estimage for the ring of G, but first we'll
    # print out where we're starting from
    if InfoLevel(InfoHAPprime) >= 1 then
      message := Concatenation("E_2 = ", String(field),
        String(IndeterminatesOfGradedAlgebraPresentation(ringQ)));
        if not IsEmpty(GeneratorsOfPresentationIdeal(ringQ)) then
          message := Concatenation(message, "/", 
            String(GeneratorsOfPresentationIdeal(ringQ)));
      fi;
      message := Concatenation(message, " x ", String(field),
        String(IndeterminatesOfGradedAlgebraPresentation(ringN)));
      if not IsEmpty(GeneratorsOfPresentationIdeal(ringN)) then
        message := Concatenation(message, "/", 
          String(GeneratorsOfPresentationIdeal(ringN)));
      fi;
      Info(InfoHAPprime, 1, message);
      Info(InfoHAPprime, 1,
        "  with generator degrees ", IndeterminateDegrees(ringQ), 
        " and ", IndeterminateDegrees(ringN), " respectively");
    fi;

    # This is Ek given in terms of the original ring indeterminates
    # We store EkOriginal in the form [ring, gens, degrees, rels]
    originalgens := Concatenation(
      IndeterminatesOfGradedAlgebraPresentation(ringQ),
      IndeterminatesOfGradedAlgebraPresentation(ringN));
    EkOriginal := rec(
      ring := PolynomialRing(field, originalgens),
      gens := originalgens,
      rels := Concatenation(GeneratorsOfPresentationIdeal(ringQ), 
        GeneratorsOfPresentationIdeal(ringN)),
      degrees := Concatenation(IndeterminateDegrees(ringQ),
        IndeterminateDegrees(ringN)));
    # Make sure that the relations are a Groebner basis
    SetTermOrdering(EkOriginal.ring, "dp");
    if not IsEmpty(EkOriginal.rels) then    
      EkOriginal.rels := SingularReducedGroebnerBasis(
        Ideal(EkOriginal.ring, EkOriginal.rels));
    fi;
    # We also store Ek given in terms of the current ring indeterminates
    # in the same format
    EkCurrent := rec(ring := EkOriginal.ring, gens := EkOriginal.gens, 
      rels := EkOriginal.rels, degrees := EkOriginal.degrees);

    # And the map between the two presentations
    EkOriginalToCurrentMap := false;


    #########################################
    # Returns the polynomial <poly> (in the original ring) modulo the 
    # current ideal, as a list with three elements: first <poly> in the 
    # current ring, then the result of modding out I in the current 
    # ring and the result in the original ring
    PolynomialModIdeal := function(poly)
      local p;

      if IsZero(poly) then
        return [poly, poly, poly];
      fi;

      if EkOriginalToCurrentMap <> false then
        # Convert the result into the current ring indeterminates
        p := ImageOfRingHomomorphism(EkOriginalToCurrentMap, poly);
        # Is it zero in the current ring indeterminates?
        if IsZero(p) then
          return [p,p,p];
        else
          # Reduce this by the ideal 
          if not IsEmpty(EkCurrent.rels) then
            p := [p, SingularPolynomialNormalForm(
              p, Ideal(EkCurrent.ring, EkCurrent.rels))];
          else
            p := [p, p];
          fi;
          # And convert back for our message as we pass it back
          p[3] := ImageOfRingHomomorphism(
            InverseRingHomomorphism(EkOriginalToCurrentMap), p[2]);
          return p;
        fi;
      else
        # Reduce this by the ideal 
        if not IsEmpty(EkCurrent.rels) then
          p := SingularPolynomialNormalForm(poly, 
            Ideal(EkOriginal.ring, EkOriginal.rels));
        else
          p := poly;
        fi;
        return [poly, p, p];
      fi;
    end;
    #########################################

    #################################################
    ## Convert EkCurrent/EkOriginal into a nice presentation with
    ## indeterminates starting at 1 and sorted degrees
    ## if display is true then some Info messages are printed
    MakeNiceRingPresentation := function(display)
      local ring, degrees, tempring, indets, rels, k, message, i,
      CurrentToTempring, TempringToRing, CurrentToRing;

      if EkOriginalToCurrentMap <> false then
        if display then
          Info(InfoHAPprime, 1, "Renaming indeterminates and sorting into increasing degree");
        fi;
        # Rephrase the result in a new ring with indeterminants 
        # [x_1, x_2, ... x_n]
        ring := PolynomialRing(field, Length(EkCurrent.gens));
        degrees := ShallowCopy(EkCurrent.degrees);

        # And create a temporary ring that avoids these, and avoids EkCurrent
        tempring := PolynomialRing(field,
          Length(IndeterminatesOfPolynomialRing(ring)),
          Concatenation(IndeterminatesOfPolynomialRing(ring), EkCurrent.gens));
          # Swap to the temporary ring and find the relations
        CurrentToTempring := HAPRingHomomorphismByIndeterminateMap(
          EkCurrent.ring, EkCurrent.rels, tempring);
        rels := ImageRelations(CurrentToTempring);
        # And now swap back to our final ring, but sort them
        # so that the degrees are of increasing order
        indets := ShallowCopy(IndeterminatesOfPolynomialRing(tempring));
        SortParallel(degrees, indets);
        TempringToRing := HAPRingHomomorphismByIndeterminateMap(
          PolynomialRing(field, indets), rels, ring);
        rels := ShallowCopy(ImageRelations(TempringToRing));

        # Print out the change of variables 
        if display and InfoLevel(InfoHAPprime) >= 2 then
          CurrentToRing := CompositionRingHomomorphism(
            CurrentToTempring, TempringToRing);
          message := "Change of indeterminates: ";
          for i in [1..Length(SourceGenerators(CurrentToRing))] do 
            message := Concatenation(
              message, String(ImageGenerators(CurrentToRing)[i]), "=", 
              String(SourceGenerators(CurrentToRing)[i]), ", ");
          od;
          Remove(message, Length(message)); # Remove the last space
          Remove(message, Length(message)); # Remove the last comma
          Info(InfoHAPprime, 2, message);
        fi;

        # And sort the relations to make them pretty
        Sort(rels);

        return GradedAlgebraPresentation(ring, rels, degrees);
      else
        rels := ShallowCopy(EkOriginal.rels);
        Sort(rels);
        return GradedAlgebraPresentation(
          PolynomialRing(field, EkOriginal.gens), rels, 
            EkOriginal.degrees);
      fi;
    end;
    #################################################

    if nsheets = "Einf" then
      sheetstodo := infinity;
      storeallsheets := false;
    else
      sheetstodo := nsheets;
      storeallsheets := true;
      sheets := [];
      sheets[k] := MakeNiceRingPresentation(false);
    fi;

    converged := false;
    while k < sheetstodo and not converged do
      # Find the derivations of all my generators
      converged := true;
      neededlength := [];
      derivations := [];
      for i in EkOriginal.gens do

        message := Concatenation(
          "  d_", String(k), "(", String(i), ")", " = ");

          if IsZero(i) then 
          message := Concatenation(message, String(zeroR));
          Add(derivations, zeroR);

          continue;
          fi;


        D := DerivationImagesOfPolynomial(i);
        if D = fail then
          # We don't have a long enough resolution.
          # What is the length of resolution that we need?
          Add(neededlength, DegreeOfPolynomial(i) + 1);
          message := Concatenation(message, 
            "need resolution of length ", String(neededlength[Length(neededlength)]));
        else
          if not IsEmpty(D) then    
            # Check that all the images of degrees less than k are zero
            for i in [1..Minimum(k, Length(D))] do
              if not IsZero(PolynomialModIdeal(D[i])[2]) then
                Error("all derivations less than k are not zero");
              fi;
            od;
          fi;
  
          # Now calculate the derivation d_k
          if Length(D) < (k+1) then
            message := Concatenation(message, "zero");
            Add(derivations, zeroR);
          else
            message := Concatenation(message, String(D[k+1]));
  
            # Now see what this ends up as in the new current indeterminants
            # and modulo the ideal
            d := PolynomialModIdeal(D[k+1]);
            if IsZero(d[1]) and not IsZero(D[k+1]) then
              # it is zero in the new indeterminates
              message := Concatenation(message, " = ", String(d[1]));
              Add(derivations, zeroR);
            else
              if d[3] <> D[k+1] then
                # It is something different modulo the ideal
                message := Concatenation(message, " = ", String(d[3]) , " mod I");
              fi;
              Add(derivations, d[2]);
            fi;
            converged := false;
          fi;
        fi;
        Info(InfoHAPprime, 1, message);
      od;

      # Check that the length of the resolution is sufficient. If not, then
      # needlength is not empty
      if not IsEmpty(neededlength) then
        if ValueOption("LargerLHSBicomplexBreak") = true then
          Error("the bicomplex is too small for convergence: you will need at least ", 
            Maximum(neededlength), " terms. Type return; to continue with a bicomplex of this size.\n");
        fi;
        if ValueOption("LargerLHSBicomplexFail") = true then
          Info(InfoHAPprime, 1, 
            " **** Not enough terms in the resolution");
          Info(InfoHAPprime, 1, 
            " **** Needs a resolution of length ", Maximum(neededlength));
          return fail;
        fi;
        # Call ourselves again with a longer length
        Info(InfoHAPprime, 1, 
          " **** Not enough terms in the resolution");
        Info(InfoHAPprime, 1, 
          " **** Repeating with resolutions of length ", Maximum(neededlength));
        return HAPPRIME_LHSSpectralSequence(
          G, N, nsheets : InitialLHSBicomplexSize := Maximum(neededlength));
      fi;

      if not converged then
        # make the derivation object and find its homology
        d := HAPDerivation(EkCurrent.ring, EkCurrent.rels, derivations);

        if EkOriginal.gens = EkCurrent.gens then
          H := HomologyOfDerivation(d);
        else
          H := HomologyOfDerivation(d, originalgens);
        fi;

        # now construct the new E_{k+1} in both the old and the new indeterminates

        EkCurrent.ring := H[1];
        EkCurrent.gens := IndeterminatesOfPolynomialRing(H[1]);
        # Make sure that the relations are a Groebner basis
        SetTermOrdering(EkCurrent.ring, "dp");
        EkCurrent.rels := SingularReducedGroebnerBasis(
          Ideal(EkCurrent.ring, H[2]));

        # The last element of H gives the mapping between the homology in terms
        # of the original ring indeterminants and the previous ones
        # We want a mapping between the original ring and the new one, so
        # combine this with the previous one.
        if EkOriginalToCurrentMap = false then
          EkOriginalToCurrentMap := H[3];
        else
          # H[3] is the map from whose image is the current ring. Take the
          # generators of the current ring (or, rather, the source of H[3] which
          # is the same generators in terms of the old current), and map them
          # back to the original current.
          # Also map the source relations of H[3] (which may have relations
          # which are lost at the image of H[3]), and map those all the way
          # back to the orginal ring. Yet more may have been lost earlier, so
          # we always concatenate these with the original relations.

#####
#####Graham  added this to test things!
tmplst:=List( SourceGenerators(H[3]),
              i->PreimageOfRingHomomorphism(EkOriginalToCurrentMap, i));
for xx in [1..Length(tmplst)] do
if IsZero(tmplst[xx]) then tmplst[xx]:=EkOriginal.rels[1]; fi;
od;
#####This change  has a net effect!
#####

          EkOriginalToCurrentMap := HAPSubringToRingHomomorphism(
          tmplst,
            Concatenation(SourceRelations(EkOriginalToCurrentMap),
              Filtered(List(SourceRelations(H[3]), 
                i->PreimageOfRingHomomorphism(EkOriginalToCurrentMap, i)), 
                  j->not IsZero(j))),

            ImagePolynomialRing(H[3]));
        fi;
        # And provide the current generators and relations in terms of the 
        # original ones, so that we can print things out nicely
        EkOriginal.gens := SourceGenerators(EkOriginalToCurrentMap);


        EkOriginal.rels := SourceRelations(EkOriginalToCurrentMap);
        EkOriginal.degrees := List(EkOriginal.gens, i->DegreeOfPolynomial(i));
        # And remember the current degrees as well
        EkCurrent.degrees := EkOriginal.degrees;

        k := k+1;
        if storeallsheets then
          sheets[k] := MakeNiceRingPresentation(false);
        fi;
      else
        k := "inf";
      fi;

      # And display what the new Ek is
      if InfoLevel(InfoHAPprime) >= 1 then
        message := Concatenation("E_", String(k), " = ");
        message := Concatenation(message, String(field), String(EkCurrent.gens));
        if not IsEmpty(EkCurrent.rels) then
          message := Concatenation(message, "/", String(EkCurrent.rels));
        fi;
        Info(InfoHAPprime, 1, message);
        if EkOriginalToCurrentMap <> false then
          message := Concatenation("E_", String(k), " = ");
          message := Concatenation(message, String(field), String(EkOriginal.gens));
          if not IsEmpty(EkOriginal.rels) then
            message := Concatenation(message, "/", String(EkOriginal.rels));
          fi;
          Info(InfoHAPprime, 1, message);
        fi;  
      fi;
    od;
    
    if storeallsheets then
      return sheets;
    else
      return MakeNiceRingPresentation(true);
    fi;

  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="HAPPRIME_CohomologyRingWithoutResolution_DTmanGradedAlgebraInt">
##  <ManSection>
##  <Func Name="HAPPRIME_CohomologyRingWithoutResolution" Arg="G"/>
##
##  <Returns>
##    <K>GradedAlgebraPresentation</K> or <K>fail</K>
##  </Returns>
##  <Description>
##  If the mod-&p; cohomology ring for <A>G</A> can be computed without 
##  building a resolution for <A>G</A> then this ring is returned, otherwise
##  this function returns <K>fail</K>.
##  <P/>
##  Current cases where the resolution for <A>G</A> is not needed are 
##  <List>
##    <Item>if <A>G</A> is the group of order two then 
##      <M>H^*(G, \mathbb{F}) = \mathbb{F}[x]</M> where <M>x</M> has degree 
##      one</Item>
##    <Item>if <A>G</A> is cyclic of order greater than two then 
##      <M>H^*(G, \mathbb{F}) = \mathbb{F}[x,y]/x^2</M> where <M>x</M> has 
##      degree one and <M>y</M> has degree two</Item>
##    <Item>if <A>G</A> can be expressed as the direct sum of other groups
##      then the cohomology rings for those groups are found (by recursive
##      calls to <Ref Func="Mod2CohomologyRingPresentation" BookName="HAPprime"/>)
##      and combined using 
##      <Ref Func="TensorProductOp" Label="for collection of algebra presentations"/></Item>
##  </List>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_CohomologyRingWithoutResolution,
  function(G)
    local R, z, DF, C;
    
    if IsCyclic(G) then
      if Order(G) = 2 then
        # The cohomology ring for C2 is k[x] with x having degree 1
        R := PolynomialRing(GF(2), 1);
       return GradedAlgebraPresentation(R, [], [1]);
      else
        # The cohomology ring for Cn (for n <> 2) is
        # k[zy]/z^2 with z having degree 1 and y degree 2
        R := PolynomialRing(GF(2), 2);
        z := IndeterminatesOfPolynomialRing(R)[1];
        return GradedAlgebraPresentation(R, [z^2], [1, 2]);
      fi;
    fi;
#return fail;   ################    
    DF := DirectFactorsOfGroup(G);
    if Length(DF) > 1 then
      Info(InfoHAPprime, 2, "Calculating mod-p cohomology ring presentation for");
      Info(InfoHAPprime, 2, "G = SmallGroup(", IdSmallGroup(G)[1], ",", 
        IdSmallGroup(G)[2], "): ", StructureDescription(G));
      Info(InfoHAPprime, 2, "by splitting group into ", Length(DF), " direct factors");
      C := List(DF, Mod2CohomologyRingPresentation);
      # if LargerLHSBicomplexFail is set, we may have some fails in the list
      if not IsHomogeneousList(C) then
        return fail;
      else
        return TensorProductOp(C);
      fi;
    fi;
    
    return fail;
  end);
#####################################################################

#####################################################################
##  <#GAPDoc Label="HAPPRIME_Polynomial2Algebra_DTmanGradedAlgebraInt">
##  <ManSection>
##  <Func Name="HAPPRIME_Polynomial2Algebra" Arg="A[, ringA], poly"/>
##
##  <Returns>
##    Algebra element
##  </Returns>
##  <Description>
##  Converts a polynomial <A>poly</A> in <A>ringA</A> (the mod-&p; cohomology 
##  ring presentation for the cohomology algebra <A>A</A>)
##  into the equivalent element in <A>A</A>.
##  If <A>ringA</A> is not provided, it is recovered using the attribute
##  <Ref Attr="PresentationOfGradedStructureConstantAlgebra"/>. The 
##  <A>ringA</A> argument is provided to allow cases where the ring and 
##  polynomial use different (but isomorphic) indeterminants from those
##  provided by <Ref Attr="PresentationOfGradedStructureConstantAlgebra"/> 
##  (as is the case in <Ref Func="HAPPRIME_LHSSpectralSequence"/>).
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_Polynomial2Algebra,
  function(arg)

    local A, ringA, poly, gens, ring, indets, terms, sum, term, prod, mons, 
      umon, m;
    
    if Length(arg) = 2 then
      A := arg[1];
      ringA := Mod2CohomologyRingPresentation(A);;
      poly := arg[2];
    elif Length(arg) = 3 then
      A := arg[1];
      ringA := arg[2];
      poly := arg[3];
    else
      Error("there must be either two or three argments");
    fi;

    # Which elements are the ring generators?
    gens := ModPRingGenerators(A);

    # The first one of these corresponds to "1", and the rest to the ring
    # indeterminates in order. So, get these
    ring := BaseRing(ringA);
    if not poly in ring then
      Error("<poly> must be a polynomial in the ring presentation of <A>");
    fi;
    indets := IndeterminatesOfPolynomialRing(ring);

    # Now decompose our polynomial
    terms := TermsOfPolynomial(poly);
    sum := Zero(A);
    for term in terms do
      # Each term is a list of pairs [mon, coeff]. Take each of these in turn
      prod := One(A);
      mons := UnivariateMonomialsOfMonomial(term[1]);
      # mon is a list of univariate monomials that make up the total monomial
      for umon in mons do
        m := IndeterminateAndExponentOfUnivariateMonomial(umon);
        if m[2] <> 0 then
          # Each m is a pair [indet, exponent].
          # What is this indet in terms of ring generators?
          prod := prod * gens[Position(indets, m[1]) + 1]^m[2];
        else
          # except if m could be [coeff, 0], which indicates a constant
          prod := m[1] * gens[1];
        fi;
      od;
      sum := sum + term[2]*prod;
    od;

    return sum;
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="HAPPRIME_Algebra2Polynomial_DTmanGradedAlgebraInt">
##  <ManSection>
##  <Func Name="HAPPRIME_Algebra2Polynomial" Arg="A, alg"/>
##
##  <Returns>
##    Polynomial
##  </Returns>
##  <Description>
##  Converts an algebra element <A>alg</A> in the mod-&p; cohomology algebra
##  <A>A</A> (given by <Ref Attr="ModPCohomologyRing" BookName="HAP"/>)
##  into the equivalent polynomial from the corresponding ring presentation
##  (given by <Ref Attr="PresentationOfGradedStructureConstantAlgebra"/>).
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_Algebra2Polynomial,
  function(A, alg)

    return Coefficients(Basis(A), alg) * ModPRingBasisAsPolynomials(A);
  end
);
#####################################################################


