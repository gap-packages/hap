#############################################################################
##
##  HAPPRIME - singular.gi
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


#####################################################################
##  <#GAPDoc Label="SingularSetNormalFormIdeal_manSingular">
##  <ManSection>
##  <Oper Name="SingularSetNormalFormIdeal" Arg="I"/>
##  <Oper Name="SingularSetNormalFormIdealNC" Arg="I"/>
##  <Returns>
##    nothing
##  </Returns>
##  <Description>
##  Sets the ideal to be used by singular for any subsequent calls to 
##  <Ref Func="SingularPolynomialNormalForm"/> to be <A>I</A>. After calling 
##  this function, the singular base ring and term ordering (see
##  <Ref Func="SingularBaseRing" BookName="singular"/> and 
##  <Ref Func="TermOrdering" BookName="singular"/>) will be set to be that of 
##  the ring containing <A>I</A>, so an additional call to 
##  <Ref Func="SingularSetBaseRing" BookName="singular"/> is not necessary.
##  <P/> 
##  The standard form of this function ensures that <A>I</A> is 
##  a reduced Gröbner basis with respect to the value of 
##  <Ref Func="TermOrdering" BookName="singular"/> for the ring containing the 
##  ideal, while the <C>NC</C> assumes that <A>I</A> is already such a Gröbner 
##  basis.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(SingularSetNormalFormIdeal,
  [IsPolynomialRingIdeal],
  function(I)
    SingularSetNormalFormIdealNC(
      Ideal(LeftActingRingOfIdeal(I), SingularReducedGroebnerBasis(I)));
  end
);
#####################################################################
InstallMethod(SingularSetNormalFormIdealNC,
  [IsPolynomialRingIdeal],
  function(I)
    local input, out;

    # Set the base ring if it is not the same or it doesn't have 
    # the same indeterminate order
    if LeftActingRingOfIdeal(I) <> SingularBaseRing or
      IndeterminatesOfPolynomialRing(LeftActingRingOfIdeal(I)) <> 
        IndeterminatesOfPolynomialRing(SingularBaseRing) then
      SingularSetBaseRing(LeftActingRingOfIdeal(I));
    fi;

    Info( InfoSingular, 2, "setting GAP_ideal ideal..." );

    # preparing the input for Singular
    input := "";

    if not HasGeneratorsOfTwoSidedIdeal(I) then
      # An ideal has no generators if the list of relations is empty.
      input := "ideal GAP_NFideal = ideal();\n";
    else
      Append( input, "ideal GAP_NFideal = " );
      Append( input, ParseGapIdealToSingIdeal( I ) );
      Append( input, ";\n" );
    fi;
    
    out := SingularCommand( input, "" );

    Info( InfoSingular, 2, "done SingularSetIdeal." );
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="SingularPolynomialNormalForm_manSingular">
##  <ManSection>
##  <Oper Name="SingularPolynomialNormalForm" Arg="poly[, I]"/>
##  <Returns>
##    Polynomial
##  </Returns>
##  <Description>
##  Returns the normal form of the polynomial <A>poly</A> after reduction
##  by the ideal <A>I</A>. The ideal can either be passed to this function,
##  in which case it is converted to a Gröbner basis (with respect to the 
##  term ordering of the ideal's ring - see 
##  <Ref Func="TermOrdering" BookName="singular"/>), or the ideal to use can
##  be set first be calling <Ref Func="SingularSetNormalFormIdeal"/>, which 
##  is more efficient for repeated use of this function (the latter function 
##  also sets the base ring and term ordering).
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(SingularPolynomialNormalForm,
  [IsPolynomial, IsPolynomialRingIdeal],
  function(poly, I)
    SingularSetNormalFormIdeal(I);
    return SingularPolynomialNormalForm(poly);
  end
);
#####################################################################
InstallOtherMethod(SingularPolynomialNormalForm,
  [IsPolynomial],
  function(poly)
    local input, out;

    Info( InfoSingular, 2, "reducing polynomial to normal form..." );

    # preparing the input for Singular
    input := "";

    Append( input, "poly GAP_poly = reduce( " );
    Append( input, ParseGapPolyToSingPoly( poly ) );
    Append( input, ", GAP_NFideal);\n" );

    out := SingularCommand( input, "string (GAP_poly)" );

    Info( InfoSingular, 2, "done SingularPolynomialNF." );

    # Fix for singular ordering bug
    return PolynomialByExtRep(FamilyObj(poly), 
      ExtRepPolynomialRatFun(ParseSingPolyToGapPoly(out)));
  end
);
#####################################################################



#####################################################################
##  <#GAPDoc Label="SingularGroebnerBasis_manSingular">
##  <ManSection>
##  <Attr Name="SingularGroebnerBasis" Arg="I"/>
##  <Returns>
##    List 
##  </Returns>
##  <Description>
##    Returns a list of relations which form a Gröbner basis for the ideal 
##    <A>I</A> given the <Ref Attr="TermOrdering" BookName="singular"/>
##    associated with the ring containing <A>I</A>. This function is the 
##    same as the <Package>singular</Package> function
##    <Ref Meth="GroebnerBasis" BookName="singular"/>, but fixes a bug in
##    that package when using unusual term ordering.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(SingularGroebnerBasis,
  [IsPolynomialRingIdeal], 
  function(I)
    local rels, fam;
    
    # An ideal has no generators if the list of relations is empty.
    # If so, return the empty list (singular doesn't check for this)
    if not HasGeneratorsOfTwoSidedIdeal(I) then
      return [];
    fi;
    
    rels := GroebnerBasis(I);
    if IsEmpty(rels) then
      return rels;
    fi;
    fam := FamilyObj(rels[1]);
    # Fix for singular ordering bug
    return List(rels, i->PolynomialByExtRep(fam, ExtRepPolynomialRatFun(i)));
  end
);
#####################################################################

#####################################################################
##  <#GAPDoc Label="SingularReducedGroebnerBasis_manSingular">
##  <ManSection>
##  <Attr Name="SingularReducedGroebnerBasis" Arg="I"/>
##  <Returns>
##    List 
##  </Returns>
##  <Description>
##    Returns a list of relations which form a reduced Gröbner basis for the 
##    ideal <A>I</A> given the <Ref Attr="TermOrdering" BookName="singular"/>
##    associated with the ring containing <A>I</A>. This function is the 
##    equivalent of the <Package>singular</Package> function
##    <Ref Meth="GroebnerBasis" BookName="singular"/> (and uses that function), 
##    but ensures that a reduced basis is returned.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(SingularReducedGroebnerBasis,
  [IsPolynomialRingIdeal], 
  function(I)
    local rels;

    # Remember the current singular options
    SingularCommand( "", "intvec GAP_optionsstore = option(get);");
    # set redSB to ask for reduced a Groebner basis
    SingularCommand( "", "option(redSB);");

    rels := SingularGroebnerBasis(I);

    # Set the options back to where they were
    SingularCommand( "", "option(set, GAP_optionsstore);");
    
    return rels;
  end 
);
#####################################################################
    
