#############################################################################
##
##  HAPPRIME - derivation.gi
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
##  <#GAPDoc Label="HAPDerivationRep_DTmanDerivationNODOC">
##  <ManSection>
##  <Filt Name="IsHAPDerivationRep" Arg="O" Type="Representation"/>
##  <Description>
##  Returns <K>true</K> if the object is in the internal representation used for 
##  a <K>HAPDerivation</K>, or <K>false</K> otherwise
##  </Returns>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareRepresentation(
  "IsHAPDerivationRep",
  IsComponentObjectRep and IsAttributeStoringRep,
  ["ring", "relations", "images"]
);
# Note this also defines the function IsHAPDerivationRep
#####################################################################

#####################################################################
# The type for a HAPDerivation is a HAPDerivation in 
# HAPDerivationRep representation, in the HAPDerivation family
HAPDerivationType := 
  NewType(NewFamily("HAPDerivationFamily"), IsHAPDerivation and IsHAPDerivationRep);
#####################################################################



#####################################################################
##  <#GAPDoc Label="HAPDerivation_DTmanDerivation_Con">
##  <ManSection Label="HAPDerivationConstructors">
##  <Heading>HAPDerivation construction functions</Heading>
##  <Oper Name="HAPDerivation" Arg="R[, I],  images"/>
##  <Oper Name="HAPDerivationNC" Arg="R, I, images"/>
##
##  <Returns>
##  <K>HAPDerivation</K>
##  </Returns>
##  <Description>
##  Construct a <K>HAPDerivation</K> object representing the derivation
##  <M>d</M> where <A>R</A> is a polynomial ring and 
##  <A>images</A> is the list of the images of the ring indeterminates under the 
##  derivation, <M>\{d(x_1), d(x_2), \ldots, d(x_n)\}</M>. 
##  An optional set of relations <A>I</A> can also be provided, which are
##  passed to <Ref Oper="KernelOfDerivation"/> when calculating the kernel
##  or homology of this derivation.
##  The function <K>HAPDerivation</K> checks that the arguments are
##  compatible, while the <C>NC</C> method performs no checks.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(HAPDerivation,
  [IsPolynomialRing, IsHomogeneousList, IsHomogeneousList and IsRationalFunctionCollection],
  function(ring, relations, images)
    local i;
    # Check that there is one image for each ring indeterminant
    if Length(IndeterminatesOfPolynomialRing(ring)) <> Length(images) then
      Error("the number of images in <images> is not the same as the number of indeterminants in <ring>.");
    fi;
    # Check that each element in <relations> is in the ring
    for i in relations do
      if not i in ring then
        Error("the entries in <relations> must all be polynomials from <ring>.");
      fi;
    od;
    # Check that each element in <images> is in the ring
    for i in images do
      if not i in ring then
        Error("the entries in <images> must all be polynomials from <ring>.");
      fi;
    od;
    return HAPDerivationNC(ring, relations, images);
  end
);
#####################################################################
InstallOtherMethod(HAPDerivation,
  [IsPolynomialRing, IsHomogeneousList and IsRationalFunctionCollection],
  function(ring, images)
    return HAPDerivation(ring, [], images);
  end
  );
#####################################################################
InstallMethod(HAPDerivationNC,
  [IsPolynomialRing, IsHomogeneousList, IsHomogeneousList and IsRationalFunctionCollection],
  function(ring, relations, images)
    local images2, indets, i;
    # The ith entry in image corresponds to the ith entry in 
    # IndeterminatesOfPolynomialRing(ring), and not necessarily to
    # the indeterminate _number_. We want to store them according to
    # indeterminant number to make later operations more efficient.
    images2 := [];
    indets := IndeterminatesOfPolynomialRing(ring);
    for i in [1..Length(images)] do
      images2[IndeterminateNumberOfUnivariateRationalFunction(indets[i])] := images[i];
    od;
    return Objectify( 
      HAPDerivationType, 
      rec(ring := ring, relations := relations, images := images2)
      );
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="DerivationRing_DTmanDerivation_Dat">
##  <ManSection>
##  <Attr Name="DerivationRing" Arg="d"/>
##
##  <Returns>
##  Polynomial ring
##  </Returns>
##  <Description>
##  Returns the ring over which the derivation <A>d</A> is defined. 
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  DerivationRing,
  [IsHAPDerivation],
  function(d)
    return d!.ring;
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="DerivationImages_DTmanDerivation_Dat">
##  <ManSection>
##  <Attr Name="DerivationImages" Arg="d"/>
##
##  <Returns>
##  List of polynomials
##  </Returns>
##  <Description>
##  A derivation <A>d</A> over a (quotient) polynomial ring is defined by 
##  a set of images. 
##  This function returns this list of images. The <M>i</M>th element of the 
##  list is the image of indeterminate number <M>i</M> in that ring family.
##  (Note that indeterminate number <M>i</M> is not necessarily the <M>i</M>th
##  indeterminate in that particular ring. 
##  See <Ref Sect="Indeterminates" BookName="ref"/> for more details.)
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  DerivationImages,
  [IsHAPDerivation],
  function(d)
    return d!.images;
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="DerivationRelations_DTmanDerivation_Dat">
##  <ManSection>
##  <Attr Name="DerivationRelations" Arg="d"/>
##
##  <Returns>
##  List of polynomials
##  </Returns>
##  <Description>
##  Returns the relations of the quotient ring over which the derivation 
##  <A>d</A> is defined. 
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  DerivationRelations,
  [IsHAPDerivation],
  function(d)
    return d!.relations;
  end
  );
#####################################################################
  
  
#####################################################################
##  <#GAPDoc Label="ViewObj_DTmanDerivationNODOC">
##  <ManSection>
##  <Meth Name="ViewObj" Arg="d" Label="for HAPDerivation"/>
##
##  <Description>
##  Prints a short description of the derivation <A>d</A>. This is the usual 
##  description printed by &GAP;.
##  </Description>
##  </ManSection>
##  <Log><![CDATA[
##  gap> View(d);
##  <Derivation>
##  ]]></Log>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  ViewObj,
  "for HAPDerivation",
  [IsHAPDerivation],
  function(obj)
    Print("<Derivation>");
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="PrintObj_DTmanDerivationNODOC">
##  <ManSection>
##  <Meth Name="PrintObj" Arg="d" Label="for HAPDerivation"/>
##
##  <Description>
##  Prints a detailed description of the derivation <A>d</A>.
##  </Description>
##  </ManSection>
##  <Log><![CDATA[
##  gap> Print(d);
##  Derivation over PolynomialRing( GF(2), ["x_1", "x_2"] ), with images:
##  d(x_1) = Z(2)^0
##  d(x_2) = Z(2)^0
##  ]]></Log>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  PrintObj,
  "for HAPDerivation",
  [IsHAPDerivation],
  function(obj)
    Print("HAPDerivation(", DerivationRing(obj), ", ", 
    DerivationRelations(obj), ", ", 
    DerivationImages(obj){
      List(IndeterminatesOfPolynomialRing(DerivationRing(obj)), 
        IndeterminateNumberOfUnivariateRationalFunction)}, ")");
  end
);
#####################################################################


#####################################################################
##  <#GAPDoc Label="Display_DTmanDerivationNODOC">
##  <ManSection>
##  <Meth Name="Display" Arg="d" Label="for HAPDerivation"/>
##
##  <Description>
##  Displays the derivation <A>d</A> in a human-readable form. 
##  </Description>
##  </ManSection>
##  <Log><![CDATA[
##  gap> Display(d);
##  Derivation over PolynomialRing( GF(2), ["x_1", "x_2"] )
##  , with images:
##  d(x_1) = Z(2)^0
##  d(x_2) = Z(2)^0
##  ]]></Log>
##  <#/GAPDoc>
#####################################################################
InstallMethod(
  Display,
  "for HAPDerivation",
  [IsHAPDerivation],
  function(obj)
    local indets, i, images;
    Print("Derivation over ");
    Display(DerivationRing(obj));
    if not IsEmpty(DerivationRelations(obj)) then
      Print(" / ", DerivationRelations(obj));
    fi;
    Print(", with images:\n");
    indets := IndeterminatesOfPolynomialRing(DerivationRing(obj));
    images := DerivationImages(obj);
    for i in [1..Length(indets)] do
      Print("d(", indets[i], ") = ");
      Display(images[IndeterminateNumberOfUnivariateRationalFunction(indets[i])]);
    od; 
  end
  );
#####################################################################



#####################################################################
##  <#GAPDoc Label="ImageOfDerivation_DTmanDerivation_Hom">
##  <ManSection>
##  <Oper Name="ImageOfDerivation" Arg="d, poly"/>
##
##  <Returns>
##  Polynomial
##  </Returns>
##  <Description>
##  Returns the image of the polynomial <A>poly</A> under the derivation
##  <A>d</A>. (<A>poly</A> must be a polynomial in the
##  derivation's ring.)
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(ImageOfDerivation, 
  [IsHAPDerivation, IsPolynomial],
  function(d, poly)
    local ring, images, oneR, zeroR, t, image, ImageOfMonomial;

    ring := DerivationRing(d);
    zeroR := Zero(ring);

    if not poly in ring then
      Error("<poly> must be a polynomial in the ring of <d>");
    fi;
    if IsOne(poly) or IsZero(poly) then
      return zeroR;
    fi;

    images := DerivationImages(d);
    oneR := One(ring);


    ##################################
    # Return (as a polynomial) the image of the 
    # monomial 
    ImageOfMonomial := function(mon)

      local coeff, n, i, j, image, im, ImageOfPart;

      if IsZero(mon) or IsOne(mon) then
        return zeroR;
      fi;

      ##########################
      ImageOfPart := function(mon)
        local exp, indetnum, unimon;
        unimon := IndeterminateAndExponentOfUnivariateMonomial(mon);
        indetnum := IndeterminateNumberOfUnivariateRationalFunction(unimon[1]);
        exp := unimon[2];
        if exp = 1 then
          return images[indetnum];
        else
          return exp*unimon[1]^(exp-1)*images[indetnum];
        fi;
      end; 
      ##########################

      image := zeroR;
      for i in UnivariateMonomialsOfMonomial(mon) do
        im := oneR;
        for j in UnivariateMonomialsOfMonomial(mon) do
          if j = i then
            im := im * ImageOfPart(i);
          else
            im := im * j;
          fi;
        od;
        image := image + im;
      od;
      return image; 
    end;
    ##################################

    image := zeroR;
    for t in TermsOfPolynomial(poly) do
      image := image + t[2] * oneR * ImageOfMonomial(t[1]);
    od;
  
    return image;
  end
);
#####################################################################



#if LoadPackage("singular") = true then
if IsPackageMarkedForLoading("singular","0") then
  #####################################################################
  ##  <#GAPDoc Label="KernelOfDerivation_DTmanDerivation_Hom">
  ##  <ManSection>
  ##  <Oper Name="KernelOfDerivation" Arg="d [, avoid]"/>
  ##
  ##  <Returns>
  ##  List 
  ##  </Returns>
  ##  <Description>
  ##  Returns a ring presentation <M>S/J</M> for the kernel of the derivation 
  ##  <A>d</A>, where <M>S</M> is a polynomial ring and <M>J</M> are a set of
  ##  relations.
  ##  This operation returns a list with the following elements:
  ##  <Enum>
  ##    <Item>the new polynomial ring <M>S</M></Item>
  ##    <Item>a basis for the ideal, as a set of relations <M>J</M></Item>
  ##    <Item>the ring isomorphism between the kernel (i.e. a subring of the 
  ##       derivation's ring) and the new ring. This is given using the
  ##       <K>HAPRingHomomorphism</K> type, for details of which see
  ##       <Ref Chap="RingHomomorphism"/></Item>
  ##  </Enum>
  ##  An optional parameter, <A>avoid</A>, can be provided which lists
  ##  indeterminates to avoid when creating the the new polynomial ring.
  ##  <P/>
  ##  <E>This function is only available if the package
  ##  <Package>singular</Package> is available</E>.
  ##  </Description>
  ##  </ManSection>
  ##  <#/GAPDoc>
  #####################################################################
  InstallOtherMethod(KernelOfDerivation, 
    "for no alternative indeterminates",
    [IsHAPDerivation],
    function(d)
      return KernelOfDerivation(d, []);
    end
  );
  #####################################################################
  InstallMethod(KernelOfDerivation, 
    "with specified alternative indeterminates",
    [IsHAPDerivation, IsHomogeneousList],
    function(d, avoid)

      local ring, ideal, p, gens, g, img, B, coeffs, A, i, K, k, oneR, zeroR, 
      indets, kernelgenerators, kernelindets, newkernelring, kernelrelations, 
      foundunity, kernelindetexps, message, Smod, l, j, h, kernelmap, 
      reductionmap, combinedmap;

      ring := DerivationRing(d);
      indets := IndeterminatesOfPolynomialRing(ring);
      ideal := DerivationRelations(d);
      p := Characteristic(CoefficientsRing(ring));
      if IsZero(p) then
        Error("can only calculate the kernel of rings with a coefficient ring of positive characteristic");
      fi;
      oneR := One(ring);
      zeroR := Zero(ring);

      # Remember the smallext power of each exponent that we know is in the 
      # kernel
      kernelindetexps := [];  
      for i in indets do
        if IsZero(ImageOfDerivation(d, i)) then
          Add(kernelindetexps, 1);
        else
          Add(kernelindetexps, p);
        fi;
      od;


      # If all the images are zero, then the kernel is the whole ring
      if ForAll(kernelindetexps, IsOne) then
        newkernelring := PolynomialRing(
          CoefficientsRing(ring), Length(indets), 
          Concatenation(indets, avoid));
        kernelmap := HAPRingHomomorphismByIndeterminateMap(
          ring, ideal, newkernelring);

        if InfoLevel(InfoHAPprime) >= 2 then
          Info(InfoHAPprime, 2, "Zero derivation: Kernel is whole ring");
          message := "Change of indeterminates: ";
          for i in [1..Length(SourceGenerators(kernelmap))] do 
            message := Concatenation(
              message, String(SourceGenerators(kernelmap)[i]), "=", 
                String(ImageGenerators(kernelmap)[i]), ", ");
          od;
          Remove(message, Length(message)); # Remove the last space
          Remove(message, Length(message)); # Remove the last comma
          Info(InfoHAPprime, 2, message);
        fi;
        return [ImagePolynomialRing(kernelmap), 
          ImageRelations(kernelmap), kernelmap];
      fi;

      # We can now start to make a list of our kernel generators by using 
      # kernelindetexps
      kernelgenerators := [];
      for i in [1..Length(indets)] do
        Add(kernelgenerators, indets[i]^kernelindetexps[i]);
      od;

      # We shall use Singular to find the kernel
      # of the derivation, by finding the kernel of a S-module
      # homomorphism

      Smod := HAPPRIME_SModule(ring, kernelindetexps);

      # Calculate the derivation of each of our S-module generators, 
      # and add it to A in S-module form
      A := [];
      for g in Smod.generators do
        img := ImageOfDerivation(d, g);
        Add(A, Smod.FromPoly(img));
      od;
      # And needs to be column-major
      A := TransposedMat(A);

      # B is made up of the generators of the ideal (again, in 
      # S-module form)
      # We need the multiple with each of the S-module generators
      # as well as the raw relations
      if IsEmpty(ideal) then
        B := [ListWithIdenticalEntries(Length(Smod.generators), zeroR)];
      else
        B := [];
        for i in ideal do
          for g in Smod.generators do
            Add(B, Smod.FromPoly(i*g));
          od;
        od;
      fi;
      # And needs to be column-major
      B := TransposedMat(B);

      # Now calculate the kernel of the module homomorphism using singular
      SingularSetBaseRing(Smod.Sring);
      K := SingularInterface("modulo", [A, B], "matrix");

      # And turn it back into GAP row-major
      K := TransposedMat(K);

      # And now the relations from our kernel
      foundunity := false;
      for coeffs in K do
        g := Smod.ToPoly(coeffs);
        if g = oneR then
          foundunity := true;
        else
          Add(kernelgenerators, g);
        fi;
      od;
      if not foundunity then
        Error("no unity element in the kernel");
      fi; 
#      Info(InfoHAPprime, 2, "Raw kernel generators are ", kernelgenerators);

      # Our kernel ideal is the old ring ideal in the derivation, plus the 
      # relations that explain how our new indeterminates relate to the old ones
      # Make sure there's no zeros in here, since that can cause problems later
      kernelrelations := Filtered(ideal, i->not IsZero(i));
      # And make sure that it's a Groebner Basis
      SetTermOrdering(ring, "dp");
      kernelrelations := SingularReducedGroebnerBasis(
        Ideal(ring, kernelrelations));


      # Now reduce this list
      # First get rid of anything in the ideal
      SingularSetNormalFormIdealNC(Ideal(ring, kernelrelations));
      i := 1;
      l := Length(kernelgenerators);
      while i <= l do
        if IsZero(SingularPolynomialNormalForm(kernelgenerators[i])) then
          # If it's in the ideal, remove it
          Remove(kernelgenerators, i);
          l := l - 1;
        else
          i := i + 1;
        fi;
      od;

      # We now want to remove any generators that are products of other ones
      # that we have.
      # Sort the generators into increasing degree so that the simplest
      # ones are tested first
      Sort(kernelgenerators);
      i := 2;
      l := Length(kernelgenerators);
      while i <= l do
        # Get the next generator
        g := kernelgenerators[i];
        # See if it can be divided by anything else in our list
        # We only need to consider factors before this in the list since
        # any factor has to have a smaller degree
        j := 1;
        while j < i do
          h := g / kernelgenerators[j];
          if IsOne(DenominatorOfRationalFunction(h)) then
            # the division was OK, so keep the answer and try again
            g := h;
          else
            # try another divisor
            j := j + 1;
          fi;
        od;
        if IsOne(g) then
          # It is entirely generated as a product of other generators,
          # so we remove it
          Remove(kernelgenerators, i);
          l := l - 1;
        else
          # Otherwise we keep it
          i := i + 1;
        fi;
      od;

      # Now assign indeterminates to each of the kernel generators
      # There will be one indeterminate for each entry in kernelgenerators
      newkernelring := PolynomialRing(
        CoefficientsRing(ring), Length(kernelgenerators), 
        Concatenation(indets, avoid));

      # And create the isomorphism between our kernel in the original ring
      # and our kernel in this ring
      # now make the map
      kernelmap := HAPSubringToRingHomomorphism(
        kernelgenerators, kernelrelations, newkernelring);

      # But the image of this map may not be in reduced form, so
      # compute the reduced form of this
      reductionmap := HAPRingReductionHomomorphism(kernelmap);

      combinedmap := CompositionRingHomomorphism(kernelmap, reductionmap);  

      if InfoLevel(InfoHAPprime) >= 2 then  
        for i in [1..Length(SourceGenerators(combinedmap))] do
          Info(InfoHAPprime, 2, 
            "Kernel generator ", SourceGenerators(combinedmap)[i], " = ", 
              ImageGenerators(combinedmap)[i]);
        od;
      fi;

      return [ImagePolynomialRing(combinedmap), 
        ImageRelations(combinedmap), combinedmap];
    end
  );
  #####################################################################



  #####################################################################
  ##  <#GAPDoc Label="HomologyOfDerivation_DTmanDerivation_Hom">
  ##  <ManSection>
  ##  <Oper Name="HomologyOfDerivation" Arg="d[, avoid]"/>
  ##
  ##  <Returns>
  ##  List 
  ##  </Returns>
  ##  <Description>
  ##  Returns a ring presentation <M>S/J</M> for the homology of the derivation 
  ##  <A>d</A>, where <M>S</M> is a polynomial ring and <M>J</M> are a set of
  ##  relations.
  ##  Returns a polynomial ring presentation for the homology <M>ker(d)/im(d)</M> 
  ##  of the derivation <A>d</A>.
  ##  This operation returns a list with the following elements:
  ##  <Enum>
  ##    <Item>the new polynomial ring <M>S</M></Item>
  ##    <Item>a basis for the ideal, as a set of relations <M>J</M></Item>
  ##    <Item>the ring isomorphism between the kernel (i.e. a subring of the 
  ##       derivation's ring) and the new ring. This is given using the
  ##       <K>HAPRingHomomorphism</K> type, for details of which see
  ##       <Ref Chap="RingHomomorphism"/></Item>
  ##  </Enum>
  ##  <P/>
  ##  An optional parameter, <A>avoid</A>, can be provided which lists
  ##  indeterminates to avoid when creating the the new polynomial ring.
  ##  </Description>
  ##  </ManSection>
  ##  <#/GAPDoc>
  #####################################################################
  InstallOtherMethod(HomologyOfDerivation, 
    "for no alternative indeterminates",
    [IsHAPDerivation],
    function(d)
      return HomologyOfDerivation(d, []);
    end
  );
  #####################################################################
  InstallMethod(HomologyOfDerivation, 
    "with a set of alternative indeterminates",
    [IsHAPDerivation, IsHomogeneousList],
    function(d, avoid)
  
      local ring, gens, imgens, g, i, kernel, kernelindetexps, reductionmap, p, 
      Smod, message, totalmap, totalgens, totalrels;
  
      ring := DerivationRing(d);
      p := Characteristic(CoefficientsRing(ring));

      # If all the images are zero, then the homology is the whole ring
      # (and is the same as the kernel)
      if IsZero(Compacted(DerivationImages(d))) then
        # Just return the kernel
        kernel := KernelOfDerivation(d, avoid);
        Info(InfoHAPprime, 2, "Zero derivation: Homology is whole ring");
        return kernel;
      fi;

      # Remember the smallest power of each exponent that we know is in the 
      # kernel
      kernelindetexps := [];  
      for i in IndeterminatesOfPolynomialRing(ring) do
        if IsZero(ImageOfDerivation(d, i)) then
          Add(kernelindetexps, 1);
        else
          Add(kernelindetexps, p);
        fi;
      od;
      # Find generators for the image. 
      # These are the images of all of the generators of the ring as an R-module
      imgens := [];
      Smod := HAPPRIME_SModule(ring, kernelindetexps);
      for g in Smod.generators do
        AddSet(imgens, ImageOfDerivation(d, g));
      od;
      # if there's a zero, remove it since it's not interesting
      i := Position(imgens, Zero(ring));
      if i <> fail then
        Remove(imgens, i);
      fi;
  
      # If the image includes one then the answer is just the trivial
      # ring
      for i in imgens do
        if IsOne(i) then
          return [Ring(Zero(i)), [ ], 
          HAPZeroRingHomomorphism(ring, DerivationRelations(d)) ];
        fi;
      od;
  
      # Calculate the kernel of the derivation
      kernel := KernelOfDerivation(d, avoid);
      Info(InfoHAPprime, 2, "Kernel relations ", kernel[2]);
  
      # The image is expressed in terms of the indeterminates of d.
      # We want to convert this into the indeterminates of the new kernel 
      # presentation. (This will always be
      # possible since the image is a subring of the kernel - provided that
      # the derivation squares to zero).
      imgens := ImageOfRingHomomorphism(kernel[3], imgens);
      Info(InfoHAPprime, 2, "Image relations ", imgens);

      # Add these to the ideal from the kernel and calculate the Grobner basis
      # to tidy things up
      imgens := SingularGroebnerBasis(
        Ideal(kernel[1], Concatenation(imgens, kernel[2])));

      # Now see if we can reduce the presentation
      reductionmap := HAPRingReductionHomomorphism(kernel[1], imgens, 
        Concatenation(IndeterminatesOfPolynomialRing(ring), avoid));
      if InfoLevel(InfoHAPprime) >= 2 then
        message := "Change of indeterminates: ";
        for i in [1..Length(SourceGenerators(reductionmap))] do 
          message := Concatenation(
            message, String(SourceGenerators(reductionmap)[i]), "=", 
              String(ImageGenerators(reductionmap)[i]), ", ");
        od;
        Remove(message, Length(message)); # Remove the last space
        Remove(message, Length(message)); # Remove the last comma
        Info(InfoHAPprime, 2, message);
      fi;

      # We want to return a mapping between the homology in the original ring 
      # indeterminates and the final reduced version.
      # Feed the image of the homology map back through and through the kernel
      # map to see what the generators of the kernel are.
      totalgens := 
        List(IndeterminatesOfPolynomialRing(ImagePolynomialRing(reductionmap)), 
          i->PreimageOfRingHomomorphism(kernel[3], 
            PreimageOfRingHomomorphism(reductionmap, i)));
      # For the relations, we keep the kernel relations, but we also want
      # any from the image that end up in the homology. So, we also feed the 
      # imgens the homology back through the kernel and define them at the
      # source (removing any zeros). This may include 
      # indeterminates that do not feature in the generators of the homology,
      # but which are killed in the homology, but that's OK
      totalrels := Concatenation(SourceRelations(kernel[3]),
        Filtered(List(SourceRelations(reductionmap), 
          i->PreimageOfRingHomomorphism(kernel[3], i)), j->not IsZero(j)));
      # now build the map
      totalmap := HAPSubringToRingHomomorphism(
        totalgens, totalrels, ImagePolynomialRing(reductionmap));

      return [
        ImagePolynomialRing(reductionmap), 
        ImageRelations(reductionmap),
        totalmap];
    end
  );
  #####################################################################
else
  InstallMethod(KernelOfDerivation, 
    [IsHAPDerivation, IsHomogeneousList],
    function(d, indetsAlt)
      Error("The package 'singular' cannot be loaded, so this 'HAPprime' function is not available");
    end
  );
  #####################################################################
  InstallMethod(HomologyOfDerivation, 
    [IsHAPDerivation, IsHomogeneousList],
    function(d, indetsAlt)
      Error("The package 'singular' cannot be loaded, so this 'HAPprime' function is not available");
    end
  );
fi;
#####################################################################




#####################################################################
##  <#GAPDoc Label="HAPPRIME_SModule_DTmanDerivationInt">
##  <ManSection>
##  <Func Name="HAPPRIME_SModule" Arg="R, exps"/>
##
##  <Returns>
##  Record
##  </Returns>
##  <Description>
##  For a polynomial ring <A>R</A>, <M>k[x_1, x_2, ..., x_n]</M>, 
##  and list of exponents <A>exps</A>, <M>[e_1, e_2, ..., e_n]</M>, returns
##  a record which represents <A>R</A> as an <M>S</M>-module, where
##  <M>S</M> is the subring of 
##  <A>R</A> given by <M>k[x_1^{e_1}, x_2^{e_2}, ..., x_n^{e_n}]</M>
##  <P/>
##  The record has the following components:
##  <List>
##    <Item><C>Rring</C> the original ring <A>R</A></Item>
##    <Item><C>Sring</C> a ring isomorphic to the subring of powers <M>S</M></Item>
##    <Item><C>Spows</C> the exponents <A>exps</A></Item>
##    <Item><C>generators</C> the list of generators of <A>R</A> as an 
##      <M>S</M>-module</Item>
##    <Item><C>FromPoly</C> a function which takes a polynomial in <A>R</A> and
##      returns the corresponding element in the <M>S</M>-module (as a 
##      vector)</Item>
##    <Item><C>ToPoly</C> a function which takes an element in the 
##      <M>S</M>-module (as a vector) and returns the corresponding polynomial
##      in <A>R</A></Item>
##  </List>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_SModule,
  function(R, exps)
  
    local indets, n, field, p, Sring, Sindets, gens, oneR, zeroR, exps2, j, i,
    gen, PolynomialToSModule, SModuleToPolynomial;
    
    if not IsPolynomialRing(R) then
      Error("<R> must be a polynomial ring");
    fi;
    indets := IndeterminatesOfPolynomialRing(R);
    n := Length(indets);
    if Length(exps) <> n or not IsHomogeneousList(exps) then
      Error("<exps> must be a list of (positive integer) exponents for the indeterminates of <R>");
    fi;

    field := CoefficientsRing(R);
    p := Characteristic(field);
    if IsZero(p) then
      Error("the S-module representation is only available for positive characteristics");
    fi;

    # We shall use a different set of indeterminates for the
    # subring. Create this.
    Sring := PolynomialRing(field, n, indets);

    Sindets := IndeterminatesOfPolynomialRing(Sring);

    # Create the set of generators
    gens := [];
    oneR := One(R);
    zeroR := Zero(R);
    # exps is the list of exponents for this generator
    exps2 := ListWithIdenticalEntries(n, 0);
    repeat
      gen := oneR;
      for j in [1..n] do
        gen := gen * indets[j]^exps2[j]; 
      od;
      Add(gens, gen);
      # Now move to the next generator
      i := 1;
      repeat
        exps2[i] := exps2[i] + 1;
        if exps2[i] = exps[i] then
          exps2[i] := 0;
          i := i + 1;
        else
          i := n+2; # this means we have successfully updated
        fi;
      until i > n;
    until i = n+1;
    
    ###################
    # Converts an element of R into an element in the S-module
    # (as a vector)
    PolynomialToSModule := function(poly)
      local coeffs, terms, i, modrep, 
        MonomialToSModuleRep;

      if IsZero(poly) then
        # note that zeroR = zeroS and oneR = oneS
        return ListWithIdenticalEntries(Length(gens), zeroR);
      fi;
      if IsOne(poly) then
        coeffs := ListWithIdenticalEntries(Length(gens), zeroR);
        coeffs[1] := oneR;
        return coeffs;
      fi;

      ##################################
      # Convert a monomial to module coefficients
      # Returns a list with as the first element the basis number, and as
      # the second element the coefficient of that basis element
      MonomialToSModuleRep := function(mon)
        local umon, gen, coeff, i, inum, exp, Spow;

        umon := UnivariateMonomialsOfMonomial(mon);

        if IsZero(umon[1]) then
          return [1, Zero(field)];
        fi;
        if IsOne(umon[1]) then
          return [1, One(field)];
        fi;

        gen := oneR;
        coeff := oneR;
        for i in umon do
          i := IndeterminateAndExponentOfUnivariateMonomial(i);
          inum := Position(indets, i[1]);
          exp := i[2];

          # How much of this indeterminate in the Sring do we have?
          Spow := QuoInt(exp, exps[inum]);

          if Spow > 0 then
            # Copy across the number of powers of p of this indeterminant in 
            # this exponent
            coeff := coeff * Sindets[inum]^Spow;
          fi;
          # And build up the generator that this corresponds to
          gen := gen * indets[inum]^(exp - Spow*exps[inum]);
        od;

        return [Position(gens, gen), coeff];
      end;
      ##################################

      coeffs := ListWithIdenticalEntries(Length(gens), zeroR);
      terms := TermsOfPolynomial(poly);
      for i in terms do
        modrep := MonomialToSModuleRep(i[1]);
        coeffs[modrep[1]] := coeffs[modrep[1]] + i[2]*modrep[2];
      od;
      
      return coeffs;
    end;
    ###################
    
    ###################
    # Converts a vector from the S-module into a polynomial element of R
    SModuleToPolynomial := function(coeffs)
      local poly, i, coeff, SMonomialToRMonomial, SPolynomialToRPolynomial;

      ##################################
      # Convert a monomial in S to a monomial in R
      SMonomialToRMonomial := function(mon)
        local Rmon, umon, i, inum;
        Rmon := oneR;
        umon := UnivariateMonomialsOfMonomial(mon);
        for i in umon do
          i := IndeterminateAndExponentOfUnivariateMonomial(i);
          inum := Position(Sindets, i[1]);
          Rmon := Rmon * indets[inum]^(i[2] * exps[inum]);
        od;
        return Rmon;
      end; 
      ##################################
      # Convert a polynomial in S to a polynomial in R
      SPolynomialToRPolynomial := function(poly)
        local Rpoly, terms, i, mon;
        if IsZero(poly) or IsOne(poly) then
          return poly;
        fi;
        Rpoly := zeroR;
        terms := TermsOfPolynomial(poly);
        for i in terms do
          mon := SMonomialToRMonomial(i[1]);
          Rpoly := Rpoly + i[2]*mon;
        od;
        return Rpoly;
      end; 
      ##################################
      poly := zeroR; 
      for i in [1..Length(gens)] do
        if not IsZero(coeffs[i]) then
          coeff := SPolynomialToRPolynomial(coeffs[i]);
          poly := poly + gens[i]*coeff;
        fi;
      od;
      return poly;
    end;
    ###################

    return rec(
    Rring := R,
    Sring := Sring,
    Spows := exps,
    generators := gens,
    FromPoly := PolynomialToSModule,
    ToPoly := SModuleToPolynomial);
  end
);
#####################################################################


