#############################################################################
##
##  HAPPRIME - poincare.gi
##  Operations and Methods for Poincare Series
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
##  <#GAPDoc Label="PoincareSeriesLHS_manPoincare">
##  <ManSection>
##  <Attr Name="PoincareSeriesLHS" Arg="G"/>
##
##  <Returns>
##    Rational function
##  </Returns>
##  <Description>
##  For a finite p-group <A>G</A>, this function calculates 
##  and returns a quotient of polynomials <M>f(x) = P(x)/Q(x)</M> 
##  (i.e. the Poincaré series) whose 
##  coefficient of <M>x^k</M> equals the rank of the vector space 
##  <M>H_k(G, \mathbb{F}_p)</M> for all <M>k</M> in the range <M>k=1</M> to 
##  <M>k=n</M>. 
##  <P/>
##  This function computes a Lyndon-Hoschild-Serre spectral sequence for 
##  the &p;-group &G;. The last sheet of
##  this sequence will have the same additive structure as the mod-<M>p</M> 
##  group cohomology ring of &G;, and thus the same Poincaré
##  series, which is returned by this function. 
##  <P/>
##  See Section <Ref Sect="ExPoincare"/> for an example and more description.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(PoincareSeriesLHS,
  "generic method",
  [IsGroup],
  function(G)
####
#### Added by Graham Ellis (22/12/2008)
    if SSortedList(Factors(Order(G)))<>[2] then
    Print("This function is currently implemented only for 2-groups.\n");
    return fail;
    fi;
####
####
    return HilbertPoincareSeries(LHSSpectralSequenceLastSheet(G));
end);
#####################################################################

if false then
#####################################################################
##  <#GAPDoc Label="PoincareSeriesAutoMem_DTmanPoincare">
##  <ManSection>
##  <Oper Name="PoincareSeriesAutoMem" Arg="G [, n]"/>
##
##  <Returns>
##    Rational function
##  </Returns>
##  <Description>
##  For a finite p-group <A>G</A>, this function calculates 
##  and returns a quotient of polynomials <M>f(x) = P(x)/Q(x)</M> 
##  (i.e. the Poincaré series) whose 
##  coefficient of <M>x^k</M> equals the rank of the vector space 
##  <M>H_k(G, \mathbb{F}_p)</M> for all <M>k</M> in the range <M>k=1</M> to 
##  <M>k=n</M>. If no value is given for <A>n</A> then increasing values of 
##  <M>n</M> are tried to find the minimum value which gives a consistent 
##  Poincaré series, defined as a the minimum value <M>n &gt; 10</M> such that
##  <C>PoincareSeries(G, n) = PoincareSeries(G, n-1) = PoincareSeries(G, n-2)</C>.
##  <P/>
##  This function uses the &HAP; function 
##  <Ref Func="PoincareSeries" BookName="HAP"/> to calculate the Poincaré
##  series, but the &HAPprime; function
##  <Ref Oper="ExtendResolutionPrimePowerGroupAutoMem" BookName="HAPprime"/> to calculate
##  and gradually extend the resolution, so should be both faster and more 
##  memory-efficient than using <K>PoincareSeries</K> by itself.
##  See Section <Ref Sect="PoincareExample"/> for an example.
##  <P/>
##  The Poincaré series calculated using this function is likely to be correct, 
##  but we have no proof that this will be the case. If a correct Poincaré
##  series is required, use <Ref Oper="PoincareSeriesLHS" BookName="HAPprime"/>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallMethod(PoincareSeriesAutoMem,
  "generic method",
  [IsGroup, IsPosInt],
  function(G, n)
    return PoincareSeries(ResolutionPrimePowerGroupAutoMem(G, n));
end);
#####################################################################
InstallMethod(PoincareSeriesAutoMem,
  "generic method",
  [IsGroup],
  function(G)
  
  local factors, numToCompare, R, 
    dims, length, Pseries, found, match, i, t;
  
  Info(InfoHAPprime, 2, "Using HAPprime version ", HAPPRIME_VersionWithSVN()); 
  Info(InfoHAPprime, 2, "Max available memory ", GasmanLimits().max); 
  t := Runtime();
  
  # If the group factors, then I can just do the Poincare series of the factors
  # and then find their product.
  factors := DirectFactorsOfGroup(G);
  if Length(factors) > 1 then  
    Info(InfoHAPprime, 1, "Splitting into ", Length(factors), " factors to find Poincaré series"); 
    factors := List(factors, g->PoincareSeriesAutoMem(g));
    return Product(factors);
  else

    # Build a resolution estimate its PoincareSeries
    length := 10;
    numToCompare := 3;
    R := ResolutionPrimePowerGroupAutoMem(G, length);
    dims := List([0..length], i->ResolutionModuleRank(R, i));

    Pseries := [];
    for i in [1..numToCompare] do
      Pseries[i] := PoincareSeries(dims, length-numToCompare+i);
    od;
    Info(InfoHAPprime, 2, "Time taken so far: ", Runtime() - t, "ms"); 

    found := false;
    while not found do
      match := List([1..(numToCompare-1)], i->Pseries[i] = Pseries[i+1]);

      if not false in match then 
        found := true;
      else
        # Find the next stage in the resolution
        R := ExtendResolutionPrimePowerGroupAutoMem(R);
        length := length + 1;
        Add(dims, ResolutionModuleRank(R, length));
        Remove(Pseries, 1);
        Pseries[numToCompare] := PoincareSeries(dims, length);
       Info(InfoHAPprime, 2, "Time taken so far: ", Runtime() - t, "ms"); 
     fi;
    od;
  fi;
  return Pseries[1];
end);
#####################################################################
InstallMethod(PoincareSeriesAutoMemStop,
  "generic method",
  [IsGroup],
  function(G)
  
  local orderG, nGgens, factors, numToCompare, R, 
    dims, length, Pseries, found, match, i, stop,
    size, radicalcomputebyes, memfree, t;

  Info(InfoHAPprime, 2, "Using HAPprime version ", HAPPRIME_VersionWithSVN()); 
  Info(InfoHAPprime, 2, "Max available memory ", GasmanLimits().max); 
  t := Runtime();
  
  orderG := Order(G);
  nGgens := Length(MinimalGeneratingSet(G));
  # If the group factors, then I can just do the Poincare series of the factors
  # and then find their product.
  factors := DirectFactorsOfGroup(G);
  if Length(factors) > 1 then  
    Info(InfoHAPprime, 1, "Splitting into ", Length(factors), " factors to find Poincaré series"); 
    factors := List(factors, g->PoincareSeriesAutoMem(g));
    return Product(factors);
  else

    # Build a resolution estimate its PoincareSeries
    length := 10;
    numToCompare := 3;
    R := ResolutionPrimePowerGroupAutoMem(G, length);
    dims := List([0..length], i->ResolutionModuleRank(R, i));

    Pseries := [];
    for i in [1..numToCompare] do
      Pseries[i] := PoincareSeries(dims, length-numToCompare+i);
    od;
    Info(InfoHAPprime, 2, "Time taken so far: ", Runtime() - t, "ms"); 

    found := false;
    stop := false;
    while not found and not stop do
      match := List([1..(numToCompare-1)], i->Pseries[i] = Pseries[i+1]);

      if not false in match then 
        found := true;
      else
        # Find the next stage in the resolution
        # How big is the current one?
        size := ResolutionModuleRank(R, ResolutionLength(R));
        radicalcomputebyes := Int(nGgens * orderG * orderG * size * size / 5); # In bytes
        GASMAN("collect");
        memfree := GasmanStatistics().full.freekb;
        memfree := memfree + GasmanLimits().max - GasmanStatistics().full.totalkb;
        memfree := memfree * 1024; # convert to Mb
        # Allow a factor of two spare
        if memfree > (radicalcomputebyes * 2) then
          R := ExtendResolutionPrimePowerGroupRadical(R);
        else
          stop := true;
        fi;

        if not stop then
          length := length + 1;
          Add(dims, ResolutionModuleRank(R, length));
          Remove(Pseries, 1);
          Pseries[numToCompare] := PoincareSeries(dims, length);
        fi;
        Info(InfoHAPprime, 2, "Time taken so far: ", Runtime() - t, "ms"); 
      fi;
    od;
    if stop then
      Info(InfoHAPprime, 1, "Stopped: about to run out of memory for Radical"); 
      return fail;
    fi;
  fi;
  return Pseries[1];
end);
#####################################################################
fi;
  
