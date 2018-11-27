gap> G:=SylowSubgroup(MathieuGroup(12),2);;
gap> FG:=GroupRing(GF(2),G);;
gap> g:=LieAlgebra(FG);;
#I  LAGUNA package: Constructing Lie algebra ...
gap> L:=LeibnizComplex(g,3);;
gap> BettiNumber(L,1);
16
