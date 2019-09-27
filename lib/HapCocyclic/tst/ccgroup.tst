#############################################################################
##
#W  ccgroup.tst                   HAPCocyclic                 Robert F. Morse
##
##

gap> ## Test creating Cc-Groups with Standard 2-Cocycle
gap> ##
gap> H := SmallGroup(256,2);;
gap> N := NormalSubgroups(H)[16];;
gap> sco := Standard2Cocycle(H,N);;
gap> ccg := CcGroup(Module(sco),sco);;
gap> IdGroup(ccg);
[ 256, 2 ]
gap> ## Create cocycles with HAP
gap> ##
gap> OG := GOuterGroup(H,N);;
gap> A := Centre(OG);;
gap> G:=ActingGroup(A);;
gap> R:=ResolutionFiniteGroup(G,3);;
gap> C:=HomToGModule(R,A);;
gap> CH:=CohomologyModule(C,2);;
gap> Elts:=Elements(ActedGroup(CH));;
gap> Length(Elts);
256
gap> lst := List(Elts{[1..5]},x->CH!.representativeCocycle(x));;
gap> ccgrps := List(lst, x->CcGroup(OG, x));;
gap> List(ccgrps,IdGroup);
[ [ 256, 3700 ], [ 256, 3 ], [ 256, 3783 ], [ 256, 1300 ], [ 256, 1553 ] ]
gap> STOP_TEST( "ccgroup.tst", 38100000 );
