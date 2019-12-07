gap> START_TEST("HAP library");
gap> G:=SylowSubgroup(MathieuGroup(12),2);;
gap> R:=ResolutionPrimePowerGroup(G,10);;
gap> A:=ModPCohomologyRing(R);
<algebra of dimension 286 over GF(2)>
gap> S:=[];;
gap> for n in [0..10] do
> B:=Filtered(Basis(A),x->A!.degree(x)=n);
> B:=Subalgebra(A,B);
> Add(S,B);
> od;
gap> List(S,Dimension);
[ 1, 75, 100, 66, 60, 83, 28, 36, 45, 55, 66 ]
gap> STOP_TEST( "tst.tst", 1000 );
