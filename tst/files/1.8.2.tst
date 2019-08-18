gap> START_TEST("HAP library");
gap> A:=[[0,-1,0,0,],[1,0,0,0,],[0,0,0,1],[0,0,-1,0]];;
gap> B:=[[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0]];;
gap> Q:=Group([A,B]);;
gap> Y:=EquivariantOrbitPolytope(Q,[1,0,0,0]);
Equivariant CW-complex of dimension 4
gap> for k in [0..2] do
> for n in [1..Y!.dimension(k)] do
> Print(Order(Y!.stabilizer(k,n)),"  ");
> od;od;
1  1  1  1  1  1  1  1 
gap>  F:=FundamentalGroupOfQuotient(Y);
<fp group on the generators [ f1, f2, f3 ]>
gap> RelatorsOfFpGroup(F);
[ f1*f2^-1*f3^-1, f1*f3^-1*f2, f1*f3*f2^-1, f1*f2*f3 ]
gap> STOP_TEST( "tst.tst", 1000 );


