gap> START_TEST("HAP library");
gap> CohomologicalData(SmallGroup(8,4),20);
Group order: 8
Group number: 4
Group description: Q8

Cohomology generators
Degree 1: a, b
Degree 4: c

Cohomology relations
1: a^2+a*b+b^2
2: b^3

Poincare series
(x^2+x+1)/(-x^3+x^2-x+1)

Steenrod squares
Sq^1(c)=0
Sq^2(c)=0
gap> a:=[[0,-1,0,0,],[1,0,0,0,],[0,0,0,1],[0,0,-1,0]];;
gap> b:=[[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0]];;
gap> G:=Group([a,b]);;
gap> rho:=IdentityMapping(G);;
gap> v:=[1,2,3,4];;
gap> A:=Mod2SteenrodAlgebra(G,10);;
gap> sw:= FundamentalMultiplesOfStiefelWhitneyClasses(rho,v,A,true);
[ [ v.1 ], [ v.2+v.3 ], [ v.5 ], [ 0*v.1 ], [ v.7 ] ]
gap> TotalStiefelWhitneyClass:= sw[1][1]+sw[2][1]+sw[3][1]+sw[4][1]+sw[5][1];
v.1+v.2+v.3+v.5+v.7
gap> STOP_TEST( "tst.tst", 1000 );
