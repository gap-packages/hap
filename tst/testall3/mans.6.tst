gap> START_TEST("HAP library");
gap> L13_1:=ThreeManifoldViaDehnSurgery([[1,2],[1,2]],13,1);;
gap> L13_2:=ThreeManifoldViaDehnSurgery([[1,2],[1,2]],13,2);;
gap> L13_1:=BarycentricSubdivision(L13_1);;
gap> L13_2:=BarycentricSubdivision(L13_2);;
gap> A13_1:=CohomologyRing(L13_1,13);;
gap> A13_2:=CohomologyRing(L13_2,13);;
gap> M13_1:=List([1..4],i->[]);;
gap> B13_1:=CanonicalBasis(A13_1);;
gap> M13_2:=List([1..4],i->[]);;
gap> B13_2:=CanonicalBasis(A13_2);;
gap> for i in [1..4] do
> for j in [1..4] do
> M13_1[i][j]:=B13_1[i]*B13_1[j];
> od;od;
gap> for i in [1..4] do
> for j in [1..4] do
> M13_2[i][j]:=B13_2[i]*B13_2[j];
> od;od;
gap> Display(M13_1);
[ [    v.1,    v.2,    v.3,    v.4 ],
  [    v.2,  0*v.1,    v.4,  0*v.1 ],
  [    v.3,    v.4,  0*v.1,  0*v.1 ],
  [    v.4,  0*v.1,  0*v.1,  0*v.1 ] ]
gap> Display(M13_2);
[ [            v.1,            v.2,            v.3,            v.4 ],
  [            v.2,          0*v.1,  (Z(13)^7)*v.4,          0*v.1 ],
  [            v.3,  (Z(13)^7)*v.4,          0*v.1,          0*v.1 ],
  [            v.4,          0*v.1,          0*v.1,          0*v.1 ] ]
gap> STOP_TEST( "tst.tst", 1000 );

