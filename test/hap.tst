gap> START_TEST("HAP library");
gap> GroupHomology(AlternatingGroup(5),2,2);
[ 2 ]
gap> GroupHomology(SmallGroup(32,3),4,2);
[ 2, 2, 2, 2, 2 ]
gap> GroupHomology(AbelianGroup([2,4,6]),4);
[ 2, 2, 2, 2, 2, 2, 2, 2 ]
gap> GroupHomology([[1,[2,3]],[2,[3,4]]],2);
[ 0, 0 ]
gap> R:=ResolutionNilpotentGroup(SmallGroup(64,135),5);;
gap> TR:=TensorWithIntegers(R);;
gap> H4:=Homology(TR,4);
[ 2, 2, 2, 2, 2, 2, 2 ]
gap> H4a:=IntegralRingGenerators(R,4);
[ [ 0, 0, 1, 1, 0, -4, -4 ], [ 0, 0, 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 0, 0, 1 ]
 ]
gap> Tensor:=ThirdHomotopyGroupOfSuspensionB(SymmetricGroup(3));
[ 2 ]
gap> Tensor2:=ThirdHomotopyGroupOfSuspensionB(SymmetricGroup(3),12);
[ 2 ]
gap> P:=PresentationOfResolution(ResolutionAbelianGroup([0,2,4],3)).relators;
[ f1*f3*f1*f3^-1, f2*f3*f2^-1*f3^-1, f1^2, f2*f1*f2^-1*f1, f2^4 ]
gap> P:=Length(P);
5
gap> R:=ResolutionAbelianGroup([0,0,0],5);;
gap> H6:=Homology(TensorWithIntegers(R),3);
[ 0 ]
gap> H7:=Cohomology(HomToIntegersModP(R,2),3);
1
gap> IsSuperperfect(PerfectGroup(120,1));
true
gap> CoefficientsOfUnivariateRationalFunction(PoincareSeries(SmallGroup(64,134)));
[ [ 1 ], [ 1, -3, 3, -1 ], 0 ]
gap> STOP_TEST( "tst.tst", 1000 );


