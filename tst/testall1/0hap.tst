##A quick test of the main older group (co)homology functions
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
gap> R:=ResolutionAbelianGroup([0,0,0],5);;
gap> H6:=Homology(TensorWithIntegers(R),3);
[ 0 ]
gap> H7:=Cohomology(HomToIntegersModP(R,2),3);
1
gap> IsSuperperfect(PerfectGroup(120,1));
true
gap> CoefficientsOfUnivariateRationalFunction(PoincareSeries(SmallGroup(16,13)));
[ [ 1, 1, 1 ], [ 1, -2, 2, -2, 1 ], 0 ]
gap> STOP_TEST( "tst.tst", 1000 );


