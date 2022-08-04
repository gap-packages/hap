gap> START_TEST("HAP library");
gap> G:=SymmetricGroup(4);;
gap> x:=(1,2,3,4,5,6,7,8);;
gap> a:=Group(x^2);;
gap> b:=Group(x);;
gap> ahomb:=GroupHomomorphismByFunction(a,b,y->y);;
gap> A:=TrivialGModuleAsGOuterGroup(G,a);;
gap> B:=TrivialGModuleAsGOuterGroup(G,b);;
gap> phi:=GOuterGroupHomomorphism();;
gap> phi!.Source:=A;;
gap> phi!.Target:=B;;
gap> phi!.Mapping:=ahomb;;
gap> Hphi:=CohomologyHomomorphism(phi,3);;
gap> Size(ImageOfGOuterGroupHomomorphism(Hphi));
8
gap> 
gap> Size(KernelOfGOuterGroupHomomorphism(Hphi));
2
gap> STOP_TEST( "tst.tst", 1000 );

