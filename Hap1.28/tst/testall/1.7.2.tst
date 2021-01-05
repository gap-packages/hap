gap> START_TEST("HAP library");
gap> F:=FreeGroup(4);;w:=F.1;;x:=F.2;;y:=F.3;;z:=F.4;;
gap> rels:=[w^8, w*x*w*(x*w*x)^-1, y^2, z*x*(x*z)^-1, z^-1*y*z*y, (x*y*x)^2];;
gap> G:=F/rels;;
gap> N2:=[];;N3:=[];;
gap> for u in GeneratorsOfGroup(G) do
> Add(N2,u^2);
> Add(N3,u^3);
> for v in GeneratorsOfGroup(G) do
> Add(N2,Comm(u,v));
> Add(N3,Comm(u,v));
> od;;od;;
gap> N2:=NormalClosure(G,Group(N2));;
gap> N3:=NormalClosure(G,Group(N3));;
gap> AbelianInvariants(N2);
[ 0, 0, 3, 3, 3, 3, 4 ]
gap> AbelianInvariants(N3);
[ 0, 2, 4 ]
gap> STOP_TEST( "tst.tst", 1000 );


