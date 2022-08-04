#D
gap> START_TEST("HAP library");
gap> ReadPackage("HAP","tst/testall1/data134.txt");;
gap> f:=function(j); return D[50][j]; end;
function( j ) ... end
gap> M:=Mapper(D,f,8,20,5);
Simplicial complex of dimension 1.

gap> #Display(GraphOfSimplicialComplex(M));
gap> GraphOfSimplicialComplex(M);
Graph on 13 vertices.

gap> STOP_TEST( "tst.tst", 1000 );


