gap> START_TEST("HAP library");
gap> Read(Concatenation(dir,"data134.txt"));
gap> f:=function(j); return D[50][j]; end;
function( j ) ... end
gap> M:=Mapper(D,f,8,20,5);
Simplicial complex of dimension 1.

gap> Display(GraphOfSimplicialComplex(M));
gap> STOP_TEST( "tst.tst", 1000 );


