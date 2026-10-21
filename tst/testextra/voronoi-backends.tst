gap> START_TEST("voronoi-backends.tst");
gap> RegisterVoronoiBackend("hap-test", request -> request);
gap> haprequest := rec(group:="unit",rank:=3,
> field:=rec(polynomial:=[-1,-1,0,1]),order:=rec(basis:=[1,2,3]),
> lattice:=rec(ideals:=["a","b","c"]),extra:=rec(example:=true));;
gap> Print(IsIdenticalObj(haprequest,VoronoiComplex(haprequest : backend := "hap-test")),"\n");
true
gap> haprequest := VoronoiComplexGL(2,-7 : backend := "hap-test");;
gap> Print([haprequest.group,haprequest.rank,haprequest.field.radicand],"\n");
[ "GL", 2, -7 ]
gap> Print([haprequest.field.kind,haprequest.order.kind,haprequest.lattice.kind],"\n");
[ "quadratic", "maximal", "free" ]
gap> haprequest := VoronoiComplexSL(2,-1 : backend := "hap-test");;
gap> Print([haprequest.group,haprequest.rank,haprequest.field.radicand],"\n");
[ "SL", 2, -1 ]
gap> haprequest := VoronoiComplexGL(3,5 : backend := "hap-test");;
gap> Print([haprequest.group,haprequest.rank,haprequest.field.radicand],"\n");
[ "GL", 3, 5 ]
gap> RegisterVoronoiBackend("hap-test-fail", request -> fail);
gap> Print(VoronoiComplex(rec() : backend := "hap-test-fail"),"\n");
fail
gap> Print(VoronoiComplexSL(2,5 : backend := "hap-test-fail"),"\n");
fail
gap> RegisterVoronoiBackend("hap-test", request -> request.rank);
gap> Print(VoronoiComplexGL(4,5 : backend := "hap-test"),"\n");
4
gap> Unbind(haprequest);
gap> Unbind(HAP_VORONOI_BACKENDS.("hap-test"));
gap> Unbind(HAP_VORONOI_BACKENDS.("hap-test-fail"));
gap> STOP_TEST("voronoi-backends.tst", 1000);
