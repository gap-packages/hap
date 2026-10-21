gap> START_TEST("voronoi-backends.tst");
gap> RegisterVoronoiBackend("hap-test", function(kind,n,d)
> return [kind,n,d]; end);
gap> Print(VoronoiComplexGL(2,-7 : backend := "hap-test"),"\n");
[ "GL", 2, -7 ]
gap> Print(VoronoiComplexSL(2,-1 : backend := "hap-test"),"\n");
[ "SL", 2, -1 ]
gap> Print(VoronoiComplexGL(3,5 : backend := "hap-test"),"\n");
[ "GL", 3, 5 ]
gap> RegisterVoronoiBackend("hap-test-fail", function(kind,n,d)
> return fail; end);
gap> Print(VoronoiComplexSL(2,5 : backend := "hap-test-fail"),"\n");
fail
gap> Print(VoronoiComplexGL(2,-7 : backend := "hap-test"),"\n");
[ "GL", 2, -7 ]
gap> Unbind(HAP_VORONOI_BACKENDS.("hap-test"));
gap> Unbind(HAP_VORONOI_BACKENDS.("hap-test-fail"));
gap> STOP_TEST("voronoi-backends.tst", 1000);
