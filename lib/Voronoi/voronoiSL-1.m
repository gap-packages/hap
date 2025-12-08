AttachSpec("~/pkg/Voronoi/VMH-ImaginaryQuadraticNumberFields/Magma/VMH-spec");
V := VoronoiAlgorithm();
index:=4;
CheckMembership := func<x| IsInSL(x)>;
Reps:=SystemOfRepresentativesFiniteIndex(V`MultFreeList,CheckMembership,index);
ComputeComplexLowIndexSubgroup("/tmp/HAP_voronoi.g", Reps, CheckMembership , V);
quit;
