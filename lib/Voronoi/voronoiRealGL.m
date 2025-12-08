AttachSpec("~/pkg/Voronoi/VMH-TotallyRealNumberFields/Magma/VMH-spec");
V := VoronoiAlgorithm();
ComputeComplexGL("/tmp/HAP_voronoi.g", V);
quit;
