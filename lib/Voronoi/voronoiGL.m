AttachSpec("~/pkg/Voronoi/VMH-ImaginaryQuadraticNumberFields/Magma/VMH-spec");
V := VoronoiAlgorithm();
ComputeComplexGL("/tmp/HAP_voronoi.g", V);
quit;
