# Run with the ArithmeticVoronoi environment pinned in ../workflows/CI.yml:
# julia --startup-file=no --project=/path/to/arithmetic_voronoi \
#     .github/scripts/voronoi-julia.jl
println("Loading the Julia provider and its dependencies...")
flush(stdout)
package_load_started = time_ns()
using ArithmeticVoronoi, Oscar, Test
import Oscar: GAP
println("Package loading: ", round((time_ns() - package_load_started) / 1e9; digits=2), " s")
flush(stdout)

# Test this HAP checkout, including its startup registration. Check before
# enabling the provider so the adapter's compatibility wrapper cannot mask
# missing entry points in HAP itself.
hap_root = realpath(joinpath(@__DIR__, "..", ".."))
@assert GAP.Packages.load(hap_root)
@assert realpath(GAP.Packages.locate_package("hap")) == hap_root
@assert GAP.Globals.IsBoundGlobal(GAP.GapObj("RegisterVoronoiBackend"))
@assert GAP.Globals.IsBoundGlobal(GAP.GapObj("VoronoiComplex"))
dispatcher = GAP.Globals.VoronoiComplex
println("HAP checkout: ", hap_root)
println("Julia ", VERSION, "; OSCAR ", pkgversion(Oscar), "; GAP ", GAP.GAP_VERSION)

@testset "Julia provider: GL2 over O_(Q(sqrt(-7)))" begin
    enable_hap_julia_backend(; progress=true, resolution_length=3)
    @test GAP.Globals.IsIdenticalObj(dispatcher, GAP.Globals.VoronoiComplex)

    request = Dict(
        "group" => "GL", "rank" => 2,
        "field" => Dict("kind" => "quadratic", "discriminant" => -7),
        "order" => Dict("kind" => "maximal"),
        "lattice" => Dict("kind" => "free"),
    )
    println("Request: ", request)
    flush(stdout)
    construct = GAP.evalstr("request -> VoronoiComplex(request : backend := \"julia\")")
    # A fresh native enumeration, without Magma, stored complexes or a C++ build.
    K = construct(GAP.GapObj(request; recursive=true))
    @test GAP.Globals.IsHapNonFreeResolution(K)
    cells = Int.(GAP.gap_to_julia(
        GAP.evalstr("K -> List([0..2], n -> K!.dimension(n))")(K)))
    @test cells == [1, 2, 1]
    println("Cell orbits in dimensions 0, 1, 2: ", cells)

    R = GAP.Globals.FreeGResolution(K, 3)
    C = GAP.Globals.TensorWithIntegers(R)
    # These invariants were independently checked with Magma/VMH;
    # neither the example nor CI needs Magma to compute them.
    for (degree, expected) in ((1, [2, 2, 2]), (2, [2, 12]))
        invariants = sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(C, degree))))
        println("H", degree, " invariant factors: ", invariants)
        @test invariants == expected
    end
end
