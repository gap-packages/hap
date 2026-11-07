# Run with the ArithmeticVoronoi revision pinned in ../workflows/CI.yml:
# julia --startup-file=no --project=/path/to/ArithmeticVoronoi.jl \
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

@testset "ArithmeticVoronoi: rational and quadratic groups" begin
    enable_hap_julia_backend(;progress=true,resolution_length=3)
    @test GAP.Globals.IsIdenticalObj(dispatcher,GAP.Globals.VoronoiComplex)
    construct=GAP.evalstr("r -> VoronoiComplex(r : backend := \"julia\")")
    # These are regression expectations from independent VMH/Magma, HAP stored
    # resolutions and classical homology. All complexes here are computed fresh.
    examples=(
        (request=voronoi_request(-7),cells=[1,2,1],homology=([2,2,2],[2,12])),
        (request=voronoi_request(5),cells=[2,4,5,4,1],homology=([0,2],[2,2,2,2])),
        (request=voronoi_request(-7;group=:SL),cells=nothing,homology=([0,4],[12])),
        (request=voronoi_request(-7;group=:PSL),cells=nothing,homology=([0,2],)),
        (request=voronoi_request(;rank=3),cells=[1,1,2,1],homology=([2],[2,2])),
        (request=voronoi_request(;rank=3,group=:SL),cells=nothing,homology=(Int[],)),
        (request=voronoi_request(-4;rank=3),cells=[1,1,3,5,4,3,2],homology=([4],)),
    )
    first_complex=nothing
    for (i,example) in enumerate(examples)
        request=example.request
        @testset "$(request)" begin
            println("Request: ",request);flush(stdout)
            started=time_ns()
            K=construct(GAP.GapObj(request;recursive=true))
            i==1 && (first_complex=K)
            @test GAP.Globals.IsHapNonFreeResolution(K)
            @test Int(GAP.evalstr("K -> K!.arithmeticVoronoi")(K).rank)==request["rank"]
            @test GAP.evalstr("K -> K!.arithmeticVoronoi")(K).groupKind==GAP.GapObj(request["group"])
            if !isnothing(example.cells)
                cells=GAP.evalstr("{K,d} -> List([0..d-1], K!.dimension)")(K,length(example.cells))
                @test Int.(GAP.gap_to_julia(cells))==example.cells
            end
            R=GAP.Globals.FreeGResolution(K,length(example.homology)+1)
            C=GAP.Globals.TensorWithIntegers(R)
            for (degree,expected) in enumerate(example.homology)
                invariants=sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(C,degree))))
                println("H",degree," invariant factors: ",invariants)
                @test invariants==expected
            end
            println("Example time: ",round((time_ns()-started)/1e9;digits=2)," s");flush(stdout)
        end
    end
    # The direct Julia interface requires no record conversion or dispatcher.
    K=hap_voronoi_complex(;group=:PSL,resolution_length=2)
    R=GAP.Globals.FreeGResolution(K,2)
    @test Int.(GAP.gap_to_julia(GAP.Globals.Homology(GAP.Globals.TensorWithIntegers(R),1)))==[6]
    # Later computations must not overwrite an earlier object's captured data.
    R=GAP.Globals.FreeGResolution(first_complex,3)
    @test sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(GAP.Globals.TensorWithIntegers(R),2))))==[2,12]
end
