HAP_KOMPLEKS123:=fail;
DIMS:=0;;   #THIS IS NOT ALLOWED! CHANGE IT!!

# Optional providers implement function(request), returning a HAP non-free
# resolution or fail for unsupported inputs. HAP passes the request record
# unchanged; providers interpret and validate its field, order, lattice, group
# and other data. No extra dependency is loaded by HAP.
BindGlobal("HAP_VORONOI_BACKENDS", rec());
BindGlobal("RegisterVoronoiBackend", function(name, provider)
    if not IsString(name) or Length(name)=0 or name="magma"
       or not IsFunction(provider) then
        Error("expected a nonempty backend name other than magma and a function");
    fi;
    HAP_VORONOI_BACKENDS.(name) := provider;
end);
BindGlobal("CallVoronoiBackend", function(name, request)
    if not IsRecord(request) then
        Error("a Voronoi request must be a record");
    fi;
    if not IsString(name) or not IsBound(HAP_VORONOI_BACKENDS.(name)) then
        Error("unknown Voronoi backend; register it with RegisterVoronoiBackend");
    fi;
    return HAP_VORONOI_BACKENDS.(name)(request);
end);
BindGlobal("VoronoiComplex", function(request)
    local backend;
    backend := ValueOption("backend");
    if backend=fail then
        Error("specify a registered Voronoi backend with : backend := name");
    fi;
    return CallVoronoiBackend(backend,request);
end);
BindGlobal("HAPQuadraticVoronoiRequest", function(kind,n,d)
    return rec(group:=kind,rank:=n,
        field:=rec(kind:="quadratic",radicand:=d),
        order:=rec(kind:="maximal"),lattice:=rec(kind:="free"));
end);

############################################################
############################################################
VoronoiComplexGL:=function(n,d)
local file,file2,m,k,G,backend;

backend:=ValueOption("backend");
if backend<>fail and backend<>"magma" then
return CallVoronoiBackend(backend,HAPQuadraticVoronoiRequest("GL",n,d));
fi;

if d<0 then
file:="~/pkg/Voronoi/VMH-ImaginaryQuadraticNumberFields/Magma/BasicData.m";
PrintTo(file, "n:=", n, "; d:=", d, "; steinitz:=1;\n");
fi;
if d>0 then
file:="~/pkg/Voronoi/VMH-TotallyRealNumberFields/Magma/BasicData.m";
PrintTo(file,"K:=QuadraticField(",d,"); n:=", n, "; V:=KMatrixSpace(K,n,1); L:=[IntegralBasis(K)[i]*V.j: i in [1..Degree(K)],j in [1..n]];"); 
fi;

file2:=DirectoriesPackageLibrary("HAP");
if d<0 then
file2:=Filename(file2,"Voronoi/voronoiGL.gi");
fi;
if d>0 then
file2:=Filename(file2,"Voronoi/voronoiRealGL.gi");
fi;

Read(file2);
G:=HAP_KOMPLEKS123!.group;
G!.bianchiInteger:=d;
SetName(G,Concatenation("GL(2,O",String(d),")"));
HAP_KOMPLEKS123!.group:=G;

for m in [0..Length(DIMS)-1] do
for k in [1..HAP_KOMPLEKS123!.dimension(m)] do
HAP_KOMPLEKS123!.boundary(m,k);
od;
od;
return HAP_KOMPLEKS123;

end;
#############################################################
#############################################################

############################################################
############################################################
VoronoiComplexSL:=function(n,d)
local file,file2,m,k,G,backend;

backend:=ValueOption("backend");
if backend<>fail and backend<>"magma" then
return CallVoronoiBackend(backend,HAPQuadraticVoronoiRequest("SL",n,d));
fi;

if d<0 then
file:="~/pkg/Voronoi/VMH-ImaginaryQuadraticNumberFields/Magma/BasicData.m";
fi;
if d>0 then
Print("For a totally real quadratic number field we'd need the order of the unit group.\n");
return fail; 
file:="~/pkg/Voronoi/VMH-TotallyRealNumberFields/Magma/BasicData.m";
fi;
PrintTo(file, "n:=", n, "; d:=", d, "; steinitz:=1;\n");

file2:=DirectoriesPackageLibrary("HAP");
file2:=Filename(file2,"Voronoi/voronoiSL.gi");

Read(file2);
G:=HAP_KOMPLEKS123!.group;
G!.bianchiInteger:=d;
SetName(G,Concatenation("SL(2,O",String(d),")"));
HAP_KOMPLEKS123!.group:=G;
for m in [0..Length(DIMS)-1] do
for k in [1..HAP_KOMPLEKS123!.dimension(m)] do
HAP_KOMPLEKS123!.boundary(m,k);
od;
od;
return HAP_KOMPLEKS123;

end;
#############################################################
#############################################################

