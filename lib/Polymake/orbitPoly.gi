#####################################################################
#####################################################################
InstallGlobalFunction(OrbitPolytope,
function(G,v,Props)
local Points,poly, tmp, x, a, b, c, d, U, V, W, X;

############################CREATE POINTS###############
if IsPermGroup(G) then
    Points := Orbit(G, v, Permuted);
else
    Points := Orbit(G, v, OnRight);
fi;
######################POINTS CREATED###################


################## CALCULATE HULL POINTS ###############

poly:=CreatePolymakeObject();
AppendPointlistToPolymakeObject(poly,Points);


################# HULL CALCULATED ###################################


if "DIMENSION" in Props or "dimension" in Props then
    tmp := Polymake(poly,"DIM");
    Print("Dimension of orbit polytope is: ", tmp, "\n");
fi;

if "VERTEX_DEGREES" in Props or "vertex_degrees" in Props then
    tmp:=Set(Polymake(poly,"VERTEX_DEGREES"));
    Print("Vertex degree in graph of polytope is: ", tmp[1], "\n");
fi;


if "VISUAL_GRAPH" in Props or "visual_graph" in Props then
    Polymake(poly,"VISUAL_GRAPH");
fi;

if "SCHLEGEL" in Props or "schlegel" in Props then
    Polymake(poly,"SCHLEGEL");
fi;


if "VISUAL" in Props or "visual" in Props then
    if IsPermGroup(G) and Length(v)=4 then
        tmp := ShallowCopy(Points);
        for x in G do
            a:=v[1^x]-v[1];
            b:=v[2^x]-v[2];
            c:=v[3^x]-v[3];
            d:=v[4^x]-v[4];
            U:=2*a-2*b;
            V:=2*c-2*d;
            W:=a+b-c-d;
            Add(tmp, [1,U,V,W]);
        od;
        poly:=CreatePolymakeObject();
        AppendPointlistToPolymakeObject(poly,tmp);
    fi;
    Polymake(poly,"VISUAL");
fi;


end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Display,
"method for displaying convex hulls of sets of points",
[IsPolymakeObject],
function(F)
    local Dim, Points,poly;
    
    Dim:=Length(Polymake(F,"F_VECTOR"));
    Points:=Polymake(F,"VERTICES");
    
    poly:=CreatePolymakeObject();
    AppendPointlistToPolymakeObject(poly,Points);
    
    Polymake(poly,"VISUAL");

end);
#####################################################################
#####################################################################

