
#############################################################
#############################################################
InstallGlobalFunction("ResolutionAffineCrystGroup",
function(G, V)
local R, rG, lG, rgens, lgens;

if not (IsAffineCrystGroupOnLeft(G) or IsAffineCrystGroupOnRight(G)) then
Print("This function must be applied to an affine crystallographic group.\n");
return fail;
fi;

if  IsAffineCrystGroupOnRight(G) then
rG:=G; rgens:=GeneratorsOfGroup(G);
lgens:=List(rgens,a->TransposedMat(a));
lG:=AffineCrystGroupOnLeft(lgens);
fi;

if  IsAffineCrystGroupOnLeft(G) then
lG:=G; lgens:=GeneratorsOfGroup(G);
rgens:=List(lgens,a->TransposedMat(a));
rG:=AffineCrystGroupOnRight(rgens);
fi;
 

R:=ResolutionBieberbachGroup(rG,V);
Apply(R!.elts,TransposedMat);
R!.group:=lG;
return R;

end);
#############################################################
#############################################################

#############################################################
#############################################################
InstallGlobalFunction("EquivariantEuclideanSpace",
function(G,V)
local R,Y,dim, Stabilizer, T;

R:=ResolutionAffineCrystGroup(G,V);

T:=Group(Identity(G));
##########
Stabilizer:=function(n,k);
return T;
end;
##########

dim:=Position(List([0..Length(R)],R!.dimension),0)-2;

return Objectify(HapEquivariantCWComplex,
            rec(
            dimension:=R!.dimension,
            boundary:=R!.boundary,
            elts:=R!.elts,
            group:=R!.group,
            stabilizer:=Stabilizer,
            properties:=
            [["dimension",dim],
            ]  ));




end);
#############################################################
#############################################################

