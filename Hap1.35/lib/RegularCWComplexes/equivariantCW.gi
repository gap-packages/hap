
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
InstallGlobalFunction("ResolutionToEquivariantCWComplex",
function(R)
local Stabilizer,T,dim,i;

#dim:=Position(List([0..Length(R)],R!.dimension),0);
#if dim=fail then
for i in [0..Length(R)] do
if R!.dimension(i)>0 then dim:=i; fi;
od;
#fi;

T:=Group(Identity(R!.group));
##########
Stabilizer:=function(n,k);
return T;
end;
##########


return Objectify(HapEquivariantCWComplex,
            rec(
            dimension:=R!.dimension,
            boundary:=R!.boundary,
            elts:=R!.elts,
            group:=R!.group,
            stabilizer:=Stabilizer,
            properties:=
            [["dimension",dim]
            ]  ));

end);
#############################################################
#############################################################

#############################################################
#############################################################
InstallGlobalFunction("EquivariantCWComplexToResolution",
function(Y)
local Stabilizer,T,dim,n,i;

for n in [0..EvaluateProperty(Y,"dimension")] do
for i in [1..Y!.dimension(n)] do
#if Order(Y!.stabilizer(n,i))>1 then
#Print("EquivariantCWComplex is not free.\n");
#return fail;
#fi;
od;
od;
dim:=EvaluateProperty(Y,"dimension");


return Objectify(HapResolution,
            rec(
            dimension:=Y!.dimension,
            boundary:=Y!.boundary,
            elts:=Y!.elts,
            group:=Y!.group,
            homotopy:=fail,
            properties:=
            [["length",dim],
             ["characteristic",0],
             ["type","resolution"]
            ]  ));

end);
#############################################################
#############################################################


#############################################################
#############################################################
InstallGlobalFunction("EquivariantEuclideanSpace",
function(G,V)
local R;

R:=ResolutionAffineCrystGroup(G,V);

return ResolutionToEquivariantCWComplex(R);

end);
#############################################################
#############################################################

#############################################################
#############################################################
InstallGlobalFunction("EquivariantOrbitPolytope",
function(arg)
local G,V,R,Y,dim, Stabilizer;

G:=arg[1];
V:=arg[2];
if Length(arg)=2 then
R:=PolytopalComplex(G,V);
else
R:=PolytopalComplex(G,V,arg[3]);
fi;

return ResolutionToEquivariantCWComplex(R);

end);
#############################################################
#############################################################

#############################################################
#############################################################
InstallGlobalFunction(EquivariantTwoComplex,
function(G)
local R,Stabilizer,T,F,S;

if IsFpGroup(G) then
R:=ResolutionAsphericalPresentation(G,2);
else
R:=ResolutionFiniteGroup(G,2);;
fi;

return ResolutionToEquivariantCWComplex(R);

end);
#############################################################
#############################################################

#############################################################
#############################################################
InstallMethod(FundamentalGroupOfQuotient,
"presentation of fundamental group of equivariant CW-complex",
[IsHapEquivariantCWComplex],

function(Y) local R, P;
R:=EquivariantCWComplexToResolution(Y);
P:=PresentationOfResolution(R);
return P.freeGroup/P.relators;
end);
#############################################################
#############################################################

#############################################################
#############################################################
InstallMethod(ChainComplexOfQuotient,
"chain complex of the quotient of an equivariant CW-complex",
[IsHapEquivariantCWComplex],

function(Y) local R, P;
R:=EquivariantCWComplexToResolution(Y);
P:=TensorWithIntegers(R);
return P;
end);
#############################################################
#############################################################


#############################################################
#############################################################
InstallGlobalFunction(RestrictedEquivariantCWComplex,
function(Y,H)
local R, X;

R:=EquivariantCWComplexToResolution(Y);
R:=ResolutionSubgroup(R,H);
return ResolutionToEquivariantCWComplex(R);

end);
#############################################################
#############################################################


