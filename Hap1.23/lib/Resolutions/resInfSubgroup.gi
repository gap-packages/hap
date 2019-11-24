#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionSubgroup,
function(R,gensK)
local 
		DimensionR, BoundaryR, HomotopyR, EltsG,
		Dimension, Boundary, Homotopy, 
		EltsK, K, sK, TransK,
		F, GhomF, EltsF, EltsKF, i, j,  
		Gword2Kword, G2K, Pair2Int, Int2Pair,
		Mult,N,G;

N:=EvaluateProperty(R,"length");
DimensionR:=R!.dimension;
BoundaryR:=R!.boundary;
HomotopyR:=R!.homotopy;
EltsG:=R!.elts;
G:=R!.group;
F:=FreeGroup(Length(GeneratorsOfGroup(G)));
GhomF:=GroupHomomorphismByImagesNC(G,F,GeneratorsOfGroup(G),
						GeneratorsOfGroup(F));

	#############################################################
	#Let's make sure all boundaries in R have been computed.
	for i in [1..N] do
	for j in [1..DimensionR(i)] do
	BoundaryR(i,j);
	od; od;
	#############################################################

EltsF:=List(EltsG,x->Image(GhomF,x));

EltsK:=[];
EltsKF:=[];
if IsList(gensK) then K:=Group(gensK); else K:=gensK; fi;
TransK:=RightTransversal(G,K);
sK:=Size(TransK);

#####################################################################
Mult:=function(i,j)
local r,x;
x:=Image(GhomF,TransK[i]*EltsG[j]);
r:=Position(EltsF,x);
if not r=fail then return r;
else
Append(EltsF,[x]);
Append(EltsG,[TransK[i]*EltsG[j]]);
return Length(EltsG);
fi;
end;
#####################################################################

#####################################################################
Dimension:=function(n);
return sK*DimensionR(n);
end;
#####################################################################

#####################################################################
G2K:=function(g)
local t,k,r,x;
t:=PositionCanonical(TransK,EltsG[g]);
x:=Image(GhomF,EltsG[g]*TransK[t]^-1);
r:=Position(EltsKF,x);
if not r=fail then k:=r;
else
Append(EltsKF,[x]);
Append(EltsK,[EltsG[g]*TransK[t]^-1]);
k:=Length(EltsK);				
fi;
return [k,t];
end;
#####################################################################

#####################################################################
Pair2Int:=function(x)
local i,t;
i:=x[1]; t:=x[2];
return SignInt(i)*((AbsoluteValue(i)-1)*sK + t);
end;
#####################################################################

#####################################################################
Int2Pair:=function(i)
local j,k, x;
j:=AbsoluteValue(i);
x:=j mod sK;
k:=(j-x)/sK;
if not x=0 then return [SignInt(i)*(k+1),x]; else
return [SignInt(i)*k,sK]; fi;
end;
#####################################################################

#####################################################################
Gword2Kword:=function(w)
local x, y, v;

v:=[];
for x in w do
y:=G2K(x[2]);
y:=[Pair2Int([x[1],y[2]]),y[1]];
Append(v,[y]);
od;
return v;
end;
#####################################################################

#####################################################################
Boundary:=function(n,i)
local x, w;
x:=Int2Pair(i);
w:=StructuralCopy(BoundaryR(n,x[1]));
Apply(w, y->[y[1],Mult(x[2],y[2])]);
return Gword2Kword(w);
end;
#####################################################################


return Objectify(HapResolution,
	  rec(
	   dimension:=Dimension,
	   boundary:=Boundary,
	   homotopy:=fail,
	   elts:=EltsK,
	   group:=K,
	   properties:=
	   [["length",EvaluateProperty(R,"length")],
	    ["characteristic",EvaluateProperty(R,"characteristic")],
	    ["type","resolution"],
	    ["reduced",false] ]));
end);
#####################################################################


