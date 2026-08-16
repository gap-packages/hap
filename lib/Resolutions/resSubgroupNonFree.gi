#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionFiniteSubgroup_NonFree,
function(arg)
local 
		R,
		DimensionR, BoundaryR, HomotopyR, EltsG,
		Dimension, Boundary, BoundaryRec, Homotopy, EltsK,
		G, K, TransK, sK, 
		Gword2Kword, G2K, G2KRec, Pair2Int, Int2Pair,
		Mult, MultRec, i,HomotopyRec;
		
Print("Warning: this function is incomplete and will return wrong answers!\n");

R:=arg[1]; G:=R!.group; K:=arg[2];   #Assuming K<=G
DimensionR:=R!.dimension;
BoundaryR:=R!.boundary;
HomotopyR:=R!.homotopy;
EltsG:=R!.elts;

EltsK:=[];

TransK:=RightTransversal(G,K);
sK:=Size(TransK);

MultRec:=List([1..Length(TransK)],i->[]);
#####################################################################
Mult:=function(i,j) local x,r;
if not IsBound(MultRec[i][j]) then

x:=TransK[i]*EltsG[j];
#r:=Position(EltsG,x);
#if r=fail then 
Add(EltsG,x); r:=Length(EltsG); 
#fi;

MultRec[i][j]:= r;


return r;
fi;

return MultRec[i][j];
end;
#####################################################################

#####################################################################
Dimension:=function(n)
return sK*DimensionR(n);
end;
#####################################################################

G2KRec:=[];
#####################################################################
G2K:=function(g)
local t,k,x;
if not IsBound(G2KRec[g]) then 

t:=PositionCanonical(TransK,EltsG[g]);

x:=EltsG[g]*TransK[t]^-1;

#k:=Position(EltsK,x); 

#if k=fail then 
Add(EltsK,x); k:=Length(EltsK); 
#fi;    

G2KRec[g]:= [k,t];

fi;
return G2KRec[g];
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
Add(v,y);
od;
return v;
end;
#####################################################################

BoundaryRec:=[];
for i in [1..EvaluateProperty(R,"length")] do
BoundaryRec[i]:=[];
od;

#####################################################################
Boundary:=function(n,ii)
local x, w, i;

if n<=0 then return []; fi;

i:=AbsInt(ii);
if not IsBound(BoundaryRec[n][i]) then
x:=Int2Pair(i);
w:=StructuralCopy(BoundaryR(n,x[1]));
Apply(w, y->[y[1],Mult(x[2],y[2])]);
BoundaryRec[n][i]:= Gword2Kword(w);
###BoundaryRec[n][i]:=AlgebraicReduction(BoundaryRec[n][i]);
fi;

if ii>0 then return BoundaryRec[n][i]; 
else return NegateWord(BoundaryRec[n][i]); fi;
end;
#####################################################################

HomotopyRec:=List([0..Length(R)],i->List([1..Dimension(i)],j->[]));;

######################
Homotopy:=function(n,e)
local x,g,pos,ae;

ae:=AbsInt(e[1]);
if not IsBound(HomotopyRec[n+1][ae][e[2]]) then

x:=Int2Pair(ae);
g:=EltsK[e[2]]*TransK[x[2]]; #Need to check the maths again here!


pos:=Position(EltsG,g);


if pos =fail then Add(EltsG,g); pos:=Length(EltsG);fi;
HomotopyRec[n+1][ae][e[2]]:= Gword2Kword( R!.homotopy(n,[x[1],pos]));
fi;

if e[1]>0 then return HomotopyRec[n+1][ae][e[2]];
else return NegateWord(HomotopyRec[n+1][ae][e[2]]);
fi;

end;
######################

return Objectify(HapResolution,
	     rec(
	     dimension:=Dimension,
	     boundary:=Boundary,
	     homotopy:=Homotopy,
	     elts:=EltsK,
	     group:=K,
             Int2Pair:=Int2Pair,
             transversal:=TransK,
	     properties:=
	     [["length",EvaluateProperty(R,"length")],
	      ["characteristic",EvaluateProperty(R,"characteristic")],
	      ["reduced",false],
	      ["type","resolution"] ]));
end);
#####################################################################


