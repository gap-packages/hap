#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionFiniteSubgroup,
function(arg)
local 
		R,gensG,gensK,
		DimensionR, BoundaryR, HomotopyR, EltsG,
		Dimension, Boundary, BoundaryRec, Homotopy, EltsK,
		G, K, TransK, sK, 
		Gword2Kword, G2K, G2KRec, Pair2Int, Int2Pair,
		Mult, MultRec, FIN, i;
		
if Length(arg)=3 then
R:=arg[1]; gensG:=arg[2]; gensK:=arg[3];
fi;
if Length(arg)=2 then
R:=arg[1]; gensG:=R!.group; gensK:=arg[2];
fi;
				#gensG and gensK originally had to be
				#generating sets. Later I allowed them to be
				#the groups G, K themselves. Sloppy!
DimensionR:=R!.dimension;
BoundaryR:=R!.boundary;
HomotopyR:=R!.homotopy;
EltsG:=R!.elts;


if IsList(gensG) then G:=Group(gensG); else G:=gensG; fi;
if IsList(gensK) then K:=Group(gensK); else K:=gensK; fi;
if Size(G)=Length(EltsG) then FIN:=true; else FIN:= false; fi;

if FIN then 
   if IsPseudoList(EltsG) then EltsK:=EltsG;
   else
   EltsK:=Elements(K);
   fi;
else
EltsK:=[];
fi;


TransK:=RightTransversal(G,K);
sK:=Size(TransK);

if FIN then
#####################################################################
Mult:=function(i,j);
return Position(EltsG,TransK[i]*EltsG[j]);
end;
#####################################################################
else
MultRec:=List([1..Length(TransK)],i->[]);
#####################################################################
Mult:=function(i,j) local x,r;
if not IsBound(MultRec[i][j]) then

x:=TransK[i]*EltsG[j];

r:=Position(EltsG,x);

if r=fail then Add(EltsG,x); r:=Length(EltsG); fi;
return r;
MultRec[i][j]:= r;
fi;

return MultRec[i][j];
end;
#####################################################################
fi;

#####################################################################
Dimension:=function(n);
return sK*DimensionR(n);
end;
#####################################################################

if FIN then
#####################################################################
G2K:=function(g)
local t,k;
t:=PositionCanonical(TransK,EltsG[g]);
k:=Position(EltsK,EltsG[g]*TransK[t]^-1);
return [k,t];
end;
#####################################################################
else
G2KRec:=[];
#####################################################################
G2K:=function(g)
local t,k,x;
if not IsBound(G2KRec[g]) then 

t:=PositionCanonical(TransK,EltsG[g]);
x:=EltsG[g]*TransK[t]^-1;

k:=Position(EltsK,x);

if k=fail then Add(EltsK,x); k:=Length(EltsK); fi;

G2KRec[g]:= [k,t];

fi;
return G2KRec[g];
end;
#####################################################################
fi;

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
#Apply(w, y->[y[1],Position(EltsG,TransK[x[2]]*EltsG[y[2]])   ]); #Changed this back but forgot why this line was ever here!!
BoundaryRec[n][i]:= Gword2Kword(w);
###BoundaryRec[n][i]:=AlgebraicReduction(BoundaryRec[n][i]);
fi;

if ii>0 then return BoundaryRec[n][i]; 
else return NegateWord(BoundaryRec[n][i]); fi;
end;
#####################################################################


######################
Homotopy:=function(n,e)
local x,g,pos;

x:=Int2Pair(e[1]);
g:=EltsK[e[2]]*TransK[x[2]]; #Need to chack the maths again here!


pos:=Position(EltsG,g);


if pos =fail then Add(EltsG,g); pos:=Length(EltsG);fi;
return Gword2Kword( R!.homotopy(n,[x[1],pos]));

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


