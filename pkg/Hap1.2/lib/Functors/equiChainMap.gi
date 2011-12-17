#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(EquivariantChainMap,
function(R,S,f)
local 
	HomotopyS, EltsQ, 
	DimensionR,BoundaryR, EltsG, Mult,
	GhomQ, 			#Let f:G-->Q
	Charact,
	map, mapgens, ChainMap, mapgensRec, 
	N,m,i,g;

N:=Minimum(EvaluateProperty(R,"length"),EvaluateProperty(S,"length"));
HomotopyS:=S.homotopy;
EltsQ:=S.elts;
	for g in S.group do
	if not g in EltsQ then Append(EltsQ,[g]);fi;
	od;
DimensionR:=R.dimension;
BoundaryR:=R.boundary;
EltsG:=R.elts;
#	for g in R.group do
#        if not g in EltsG then Append(EltsG,[g]);fi;
#        od;


mapgensRec:=[];
for m in [0..N] do
mapgensRec[m+1]:=[];
for i in [1..DimensionR(m)] do
mapgensRec[m+1][i]:=[];
for g in [1..Length(R.elts)] do
mapgensRec[m+1][i][g]:=0;
od;
od;
od;

#####################################################################
GhomQ:=function(i);
return Position(EltsQ,Image(f,EltsG[i]));
end;
#####################################################################

#####################################################################
Mult:=function(i,j);
return Position(EltsQ,EltsQ[i]*EltsQ[j]);
end;
#####################################################################

#####################################################################
mapgens:=function(x,m)
local z,u,a,y;

if mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]=0 then

if x[2]>1 then
y:=ShallowCopy(mapgens([x[1],1],m));
Apply(y,b->[b[1],Mult(GhomQ(x[2]),b[2])]);

	if x[1]>0 then mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]:=y;
	else
	mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]:=NegateWord(y);
	fi;
return y;
fi;
   if m=0 then 
   u:=[[SignInt(x[1]),GhomQ(x[2])]]; 
      if x[1]>0 then
      mapgensRec[m+1][x[1]][x[2]]:=u;
      else
      mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]:=NegateWord(u);
      fi;
   return u;
   fi;

   if m>0 then y:=StructuralCopy(BoundaryR(m,x[1]));
   z:=map(y,m-1);
   u:=[];
      for a in z do
            u:=AddWords(HomotopyS(m-1,a),u);
      od;
      Apply(u,t->[t[1],Mult(GhomQ(x[2]),t[2])]);
      if x[1]>0 then
            mapgensRec[m+1][x[1]][x[2]]:=u;
      else
      mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]:=NegateWord(u);
      fi;
      return  u;
   fi;

else if x[1]>0 then return mapgensRec[m+1][AbsoluteValue(x[1])][x[2]];
else return NegateWord(mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]); fi;
fi;

end;
#####################################################################

#####################################################################
map:=function(w,m)
local a, u,v,x,y,z;
v:=[];

   for x in w do
   v:=AddWords(mapgens(x,m),v);
   od;

return AlgebraicReduction(v);
end;
#####################################################################

Charact:=Maximum(EvaluateProperty(R,"characteristic"),
                 EvaluateProperty(R,"characteristic"));

if EvaluateProperty(R,"characteristic")>0
   and
   EvaluateProperty(R,"characteristic")>0
   and
   not EvaluateProperty(R,"characteristic")=
         EvaluateProperty(R,"characteristic")
then Charact:=1;
fi;

return rec(
	    source:=R,
	    target:=S,
	    mapping:=map,
	    properties:=
	    [["type","equivariantChainMap"],
	     ["characteristic",Charact]  ]);
end);
#####################################################################

