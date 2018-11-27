#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ModularEquivariantChainMap,
function(R,S,f)
local 
	HomotopyS, EltsQ, 
	DimensionR,BoundaryR, EltsG, Mult,
	GhomQ, 			#Let f:G-->Q
	GhomQlst,
	Charact,
	map, mapgens, ChainMap, mapgensRec, 
	Multmat,
	N,m,i,j,g;

if (not "solutionMatBoundaryMatrices" in NamesOfComponents(S) )
and (not "solutionMatBoundaryMatrices" in NamesOfComponents(S)  )
then 
Print("This function can only be applied to resolutions constructed using ResolutionPrimePowerGroup().\n");
return fail;
fi;

if not
EvaluateProperty(R,"characteristic")= EvaluateProperty(S,"characteristic")
then
Print("This function must be applied to two resolutions of equal characteristic.\n");
return fail; 
fi;

Charact:=EvaluateProperty(R,"characteristic");

N:=Minimum(EvaluateProperty(R,"length"),EvaluateProperty(S,"length"));
HomotopyS:=S!.homotopy;
EltsQ:=S!.elts;
DimensionR:=R!.dimension;
BoundaryR:=R!.boundary;
EltsG:=R!.elts;

mapgensRec:=[];
for m in [0..N] do
mapgensRec[m+1]:=[];
for i in [1..DimensionR(m)] do
mapgensRec[m+1][i]:=[];
for g in [1..Length(R!.elts)] do
mapgensRec[m+1][i][g]:=0;
od;
od;
od;

#####################################################################
GhomQ:=function(i);
return Position(EltsQ,Image(f,EltsG[i]));
end;
#####################################################################
GhomQlst:=List([1..Order(R!.group)],GhomQ);
#####################################################################
GhomQ:=function(i);
return GhomQlst[i];
end;
#####################################################################

#####################################################################
Mult:=function(i,j);
return Position(EltsQ,EltsQ[i]*EltsQ[j]);
end;
#####################################################################
if Order(S!.group)<1000 then
Multmat:=[];
for i in [1..Order(S!.group)] do
Multmat[i]:=[];
for j in [1..Order(S!.group)] do
Multmat[i][j]:=Mult(i,j);
od;
od;
#####################################################################
Mult:=function(i,j);
return Multmat[i][j];
end;
#####################################################################
fi;


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
##########################
      for a in Collected(z) do
      Append(u,MultiplyWord(a[2] mod Charact,
List(HomotopyS(m-1,a[1]), t->[t[1],Mult(GhomQ(x[2]),t[2])])
));
      od;
#########################

u:=AlgebraicReduction(u,Charact);

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

#v:=Concatenation(List(w,x->mapgens(x,m)));

v:=Collected(w);
Apply(v,x->MultiplyWord(x[2] mod Charact,  mapgens(x[1],m)));
v:= Concatenation(v);

return AlgebraicReduction(v,Charact);
end;
#####################################################################


return Objectify(HapEquivariantChainMap,
	   rec(
	    source:=R,
	    target:=S,
	    mapping:=map,
	    properties:=
	    [["type","equivariantChainMap"],
	     ["characteristic",Charact]  ]));
end);
#####################################################################

