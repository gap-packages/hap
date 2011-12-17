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
	QisFinite,
	N,m,i,g;

N:=Minimum(EvaluateProperty(R,"length"),EvaluateProperty(S,"length"));
HomotopyS:=S!.homotopy;
EltsQ:=S!.elts;
QisFinite:=false;
if IsFinite(S!.group) then
	if Order(S!.group)=Length(S!.elts) then QisFinite:=true; fi;
fi;
if QisFinite then
	for g in S!.group do
	if not g in EltsQ then Append(EltsQ,[g]);fi;
	od;
fi;
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

if QisFinite then
#####################################################################
GhomQ:=function(i);
return Position(EltsQ,Image(f,EltsG[i]));
end;
#####################################################################
else
#####################################################################
GhomQ:=function(i)
local p,Eltq;
Eltq:=Image(f,EltsG[i]);
p:= Position(EltsQ,Eltq);
if p=fail then Append(EltsQ,Eltq);
p:=Length(EltsQ); fi;
return p;
end;
#####################################################################
fi;

if QisFinite then
#####################################################################
Mult:=function(i,j);
return Position(EltsQ,EltsQ[i]*EltsQ[j]);
end;
#####################################################################
else
#####################################################################
Mult:=function(i,j)
local p,Eltq;
Eltq:=EltsQ[i]*EltsQ[j];
p:= Position(EltsQ,Eltq);
if p=fail then Append(EltsQ,Eltq);
p:=Length(EltsQ); fi;
return p;
end;
#####################################################################
fi;

Charact:=Maximum(EvaluateProperty(R,"characteristic"),
                 EvaluateProperty(S,"characteristic"));

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
if Charact=0 then 
 for a in Collected(z) do
      Append(u,MultiplyWord(a[2],
List(HomotopyS(m-1,a[1]), t->[t[1],Mult(GhomQ(x[2]),t[2])])
));
      od;
else
      for a in Collected(z) do
      Append(u,MultiplyWord(a[2] mod Charact,
List(HomotopyS(m-1,a[1]), t->[t[1],Mult(GhomQ(x[2]),t[2])])
));
      od;
fi;
#########################
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

v:=Collected(w);
if Charact=0 then
Apply(v,x->MultiplyWord(x[2],  mapgens(x[1],m)));
else
Apply(v,x->MultiplyWord(x[2] mod Charact,  mapgens(x[1],m)));
fi;
v:= Concatenation(v);
return v;
end;
#####################################################################


if EvaluateProperty(R,"characteristic")>0
   and
   EvaluateProperty(S,"characteristic")>0
   and
   not EvaluateProperty(R,"characteristic")=
         EvaluateProperty(S,"characteristic")
then Charact:=1;
fi;

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

