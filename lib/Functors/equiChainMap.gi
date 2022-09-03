#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(EquivariantChainMap,
function(R,S,f)
local 
	HomotopyS, EltsQ, 
	DimensionR,BoundaryR, EltsG, Mult,
	GhomQ, 			#Let f:G-->Q
	GhomQlst,
	Charact,
	map, mapgens, ChainMap, mapgensRec, 
	QisFinite,
	Multmat,
        IDEL,
	N,m,i,j,g,AlgRed;

N:=Minimum(EvaluateProperty(R,"length"),EvaluateProperty(S,"length"));
HomotopyS:=S!.homotopy;
EltsQ:=S!.elts;
QisFinite:=false;
if IsFinite(S!.group) then
	if Order(S!.group)=Length(S!.elts) then QisFinite:=true; fi;
fi;
if QisFinite then
	for g in S!.group do
	if not g in EltsQ then Add(EltsQ,g);fi;
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
od;
od;

if QisFinite then
#####################################################################
GhomQ:=function(i);
#return Position(EltsQ,Image(f,EltsG[i]));
return Position(EltsQ,ImageElm(f,EltsG[i]));   #Added Jan 2012
end;
#####################################################################
if IsFinite(R!.group) then
GhomQlst:=List([1..Order(R!.group)],GhomQ);
#####################################################################
GhomQ:=function(i);
return GhomQlst[i];
end;
#####################################################################
fi;
else
#####################################################################
GhomQ:=function(i)
local p,Eltq;
#Eltq:=Image(f,EltsG[i]);
Eltq:=ImageElm(f,EltsG[i]);   #Added Jan 2012
p:= Position(EltsQ,Eltq);
if p=fail then Add(EltsQ,Eltq); ##Changed!
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
if N*Order(S!.group)<200  then
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
else
#####################################################################
Mult:=function(i,j)
local p,Eltq;
Eltq:=EltsQ[i]*EltsQ[j];
p:= Position(EltsQ,Eltq);
if p=fail then Add(EltsQ,Eltq); ##changed!
p:=Length(EltsQ); fi;
return p;
end;
#####################################################################
fi;

#if IsBound(EltsQ!.mult) then Mult:=EltsQ!.mult; fi;

Charact:=Maximum(EvaluateProperty(R,"characteristic"),
                 EvaluateProperty(S,"characteristic"));

if not IsPrime(Charact) then AlgRed:=AlgebraicReduction;
else
AlgRed:=function(v); return AlgebraicReduction(v,Charact); end;
fi;


IDEL:=Position(EltsG,Identity(R!.group));   ##8/02/2012
if IDEL=fail then Add(EltsG,Identity(R!.group)); IDEL:=Length(EltsG);fi;
#####################################################################
mapgens:=function(x,m)
local z,u,a,y;

if not IsBound(mapgensRec[m+1][AbsoluteValue(x[1])][x[2]]) then

if not x[2]=IDEL then
y:=ShallowCopy(mapgens([x[1],IDEL],m));
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
u:=AlgRed(u);
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

return AlgRed(v);
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
            originalHom:=f,
	    properties:=
	    [["type","equivariantChainMap"],
	     ["characteristic",Charact]  ]));
end);
#####################################################################

