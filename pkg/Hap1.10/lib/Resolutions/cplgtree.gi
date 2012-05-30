SL2ZTree:=function(m,p)
local Elts,H,K,G,Gamma,
      Id,ID,Idcoset,
      Boundary,Dimension,Action,Stabilizer,Homotopy,StabGrps,pos;
Elts:=[[[1,0],[0,1]]];
#StabGrps:=[];
if p=0 then
  H:=Group([[0,-1],[1,0]]);
  K:=Group([[0,-1],[1,1]]);
  G:=SL(2,Integers);
  Gamma:=Group([[-1,0],[0,-1]]);
  Append(Elts,Elements(H));
  Append(Elts,Elements(K));
  Append(Elts,Elements(Gamma));
  Elts:=SSortedList(Elts);
else
if m=1 then H:=SL(2,Integers);
else H:=SL2Z(1/m);fi;
K:=SL2Z(p);
G:=SL2Z(1/(m*p));
ID:=Group(One(G));
Gamma:=CongruenceSubgroupGamma0(p);
Append(Elts,GeneratorsOfGroup(H));
Append(Elts,GeneratorsOfGroup(K));
Append(Elts,GeneratorsOfGroup(Gamma));
Elts:=SSortedList(Elts);
fi;
SetName(Gamma,"Gamma");
Id:=Position(Elts,[[1,0],[0,1]]);
#######################
Boundary:=function(n,k)
local w,t;
if not n=1 then return [];fi;
t:=pos(Elts,CanonicalRightCountableCosetElement(Gamma,Elts[Id]^-1)^-1); 
w:=[[1,t],[-2,t]];
if k>0 then return w;
else return NegateWord(w);
fi;
end;
#######################
Dimension:=function(n)
if not n in [0,1] then return 0;fi;
if n=0 then return 2;fi;
if n=1 then return 1;fi;
end;
#######################

#######################
Action:=function(n,k,l);
return 1;
end;
#######################
StabGrps:=[];
Add(StabGrps,[H,K]);
Add(StabGrps,[Gamma]);
#######################
Stabilizer:=function(n,k);
return StabGrps[n+1][k];
end;
#######################
pos:=function(Elts,g)
if Position(Elts,g)=fail then Add(Elts,g);return Position(Elts,g);
else return Position(Elts,g);
fi;
end;
Idcoset:=pos(Elts,CanonicalRightCountableCosetElement(Gamma,Elts[Id]^-1)^-1);
#######################
Homotopy:=function(n,w)
local d,path,i,h,k,g,pk,r,t;
k:=w[1];
g:=w[2];
pk:=AbsInt(k);
d:=SL2ZmElementsDecomposition(Elts[g],p);
if p=0 then
r:=[];
Add(r,d[1]);
for i in [2..Length(d)] do
    if d[i]*d[i-1] in H then
	r[Length(r)]:=r[Length(r)]*d[i];
    else
	Add(r,d[i]);
    fi;
od;
d:=r;
r:=[];
Add(r,d[1]);
for i in [2..Length(d)] do
    if d[i]*d[i-1] in K then
	r[Length(r)]:=r[Length(r)]*d[i];
    else
	Add(r,d[i]);
    fi;
od;
d:=StructuralCopy(r);
fi;
#Print("d=",d,"\n");
if d[1] in K and not d[1] in Gamma then 
r:=[[[1,0],[0,1]]];
Append(r,d);
d:=StructuralCopy(r);
#Print("d=",d,"\n");
fi;
#Print("d=",d,"\n");
h:=[[1,0],[0,1]];
path:=[];
if d[Length(d)] in Gamma and Length(d)>1 then Remove(d,Length(d));fi;

if pk=1 then 
    #if g in Gamma then return [1,Id];fi;
    if d[Length(d)] in H then Remove(d,Length(d));fi;
    for i in [1..Length(d)] do
	#d[i]:=CanonicalRightCosetElement(Gamma,d[i]^-1)^-1;
	h:=h*d[i];
	t:=CanonicalRightCountableCosetElement(Gamma,h^-1)^-1;
	Add(path,[(-1)^(i),pos(Elts,t)]);
	#if h in H and i<Length(d) then path:=[];fi;
    od;
else
    #Print("d=",d,"\n");
    if Elts[g] in Gamma then Add(path,[-1,Idcoset]);fi;
    if d[Length(d)] in K then Remove(d,Length(d));fi;
    for i in [1..Length(d)] do
	#d[i]:=CanonicalRightCosetElement(Gamma,d[i]^-1)^-1;
	h:=h*d[i];
	#Print("h=",h,"\n");
	t:=CanonicalRightCountableCosetElement(Gamma,h^-1)^-1;
	Add(path,[(-1)^(i),pos(Elts,t)]);
	#h:=h*d[i];
	#if h in K and i<Length(d) then path:=[];fi;
    od;
fi;
if k>0 then
return path;
else return NegateWord(path);
fi;
end;
#######################


return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=Homotopy,
            elts:=Elts,
            group:=G,
            stabilizer:=Stabilizer,
            action:=Action,
            properties:=
            [["length",10000],
             ["characteristic",0],
             ["type","resolution"]]  ));

end;
###################################################################
TreeOfResolutionsToSL2Zcomplex:=function(D,G)
local RH,RK,RGamma,
      H,K,Gamma,
      NameH,
      i,j,m,p,
      C,Resolutions,NamesOfGroups;
RH:=D[1];
RK:=D[2];
RGamma:=D[3];
H:=RH!.group;
K:=RK!.group;
Gamma:=RGamma!.group;
NameH:=H!.Name;
if H=SL(2,Integers) then m:=1;
else
i:=Position(NameH,'/');
j:=Position(NameH,']');
m:=Int(NameH{[i+1..j-1]});
fi;
p:=Gamma!.LevelOfCongruenceSubgroup;
C:=SL2ZTree(m,p);
NamesOfGroups:=[Name(H),Name(K),Name(Gamma)];
Resolutions:=[RH,RK,RGamma];
C!.resolutions:=[Resolutions,NamesOfGroups];
return C;
end;
####################################################################