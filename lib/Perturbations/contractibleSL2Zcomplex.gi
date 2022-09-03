#(C) This function was written by Bui Anh Tuan

################################################
################################################
InstallGlobalFunction(ContractibleSL2ZComplex,
function()
local
	C, 
        G, StabilizerGroups, Stabilizer,
        lnth,
        dims,Dimension,
        Boundary,
        boundaryList,
        Elts,
	Rot,Stab,
        RotSubGroups,Action, ActionRecord,
        TransMat,
        St0,St1, x, n,k,s,BI,SGN,tmp, LstEl , 
        bool, name,
	GeneratorsRepresentation,SimplifyGeneratorsRepresentation,EdgeFinder,
	pos,AddList,Chomotopy,Path2Gindex,FinalHomotopy,Homotopy,Sign,
        EdgeFinder1,RefineEdge,Edge;

bool:=ReadPackage("HAP","lib/Perturbations/Gcomplexes/SL2Z");

if bool = false then
Print("Complex failed to load.\n");
return fail;
fi;

if HAP_GCOMPLEX_SETUP[1] then 
TransMat:=function(x); return x^-1; end;
else
TransMat:=function(x); return x; end;
fi;


C:=StructuralCopy(HAP_GCOMPLEX_LIST);
lnth:=Length(C)-1;

dims:=List([1..lnth+1],n->Length(C[n]));

###################
Dimension:=function(n);
if n>lnth then return 0; fi;
return dims[n+1];
end;
###################

Elts:=[[[1,0],[0,1]]];
StabilizerGroups:=[];
RotSubGroups:=[];
boundaryList:=[];


#######
for n in [1..lnth+1] do
boundaryList[n]:=[];
StabilizerGroups[n]:=[];
RotSubGroups[n]:=[];
  for k in [1..Dimension(n-1)] do
  Append(Elts,Elements(C[n][k].TheMatrixStab));
  Add(StabilizerGroups[n],C[n][k].TheMatrixStab);
  Add(RotSubGroups[n],C[n][k].TheRotSubgroup);
  od;
od;
####

Elts:=SSortedList(Elts);

#######
for n in [1..lnth+1] do
boundaryList[n]:=[];
for k in [1..Dimension(n-1)] do
tmp:=C[n][k].BoundaryImage;
BI:=tmp.ListIFace;
SGN:=tmp.ListSign;
LstEl:=List(tmp.ListElt,w->TransMat(w));
Append(Elts,Difference(LstEl,Elts));
for s in [1..Length(BI)] do
BI[s]:=[SGN[s]*BI[s],Position(Elts,LstEl[s])];
od;
Add(boundaryList[n],BI);
od;
od;
####

ActionRecord:=[];
for n in [1..lnth+1] do
ActionRecord[n]:=[];
for k in [1..Dimension(n-1)] do
ActionRecord[n][k]:=[];
od;
od;

G:=Group(Elts);

St0:=StabilizerGroups[1][1];
St1:=StabilizerGroups[2][1];
####################
Boundary:=function(n,k)
local b;
if k>0 then
b:=boundaryList[n+1][k];
#Print("b",b,"\n\n");
Apply(b,x->[x[1],  pos(CanonicalRightCosetElement(St0,Elts[x[2]]^-1)^-1)] );
#Print("b1",b,"\n\n");
return b;
#return boundaryList[n+1][k];
else
b:=boundaryList[n+1][-k];
Apply(b,x->[x[1],  pos(CanonicalRightCosetElement(St0,Elts[x[2]]^-1)^-1)] );
return NegateWord(b);
#return NegateWord(boundaryList[n+1][-k]);
fi;
end;
####################

####################
Stabilizer:=function(n,k);
return StabilizerGroups[n+1][k];
end;
####################



####################
Action:=function(n,k,g)
local id,r,u,H,abk,ans;

abk:=AbsInt(k);

if not IsBound(ActionRecord[n+1][abk][g]) then
H:=StabilizerGroups[n+1][abk];

if Order(H)=infinity then ActionRecord[n+1][abk][g]:=1;
#So we are assuming that any infinite stabilizer group acts trivially!!
else
######
id:=CanonicalRightCosetElement(H,Identity(H));
r:=CanonicalRightCosetElement(H,Elts[g]^-1);
r:=id^-1*r;
u:=r*Elts[g];

if u in RotSubGroups[n+1][abk] then  ans:= 1;
else ans:= -1; fi;

ActionRecord[n+1][abk][g]:=ans;
fi;
######
fi;

return ActionRecord[n+1][abk][g];
end;

####################SL2Z-homotopy########################################
GeneratorsRepresentation:=function(g)
local
        S,T,
        q,i,j,
        index;
#g:=CanonicalRightCosetElement(St0,g^-1)^-1;
if Determinant(g)<>1 then return "input is not in SL(2,Z)";
else
S:=[[0,-1],[1,0]];
T:=[[1,1],[0,1]];
i:=0;
index:=[];
while g[2][1]<>0 do
if g[1][1]*g[2][1]<0 and IsInt(g[1][1]/g[2][1])=false then
        q:=Int(g[1][1]/g[2][1])-1;
else q:=Int(g[1][1]/g[2][1]);
fi;
i:=i+1;
index[i]:=q;
i:=i+1;
index[i]:=1;
g:=S*((T^(-q))*g);
od;
if g[1][1]=1 then
        i:=i+1;
        index[i]:=g[1][2];
fi;
if g[1][1]=-1 then
        i:=i+1;
        index[i]:=-g[1][2];
fi;
fi;
return index;
end;
############
SimplifyGeneratorsRepresentation:=function(index)
local
        i,j,k,p,
        temp,
        S,T,Y;
p:=0;
S:=[[0,-1],[1,0]];
T:=[[1,1],[0,1]];
Y:=[[0,-1],[1,1]];
temp:=[[1,0],[0,1]];
i:=Length(index);
while i>0 do
        if i mod 2 =1 then
                temp:=T^(index[i])*temp;
        else temp:=S*temp;
        fi;
        for j in [1..6] do
                if temp=Y^j then
                        if index[Length(index)]=0 then p:=1;fi;
                        for k in [i..Length(index)] do
                                Remove(index);
                        od;
                        if p=1 then Add(index,0);fi;
                fi;
        od;
        i:=i-1;
od;
if index=[0] then return [];fi;
return index;
end;
#############
EdgeFinder:=function(index)
local
        edge,S,T,
        id,m,
        sign;
if index=[] then return []; fi;
S:=[[0,-1],[1,0]];
T:=[[1,1],[0,1]];
id:=[[1,0],[0,1]];
if index=[1] then return id;
fi;
if index=[0] then return id;
fi;
if index=[-1] then return T^-1;
fi;
if index[1]=0 then Remove(index,1);fi;
if (Length(index) mod 2)=1 then
        sign:=index[1];
        index[1]:=index[1]-SignInt(index[1]);
        edge:=T^SignInt(sign)*EdgeFinder(index);
else
        Remove(index,1);
        edge:=S*EdgeFinder(index);
fi;
return edge;
end;
###################
Edge:=function(g)
g:=CanonicalRightCosetElement(StabilizerGroups[1][1],g^-1)^-1;
return EdgeFinder(g);
end;
##################
RefineEdge:=function(g)
local S;
S:=[[0,-1],[1,0]];
if g=[] then return [];fi;
if g[1][1]<=0 then
        g:=g*S^2;
        if g[1][1]<g[1][2] then
                g:=g*S;
        fi;
fi;
return g;
end;
#################
AddList:=function(g,h)
Add(g,h);
return g;
end;
#################
Chomotopy:=function(g)
local
        index,n,i,
        Y,S,T,d,h,
        path,edge,K;
if g=[] then return [];fi;
S:=[[0,-1],[1,0]];
T:=[[1,1],[0,1]];
Y:=[[0,-1],[1,1]];
path:=[];
for i in [1..6] do
        if g=Y^i then return [];
        fi;
od;
K:=Group(Y);
g:=CanonicalRightCosetElement(K,g^-1)^-1;
d:=SimplifyGeneratorsRepresentation(GeneratorsRepresentation(g));
#edge:=RefineEdge(EdgeFinder(d));
edge:=EdgeFinder(d);
edge:=CanonicalRightCosetElement(St1,edge^-1)^-1;
h:=[Sign(g,edge),edge];
return AddList(Chomotopy(edge*S*edge^-1*g),h);
end;
#################
pos:=function(g)
if Position(Elts,g)=fail then
        Add(Elts,g);
fi;
return Position(Elts,g);
end;
##################
Path2Gindex:=function(path)
local g,i;
g:=[];
for i in [1..Length(path)] do
        Add(g,[path[i][1],pos(path[i][2])]);
        #Add(g,[1,pos(CanonicalRightCosetElement(St1,path[i]^-1)^-1)]);
        #Add(g,[1,pos(path[i])]);
od;
#if Sign(g)=1 then return g;
#else
#       return NegateWord(g);
#fi;
return g;
end;

##################
#Sign:=function(index)
#local g;
#if index=[] then return 1;
#else
#g:=Elts[index[1][2]];
#if g=CanonicalRightCosetElement(St1,[[1,0],[0,1]]^-1)^-1 then return -1;fi;
#if g=CanonicalRightCosetElement(St1,[[1,1],[0,1]])^-1 then return 1;fi;
#if g=CanonicalRightCosetElement(St1,[[0,-1],[1,1]]^-1)^-1 then return -1;fi;
#fi;
#end;
##################
Sign:=function(g,h)
if g^-1*h in St0 then return 1;
else return -1;fi;
end;
##################
FinalHomotopy:=function(n,p)
local
        k,g;

if n<>0 then return [];
else
        if IsList(p[1]) then return Path2Gindex(Chomotopy(p));
        else
                k:=p[1];
                if AbsInt(k)<>1 then return "Number of Generators is 1";
                else
                        g:=p[2];
                        if Elts[g]=[] then return [];fi;
                        if k>0 then return Path2Gindex(Chomotopy(Elts[g]));
                        else
                            return NegateWord(Path2Gindex(Chomotopy(Elts[g])));
                        fi;
                fi;
        fi;
fi;
end;
####################END: SL2Z-homotopy###################################
Homotopy:=function(n,g)
if name="SL2Z" then return FinalHomotopy(n,g);
else
        return fail;
fi;
end;
#########################################################################



G:=SL(2,Integers); 

return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=FinalHomotopy,
            elts:=Elts,
            group:=G,
            stabilizer:=Stabilizer,
            action:=Action,
	    edge:=Edge,
	    gens:=GeneratorsRepresentation,
            properties:=
            [["length",Maximum(1000,lnth)],
             ["characteristic",0],
             ["type","resolution"],
             ["reduced",true]]  ));

end);
################################################
################################################


