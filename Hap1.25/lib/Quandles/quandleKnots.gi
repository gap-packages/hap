#####################################################################
#####################################################################
InstallMethod( ViewObj,
 "for HapQuandlePresentation",
 [IsHapQuandlePresentation],
 function(R)
 Print("Quandle presentation of ", Length(R!.generators), " generators and ",Length(R!.relators), " relators.\n");  
end);

#####################################################################
#####################################################################
InstallMethod( PrintObj,
 "for HapQuandlePresentation",
 [IsHapQuandlePresentation],
 function(R)
 Print("Quandle presentation of ", Length(R!.generators), " generators and ",Length(R!.relators), " relators.\n");
end);


#####################################################################
#####################################################################
InstallGlobalFunction(PresentationKnotQuandle,
function(gaussCode)
local generators,relators,genAndRel,i,j,beg,arcs,arc1,arc2,n,a,b,c,sign;

#A relation w=v in a presentation will be stored as a list [v,w].
#If v is an integer then it represents generators[v];
#if v=[x,y] is a list then it represents the word x*y.

# List the generators
generators:=[1..Length(gaussCode[2])];

# Fill in arcs for 1..lengthGene-1
beg:=1;
while gaussCode[1][1][beg]>0 do beg:=beg+1; od; i:=beg;

arcs:=[];
for j in [1..Length(generators)-1] do

arc1:=[gaussCode[1][1][i]];
while gaussCode[1][1][i+1]>0 do i:=i+1; Add(arc1,gaussCode[1][1][i]); od;
i:=i+1;
Add(arc1,gaussCode[1][1][i]);
Add(arcs,[j,arc1]);

od;

# Fill in arcs for lengthGene
arc1:=[];
for j in [i..Length(gaussCode[1][1])] do Add(arc1,gaussCode[1][1][j]); od;
for j in [1..beg] do Add(arc1,gaussCode[1][1][j]); od;
Add(arcs,[Length(generators),arc1]);


# Fill in the relators
relators:=[];

for arc1 in arcs do

c:=arc1[1];
n:=-arc1[2][1];
for arc2 in arcs do
if n in arc2[2] then b:=arc2[1]; fi;
if (-n in arc2[2]) and (not arc1=arc2) then a:=arc2[1]; fi;
od;

sign:=gaussCode[2][n];
if sign<0 then Add(relators, [[a,b],c]);
else Add(relators,[[c,b],a]); fi;

od;

return Objectify(HapQuandlePresentation, rec(generators:=generators,relators:=relators));

end);

#####################################################################
#####################################################################
InstallGlobalFunction(PD2GC,
function(PD)
local GC,n,currentCross,cross;

GC:=[[[]],[]];
Add(GC[1][1],-1);
n:=PD[1][3];
currentCross:=PD[1];
while not n=PD[1][1] do
for cross in PD do
if n in cross and not cross=currentCross then
if Position(cross,n)=1 then Add(GC[1][1],-Position(PD,cross)); n:=PD[Position(PD,cross)][3]; currentCross:=cross; break; fi;
if Position(cross,n)=2 then Add(GC[1][1],Position(PD,cross)); n:=PD[Position(PD,cross)][4]; currentCross:=cross; break; fi;
if Position(cross,n)=4 then
Add(GC[1][1],Position(PD,cross)); n:=PD[Position(PD,cross)][2];
GC[2][Position(PD,cross)]:=1;
currentCross:=cross; break;
fi;
fi;
od;
od;

for n in [1..Length(PD)] do
if not IsBound(GC[2][n]) then GC[2][n]:=-1; fi;
od;
return GC;
end);

#####################################################################
#####################################################################
InstallGlobalFunction(PlanarDiagramKnot,
function(n,k);
if not IsBound(Cedric_PlanarDiagram[n][k]) then return fail; fi;
return Cedric_PlanarDiagram[n][k];
end);

#####################################################################
#####################################################################
InstallGlobalFunction(GaussCodeKnot,
function(n,k);
if not IsBound(Cedric_PlanarDiagram[n][k]) then return fail; fi;
return PD2GC(Cedric_PlanarDiagram[n][k]);
end);

#####################################################################
#####################################################################
InstallGlobalFunction(PresentationKnotQuandleKnot,
function(n,k);
if not IsBound(Cedric_PlanarDiagram[n][k]) then return fail; fi;
return PresentationKnotQuandle(PD2GC(Cedric_PlanarDiagram[n][k]));
end);

#####################################################################
#####################################################################
InstallGlobalFunction(Cedric_IsHomomorphism,
function(mapping,relations,multT)
local x,y,z,re;

for re in relations do
x:=re[1][1];
y:=re[1][2];
z:=re[2];
if not mapping[z]=multT[mapping[x]][mapping[y]] then return false; fi;
od;
return true;
end);

#####################################################################
#####################################################################
InstallMethod(NumberOfHomomorphisms,
"for a finite quandle presentation  and a finite quandle",
[IsHapQuandlePresentation,IsMagma],
function(genRelQ,finiteQ)
local multTab,orderFiniteQ,nbGen,nbHom,phi,i,j,bool,q,T;

if not IsQuandle(finiteQ) then TryNextMethod(); fi;

if IsConnected(finiteQ) and (not Size(finiteQ)=1) then T:=NumberOfHomomorphisms_connected(genRelQ,finiteQ); return Sum(T); fi;

multTab:=MultiplicationTable(finiteQ);
orderFiniteQ:=Size(finiteQ);
nbGen:=Length(genRelQ!.generators);
nbHom:=0;

phi:=ListWithIdenticalEntries(nbGen,1);

while not phi=ListWithIdenticalEntries(nbGen,orderFiniteQ) do

if Cedric_IsHomomorphism(phi,genRelQ!.relators,multTab) then nbHom:=nbHom+1; fi;

if phi[nbGen]<orderFiniteQ then phi[nbGen]:=phi[nbGen]+1;
else
i:=nbGen-1;
while phi[i]=orderFiniteQ do i:=i-1; od;
phi[i]:=phi[i]+1;
for j in [i+1..nbGen] do phi[j]:=1; od;
fi;

od;

if Cedric_IsHomomorphism(ListWithIdenticalEntries(nbGen,orderFiniteQ),genRelQ!.relators,multTab) then nbHom:=nbHom+1; fi;

return nbHom;
end);


#####################################################################
#####################################################################
InstallMethod(PartitionedNumberOfHomomorphisms,"for a KnotQuandle and a Connected Quandle",[IsRecord,IsMagma],
function(genRelQ,finiteQ)

if (not IsQuandle(finiteQ)) or (not IsConnected(finiteQ)) then TryNextMethod(); fi;

return NumberOfHomomorphisms_connected(genRelQ,finiteQ);
end);


#####################################################################
#####################################################################
InstallGlobalFunction(KnotInvariantCedric,
function(genRelQ,n,m)
local N,i,CQ,Q;

N:=[];

if n=1 then Add(N,[NumberOfHomomorphisms(genRelQ,ConnectedQuandle(1,1))]); n:=n+1; fi;

for i in [n..m] do
CQ:=ConnectedQuandles(i);
for Q in CQ do
Add(N,PartitionedNumberOfHomomorphisms(genRelQ,Q));
od; od;
return N;
end);
