#(C) 2009 Graham Ellis

###############################################################
###############################################################
InstallGlobalFunction(PersistentHomologyOfQuotientGroupSeries_Int,
function(arg)
		#S is a chain of normal subgroups with S[1]=G.   
local
	S,deg,Resol,RESLST,HOMS,A,B,
	N,P,Q,RP,RQ,T,Persistence,hom,eqmap,map,i,j;

#####INPUT#################
S:=arg[1]; #A normal series of subgroups
deg:=arg[2]; #The degree in which homology will be calculated.
if Length(arg)>2 then Resol:=arg[3]; else 
Resol:=ResolutionFiniteGroup;
fi;
###########################

HOMS:=NormalSeriesToQuotientHomomorphisms(S);;
if Order(Range(HOMS[Length(HOMS)]))=1 then
HOMS:=HOMS{[1..Length(HOMS)-1]};
fi;
RESLST:=[];

RP:=Resol(Source(HOMS[1]),deg+1);
Add(RESLST,RP);

for i in [1..Length(HOMS)] do
RQ:=Resol(Range(HOMS[i]),deg+1);
Add(RESLST,RQ);
od;

Persistence:=List( [1..Length(HOMS)+1], i->List([1..Length(HOMS)+1],j->0));

for i in [1..Length(HOMS)] do
for j in [i..Length(HOMS)] do
if i=j then hom:=IdentityMapping(Source(HOMS[i]));
else
hom:=Product(HOMS{[i..j-1]});
fi;
eqmap:=EquivariantChainMap(RESLST[i],RESLST[j],hom);
map:=TensorWithIntegers(eqmap);
map:=Homology(map,deg);
#map:=Image(map);
if Size(Image(map))=1 then A:=[];
else A:=AbelianInvariants(Image(map));
fi;
if Size(Range(map)/Image(map))=1 then B:=[];
else B:=AbelianInvariants(Range(map)/Image(map));
fi;

Persistence[i][j]:=[A,B];
od;
od;

return Persistence;

end);
###############################################################
###############################################################

###############################################################
###############################################################
InstallGlobalFunction(PersistentCohomologyOfQuotientGroupSeries,
function(arg)
                #S is a chain of normal subgroups with S[1]=G.
local
        S,deg,Resol,RESLST,HOMS,A,B,
        N,P,Q,RP,RQ,T,prime,Persistence,hom,eqmap,map,LN,i,j;

#####INPUT#################
S:=arg[1]; #A normal series of subgroups
deg:=arg[2]; #The degree in which homology will be calculated.
if Length(arg)>2 then 
prime:=arg[3];
else prime:=PrimePGroup(S[1]); fi;
if Length(arg)>3 then
Resol:=arg[4]; else
  if not prime=0 then
  Resol:=ResolutionPrimePowerGroup;
  else
  Resol:=ResolutionFiniteGroup;
  fi;
fi;
###########################

HOMS:=NormalSeriesToQuotientHomomorphisms(S);;
if Order(Range(HOMS[Length(HOMS)]))=1 then
HOMS:=HOMS{[1..Length(HOMS)-1]};
fi;
RESLST:=[];

RP:=Resol(Source(HOMS[1]),deg+1);
Add(RESLST,RP);

for i in [1..Length(HOMS)] do
RQ:=Resol(Range(HOMS[i]),deg+1);
Add(RESLST,RQ);
od;

LN:=Length(HOMS)+1;
Persistence:=List( [1..LN], i->List([1..LN],j->0));

for i in [1..LN] do
for j in [i..LN] do
if i=j  then 
        if i=LN then
	hom:=IdentityMapping(Image(HOMS[LN-1]));
        else
	hom:=IdentityMapping(Source(HOMS[i]));
	fi;
else
 	hom:=Product(HOMS{[i..j-1]});

fi;
eqmap:=EquivariantChainMap(RESLST[i],RESLST[j],hom);
if not prime=0 then
map:=HomToIntegersModP(eqmap,prime);
map:=Cohomology(map,deg);
Persistence[i][j]:=Log(Order(Image(map)),prime);

else

map:=HomToIntegers(eqmap);
map:=Cohomology(map,deg);
A:=Image(map);
B:=Range(map)/A;
  if Order(A)>1 then A:=AbelianInvariants(A); else A:=[]; fi;
  if Order(B)>1 then B:=AbelianInvariants(B); else B:=[]; fi;
Persistence[i][LN-j]:=[A,B];

fi;
od;
od;

Persistence:=Persistence{[2..Length(Persistence)]};
Persistence:=TransposedMat(Persistence);
Persistence:=Persistence{[2..Length(Persistence)]};
Persistence:=TransposedMat(Persistence);

return Persistence;

end);
###############################################################
###############################################################


###############################################################
###############################################################
InstallGlobalFunction(NormalSeriesToQuotientHomomorphisms,
function(SS)
local
        S,
        HOMS,N,P,G,iso,T,hom,hom1,i;

if Size(SS[1])<=Size(SS[Length(SS)])
then S:=SS;
else
S:=Reversed(SS);
fi;

T:=[];
for i in [1..Length(S)] do
T[i]:=NaturalHomomorphismByNormalSubgroup(S[Length(S)],S[i]);
od;

HOMS:=[T[1]];
for i in [1..Length(S)-1] do
P:=Image(T[i]);
N:=List(GeneratorsOfGroup(S[i+1]),x->Image(T[i],x));
Add(N,Identity(P));
N:=Subgroup(P,N);
hom:=NaturalHomomorphismByNormalSubgroup(P,N);
iso:=IsomorphismGroups(Image(hom),Image(T[i+1]));
hom1:=GroupHomomorphismByImages(P,Image(T[i+1]),
GeneratorsOfGroup(P),
List(GeneratorsOfGroup(P), x-> Image(iso,Image(hom,x))));
Add(HOMS,hom1);
od;

return HOMS;

end);
###############################################################
###############################################################


###############################################################
###############################################################
InstallGlobalFunction(PersistentHomologyOfQuotientGroupSeries,
function(arg)
                #S is a chain of normal subgroups with S[1]=G.
local
        S,deg,prime,Resol,
        HOMS,RP,RQ,T,Persistence,MAPS,eqmap,map,i;

#####INPUT#################
S:=arg[1]; #A normal series of subgroups
deg:=arg[2]; #The degree in which homology will be calculated.
if Length(arg)>2 then
prime:=arg[3]; else prime:=PrimePGroup(S[1]);
fi;
if Length(arg)>3 then Resol:=arg[4]; else
Resol:=ResolutionPrimePowerGroup;
fi;
###########################

if prime=0 then
   if Length(arg)>3 then
   return PersistentHomologyOfQuotientGroupSeries_Int(S,deg,Resol);
   else
   return PersistentHomologyOfQuotientGroupSeries_Int(S,deg);
   fi;
fi;

HOMS:=NormalSeriesToQuotientHomomorphisms(S){[1..Length(S)-1]};
MAPS:=[];

if Length(HOMS)=1 and Resol=ResolutionPrimePowerGroup then

  RP:=Resol(Source(HOMS[1]),deg);
  MAPS:=[
  IdentityMapping(GF(prime)^RP!.dimension(deg))
  ];

else

  RP:=Resol(Source(HOMS[1]),deg+1);
  MAPS:=[
  IdentityMapping(HomologyVectorSpace(TensorWithIntegersModP(RP,prime),deg))
  ];

fi;

for i in [2..Length(HOMS)] do
RQ:=Resol(Image(HOMS[i]),deg+1);
eqmap:=EquivariantChainMap(RP,RQ,HOMS[i]);
map:=TensorWithIntegersModP(eqmap,prime);
map:=HomologyVectorSpace(map,deg);
Add(MAPS,map);
RP:=RQ;
od;

MAPS:=LinearHomomorphismsPersistenceMat(MAPS);
MAPS:=MAPS{[2..Length(MAPS)]};
MAPS:=TransposedMat(MAPS);
MAPS:=MAPS{[2..Length(MAPS)]};
MAPS:=TransposedMat(MAPS);

return MAPS;
end);
###############################################################
###############################################################


###############################################################
###############################################################
InstallGlobalFunction(LinearHomomorphismsPersistenceMat,
function(HOMS)
local
	mat,x,y,L;

L:=Length(HOMS)+1;
mat:=IdentityMat(L)*0;

for x in [1..L] do
for y in [x..L] do
if y=x then
  if x=L then
  mat[x][y]:=IdentityMapping(Range(HOMS[x-1]));
  else
  mat[x][y]:=IdentityMapping(Source(HOMS[x]));
  fi;
else
mat[x][y]:= Product(List([x..y-1],i->HOMS[i]));
fi;
mat[x][y]:=Dimension(Image(mat[x][y]));
od;
od;

return mat;

end);
###############################################################
###############################################################


###############################################################
###############################################################
InstallGlobalFunction(PersistentHomologyOfSubGroupSeries,
function(arg)
local
	S,deg,prime,hom,HOMS,MAPS,RP,RQ,Resol,eqmap,map,i;

#####INPUT#################
S:=arg[1];
if Order(S[1])>Order(S[Length(S)]) then
S:=Reversed(S); 
fi; #A normal series of subgroups
if Order(S[1])=1 then S:=S{[2..Length(S)]};fi;
deg:=arg[2]; #The degree in which homology will be calculated.
if Length(arg)>2 then
prime:=arg[3]; else prime:=PrimePGroup(S[1]);
fi;
if Length(arg)>3 then Resol:=arg[4]; 
else Resol:=ResolutionPrimePowerGroup;
fi;
###########################

Add(S,S[Length(S)]);
HOMS:=[];

for i in [1..Length(S)-1] do
hom:=GroupHomomorphismByFunction(S[i],S[i+1],x->x);;
Add(HOMS,hom);
od;


MAPS:=[];

RP:=Resol(Source(HOMS[1]),deg+1);
for i in [1..Length(HOMS)] do
RQ:=Resol(Image(HOMS[i]),deg+1);
eqmap:=EquivariantChainMap(RP,RQ,HOMS[i]);
map:=TensorWithIntegersModP(eqmap,prime);
map:=HomologyVectorSpace(map,deg);
Add(MAPS,map);
RP:=RQ;
od;


MAPS:=LinearHomomorphismsPersistenceMat(MAPS);
MAPS:=MAPS{[2..Length(MAPS)]};
MAPS:=TransposedMat(MAPS);
MAPS:=MAPS{[2..Length(MAPS)]};
MAPS:=TransposedMat(MAPS);

return MAPS;

end);
###############################################################
###############################################################


###############################################################
###############################################################
InstallGlobalFunction(PersistentHomologyOfPureCubicalComplex,
function(arg)
local
	M,LM,deg,prime,toggle,triv,ChnCmps,vsmaps,A,S,TS,TM,L,bool,x,y,i,F;

##Input###################
deg:=arg[2];
prime:=arg[3];
toggle:=IsList(arg[1]);
##########################

##IF STARTS
if toggle then 
LM:=arg[1];
M:=LM[1];

  if deg=0 then return
  PersistentHomologyOfPureCubicalComplex_Alt(LM,deg,prime);
  fi;

else

M:=arg[1];

  if deg=0 then return
  PersistentHomologyOfPureCubicalComplex_Alt(M,deg,prime);
  fi;

fi;
##IF ENDS

F:=HomotopyEquivalentMinimalPureCubicalSubcomplex;
triv:=List([1..Dimension(M)+1],i->0); triv[1]:=1;

S:=ContractibleSubcomplexOfPureCubicalComplex(M);
if not toggle then LM:=[M]; fi;
L:=[[M,S]];

bool:=true;
i:=1;
while bool do
i:=i+1;
if toggle then
  TM:=LM[i];
else 
  TM:=ThickenedPureCubicalComplex(L[Length(L)][1]);
fi;
TS:=HomotopyEquivalentMaximalPureCubicalSubcomplex(TM,L[Length(L)][2]);
if toggle then 
  bool:=i<Length(LM);
  else
  bool:=not Bettinumbers(ContractedComplex(TM))=triv;
fi;
Add(L,[TM,TS]); 
od;

vsmaps:=[];

L[1][1]:=F(L[1][1],L[1][2]);
for x in [1..Length(L)-1] do
L[x+1][1]:=F(L[x+1][1],PureCubicalComplexUnion(L[x][1],L[x+1][2]));
od;


L[1]:=ChainComplexOfPair(L[1][1],L[1][2]);
#################
for x in [1..Length(L)-1] do
L[x+1]:=ChainComplexOfPair(L[x+1][1],L[x+1][2]);
vsmaps[x]:= 
TensorWithIntegersModP( ChainMapOfCubicalPairs(L[x],L[x+1]) , prime);
vsmaps[x]:=HomologyVectorSpace(vsmaps[x],deg);
Unbind(L[x]);
od;
#################

return LinearHomomorphismsPersistenceMat(vsmaps);

end);
###############################################################
###############################################################

###############################################################
###############################################################
InstallGlobalFunction(PersistentHomologyOfPureCubicalComplex_Alt,
function(arg)
local
        M,Bettis,triv,Persist,deg,prime,TM,L,SL,bool,CollapseMat,
        Collapses_dim2,Fun,x,y;

##Input###################
deg:=arg[2];
prime:=arg[3];
if IsList(arg[1]) then L:=arg[1];
M:=L[1];
else 
M:=arg[1];
L:=[M];
fi;
deg:=Minimum(deg,Dimension(M));
##########################

if not deg in [0,Dimension(M)-1,Dimension(M)] then 
return fail; 
fi;

if IsList(arg[1]) then
	Bettis:=List(L,x->[Bettinumbers(x,0),Bettinumbers(x,deg)]);
else
	triv:=List([1..Dimension(M)+1],i->0);
	triv[1]:=1;
	Bettis:=[Bettinumbers(ContractedComplex(M))];

	bool:=true;
	while bool do
	TM:=ThickenedPureCubicalComplex(L[Length(L)]);
	Add(Bettis,Bettinumbers(ContractedComplex(TM)));
	bool:= not Bettis[Length(Bettis)]=triv;;
	Add(L,TM); 
	od;
	Bettis:=List(Bettis,b->[b[1],b[deg+1]]);
fi;


Persist:=List([1..Length(L)],x->List([1..Length(L)],y->0));;

###############################
if deg=0 then
for x in [1..Length(L)] do
for y in [x..Length(L)] do
Persist[x][y]:=Bettis[y][1];
od;
od;
return Persist;
fi;
###############################

###############################
if deg=Dimension(M) then
return Persist;
fi;
###############################

###############################
Collapses_dim2:=function(X,Y)
##Number of cycles killed going from position X to position Y
local
	Acomp,A,Bcomp,TP, rows, cols, cnt, bettizero, Pathcomps, k, i, j;

Acomp:=PureCubicalComplex(FrameArray(L[Y]!.binaryArray));
A:=Acomp!.binaryArray;
Bcomp:=PureCubicalComplex(FrameArray(L[X]!.binaryArray));
rows:=Length(A);
cols:=Length(A[1]);
cnt:=0;
bettizero:=Bettis[X][1];

for k in [1..bettizero] do
P:=PathComponentOfPureCubicalComplex(Bcomp,k);

TP:=HomotopyEquivalentMaximalPureCubicalSubcomplex(Acomp,P);
TP:=TP!.binaryArray;

for i in [2..rows-1] do
for j in [2..cols-1] do
if A[i][j]=1 and TP[i][j]=0 then
if 4=Sum([TP[i-1][j]+TP[i+1][j]+TP[i][j-1]+TP[i][j+1]]) then cnt:=cnt+1;fi;;
fi;
od;
od;

od;

return cnt;
end;
###############################

#SL:=[ContractedComplex(L[1])];
#for x in [2..Length(L)] do
#Add(SL,HomotopyEquivalentMinimalPureCubicalSubcomplex(L[x],SL[x-1]));
#SL[x-1]:=HomotopyEquivalentMaximalPureCubicalSubcomplex(SL[x],SL[x-1]);
#od;
#L:=SL;


###############################
if deg=Dimension(M)-1 then

###############
Fun:=function(x,y);
#Ths function implements some partial recurrence in CollapseMat.
if CollapseMat[x+1][y]=0 then return CollapseMat[x][y-1]; fi;
if CollapseMat[x]=Bettis[x][2] then return CollapseMat[x][y-1]; fi;
if not y=Length(Bettis) then
   if Bettis[y+1]=0 then return Bettis[x][2]-Sum(CollapseMat[x]); fi;
fi;
return Collapses_dim2(x,y);
end;
###############

CollapseMat:=[];
for x in Reversed([1..Length(L)]) do
CollapseMat[x]:=List([1..Length(L)],a->0);
#CollapseMat[x]:=[];
if not x=Length(L) then 
CollapseMat[x][x+1]:=Collapses_dim2(x,x+1);
fi;
	for y in [x+2..Length(L)] do
	CollapseMat[x][y]:=Fun(x,y);;
	od;
od;

for x in [1..Length(L)] do
for y in [x..Length(L)] do
Persist[x][y]:=StructuralCopy(Bettis[x][2]);
if not x=y then
Persist[x][y]:=Persist[x][y]-CollapseMat[x][y];
fi;
od;
od;
return Persist;
fi;
###############################
end);
###############################################################
###############################################################


#########################################################
#########################################################
InstallGlobalFunction(BarCode,
function(P)
local PT,B,newrow,cols,rows,i,j,k,m;

PT:=MutableCopyMat(TransposedMat(P));
Add(PT,PT[1]*0);
PT:=TransposedMat(PT);
PT:=Concatenation([PT[1]*0],PT);
for i in Reversed([2..Length(PT)]) do
PT[i]:=PT[i]-PT[i-1];
od;

B:=[];
rows:=Length(P);
cols:=Length(P[1]);

for i in [1..rows] do
for j in Reversed([i..cols]) do
for k in [1..PT[i+1][j]-PT[i+1][j+1]] do
newrow:=List([1..cols],a->0);
for m in [i..j] do
newrow[m]:=1;
od;
Add(B,newrow);
od;

od;
od;

return B;
end);
###########################################################
###########################################################

###########################################################
###########################################################
InstallGlobalFunction(BarCodeDisplay,
function(arg)
local P,B,i,j,Disp,tmpDir,barcodedot,barcodegif;

###############
P:=arg[1];
if Length(arg)=2 then
Disp:=arg[2];
else
Disp:=DISPLAY_PATH;
fi;
###############

B:=BarCode(P);

tmpDir:=DirectoryTemporary();
barcodedot:=Filename(tmpDir,"tmpIn.log");
barcodegif:=Filename(tmpDir,"basic.gif");


PrintTo(barcodedot,"digraph finite_state_machine {\n\n");
AppendTo(barcodedot, "rankdir=LR;\n\n");
AppendTo(barcodedot,
"node [style=filled,shape=point]\n\n");

for i in [1..Length(B)] do
for j in [1..Length(B[1])] do
if B[i][j]=1 then
AppendTo(barcodedot,"node [color=black,fontcolor=black];\n",i,j,"\n");
else
AppendTo(barcodedot,"node [color=white,fontcolor=white];\n",i,j,"\n");
fi;
od;
od;

for i in [1..Length(B)] do
for j in [1..Length(B[1])-1] do
if B[i][j]=1 and B[i][j+1]=1 then
AppendTo(barcodedot,i,j,"->",i,j+1," [color=black,arrowhead=none];\n");
else
AppendTo(barcodedot,i,j,"->",i,j+1," [color=white,arrowhead=none];\n");
fi;
od;
od;

AppendTo(barcodedot,"}");

Exec(Concatenation("dot -Tgif ",barcodedot ,">", barcodegif));
Exec(Concatenation(Disp," ",barcodegif));
Sleep(2);
Exec(Concatenation("rm -r ",barcodegif{[1..Length(barcodegif)-9]}));

end);
#################################################################
#################################################################

