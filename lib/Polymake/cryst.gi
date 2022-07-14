
############################################################
############################################################
InstallGlobalFunction(CrystallographicComplex,
function(arg)
local G,K,pt,AG,SAG,F,V,FL,Boundaries,ind,inv,
Cells,VCells,GCells,PGCells,setV,
Dimn,BNDS, Bndy, REPS, NREPS,OREPS, Elts,Canonical,
STAB, STABREC, ACTION, Action, StandardWord, StabilizerBasis,
n,Y,k,x,y,i;

#In constructing K we'll use right actions, and at the end of the 
#construction we'll convert to left actions.
#The Stabilizer and Action and StandardWord functions need care -- 
#they will only ever be used with left actions.

G:=arg[1];
G:=StandardAffineCrystGroup(G);
if Length(arg)>1 then pt:=arg[2]; else 
pt:=List([1..Length(One(G))-1], i->1/Primes[i]);
fi;
####CREATE FUNDAMENTAL DOMAIN AND CORRESPONDING CW-COMPLE####
AG:=AffineCrystGroupOnRight(GeneratorsOfGroup(G));;
SAG:=StandardAffineCrystGroup(AG);;
F:=FundamentalDomainStandardSpaceGroup(pt,SAG);;
V:=Polymake(F,"VERTICES");;
setV:=Set(V);
FL:=PolymakeFaceLattice(F,true);
FL:=FL{[1..Length(FL)-1]};
ind:=List(FL,h->Minimum(Flat(h)));

for n in [1..Length(ind)] do
FL[n]:=FL[n]-ind[n]+1;
od;
ind:=List(FL,h->Length(h));

##########################################################
Boundaries:=[];
Boundaries[1]:=List([1..ind[1]],x->[1,0]);
for n in [2..Length(ind)] do
Boundaries[n]:=List([1..ind[n]],x->[]);
for x in [1..Length(FL[n-1])] do
for i in FL[n-1][x] do
Add(Boundaries[n][i],x);
od;
od;
Boundaries[n]:=List(Boundaries[n],b->Concatenation([Length(b)],SSortedList(b)));
od;
Boundaries[n+1]:=[Concatenation([ind[n]],[1..ind[n]])];
Boundaries[n+2]:=[];
##########################################################

Y:=RegularCWComplex(Boundaries);
###FUNDAMENTAL DOMAIN AND CW-COMPLEX CREATED##############

#return Y;

Cells:=[];  #THIS RECORDS INTEGER VERTICES FOR EACH CELL
Cells[1]:=List([1..Y!.nrCells(0)],i->[i]);
for k in [1..Dimension(Y)] do
Cells[k+1]:=List(Y!.boundaries[k+1], x->
SSortedList(  Concatenation( List(x{[2..Length(x)]}, i->Cells[k][i])  ))
);
od;

VCells:=[]; #THIS RECORDS VECTOR VERTICES FOR EACH CELL
for k in [1..Length(Cells)] do
VCells[k]:= List(Cells[k], S -> SortedList(List(S,i->V[i])));
od;

############################
inv:=function(X)
local st;
st:=Set(List(X,v->V[v]));
st:=OrbitPartInVertexSetsStandardSpaceGroup(SAG,st,setV);
return st;
end;
############################

GCells:=[]; #THIS RECORDS CELL G-ORBITS USING INTEGER CELL VERTICES
for k in [1..Length(Cells)] do
GCells[k]:=Classify(  Cells[k] , inv);
od;

PGCells:=[]; #THIS RECORDS CELL G-ORBITS USING POSITION IN Y
for k in [1..Length(GCells)] do
PGCells[k]:=[];
for x in GCells[k] do
y:=List(x, i->Position(Cells[k],i));
Add(PGCells[k],y);
od;
od;

OREPS:=[];  #THIS CONVERTS A CELL NUMBER IN Y TO ITS ORBIT REP NUMBER IN Y
for k in [1..Length(PGCells)] do
OREPS[k]:=[];
for x in PGCells[k] do
for i in x do
OREPS[k][i]:=x[1];
od;
od;
od;

NREPS:=[]; #THIS RECORDS THE POSITION IN Y OF AN ORBIT REP IN THE G-COMPLEX
for k in [1..Length(Cells)] do
REPS:=List(GCells[k],x->x[1]);
NREPS[k]:=List(REPS,x->Position(Cells[k],x));
od;

Elts:=[One(SAG)];

###############################################
Canonical:=function(n,k)
local kk, ll, g, pos;
#inputs a dimension n and cell number k in Y
#returns a pair [kk,g] with kk cell number of the orbit
#rep and g the number in Elts such that g.kk=k

kk:=OREPS[n+1][k];
ll:=Position(NREPS[n+1],kk);
g:=RepresentativeActionOnRightOnSets(SAG,VCells[n+1][kk], VCells[n+1][k]);
if g=fail then Print("Error!!!\n"); return fail;fi; #THIS WON'T HAPPEN
pos:=Position(Elts,g);
if pos=fail then Add(Elts,g); pos:=Length(Elts); fi;
return [ll,pos];
end;
###############################################


########################################
Dimn:=function(n);
if n<0 or n>Dimension(Y) then return 0; fi;
return Length(GCells[n+1]);
end;
########################################

BNDS:=[];
for n in [2..Length(NREPS)] do
BNDS[n-1]:=[];
for k in NREPS[n] do
x:=1*Y!.boundaries[n][k]{[2..Y!.boundaries[n][k][1]+1]};
y:=List(x, i->Canonical(n-2,i));
for i in [1..Length(Y!.orientation[n][k])] do
y[i]:=[Y!.orientation[n][k][i]*y[i][1],y[i][2]];
od;
Add(BNDS[n-1], y );
od;
od;

########################################
Bndy:=function(n,k);
if n<1 or n>Dimension(Y) then return []; fi;
return BNDS[n][k];
end;
########################################

########################################
STAB:=function(n,kk)
local st,k,  gens;
k:=NREPS[n+1][kk];
st:=VCells[n+1][k];
st:=StabilizerOnSetsStandardSpaceGroup(SAG,st);
st:=List(st,a->TransposedMat(a));
st:=Group(st);
return st;
end;
########################################

STABREC:=[];
for n in [1..Length(NREPS)] do
STABREC[n]:=[];
for k in [1..Length(NREPS[n])] do
STABREC[n][k]:=STAB(n-1,k);
od;
od;

########################################
STAB:=function(n,k);
return STABREC[n+1][k];
end;
########################################


#####################################################################
StandardWord:=function(k,bnd)
local w,x,c,pos;
w:=[];
for x in bnd do
c:=CanonicalRightCosetElement(STAB(k,AbsInt(x[1])), Elts[x[2]]^-1 )^-1;
pos:=Position(Elts,c);
if pos=fail then Add(Elts,c); pos:=Length(Elts); fi;
Add(w,[x[1]*ACTION(k,AbsInt(x[1]),x[2]),pos]);
#Add(w,[x[1],pos]);
od;
return AlgebraicReduction(w);
end;
#####################################################################

####################################################
Action:=function(uu,V)
local L, M, T, d, u;

u:=TransposedMat(uu);


d:=Length(V);

L:=1*u{[1..d]};
M:=List(L,r->r{[1..d]});
T:=1*u[d+1]{[1..d]}  ;

return V*M + T;
end;
####################################################


###############################################################
# This describes how the group G acts on the orientation.
ACTION:=function(n,k,h)
local bas, Gbas, mat,id,r,u,H,A,B,p,kk,cg,d;

if n=0 then return 1; fi;
if Elts[h]=One(SAG) then return 1; fi;
H:=STAB(n,AbsInt(k));


id:=CanonicalRightCosetElement(H,Identity(SAG));
r:=CanonicalRightCosetElement(H,Elts[h]^-1);
r:=id^-1*r;
u:=r*Elts[h];

if not u in H then Print("ERROR\n"); fi;

kk:=NREPS[n+1][AbsInt(k)];
A:=VCells[n+1][kk];
cg:=Sum(A)/(Length(A));
A:=List(A,v->v-cg);
bas:=SemiEchelonMat(A).vectors;
Gbas:=List(bas,V->Action(u,V+cg)-cg);
mat:=List(Gbas,b->SolutionMat(bas,b));
d:=Determinant(mat);
if not AbsInt(d)=1 then Print("ERROR\n"); fi;
return Determinant(mat);

end;
###############################################################



#So far we have constructed a complex of RIGHT modules.
#But we've already fixed STAB to left module convention.
Elts:=List(Elts,a->TransposedMat(a));
SAG:=GeneratorsOfGroup(SAG);
SAG:=List(SAG,a->TransposedMat(a));
SAG:=Group(SAG);


K:=    Objectify(HapNonFreeResolution,
           rec(
            dimension:=Dimn,
            boundary:=Bndy,
            homotopy:=fail,
            elts:=Elts,
            group:=SAG,
            stabilizer:=STAB,
            action:=ACTION,
            standardWord:=StandardWord,
            properties:=
             [["type","non-free resolution"],
              ["length",1000],
              ["characteristic", 0] ]));

RecalculateIncidenceNumbers(K);

return K;
end);
###############################################################
###############################################################

###############################################################
###############################################################
InstallGlobalFunction(ResolutionSpaceGroup,
function(arg)
local G, N, V, R, K, pos;

G:=arg[1];
N:=arg[2];
if Length(arg)>2 then V:=arg[3]; fi;

if Length(arg)=2 then
K:=CrystallographicComplex(G);
else
K:=CrystallographicComplex(G,V);
fi;

R:=FreeGResolution(K,N);
#pos:=PositionProperty(R!.properties,x->x[1]="characteristic");
#R!.properties[pos][2]:=2;
return R;
end);
###############################################################
###############################################################
