
##########################################################
##########################################################
InstallGlobalFunction(UniversalCover,
function(X)
local U, G, dim, Elts, FreeElts, Boundary,PseudoBoundary, gamma, epi, n, i;

if not IsHapRegularCWComplex(X) then
Print("The function applies only to regular CW complexes.\n");
return fail;
fi;
if Length(PiZero(X)[1])>1 then
Print("The function applies only to path connected regular CW complexes.\n");
return fail;
fi;

OrientRegularCWComplex(X);
dim:=EvaluateProperty(X,"dimension");
G:=FundamentalGroupOfRegularCWComplex(X,"nosimplify");
epi:=EpimorphismFromFreeGroup(G);
gamma:=G!.edgeToWord;
Elts:=[One(G)];
FreeElts:=[PreImagesRepresentative(epi,One(G))];
PseudoBoundary:=List([1..dim],i->[]);


##################################
Boundary:=function(n,k)
local vts,g, fg, ng, B, bb, bbb, xx, ii, j,kk, pos, bool, indx, bnd,ornt;

if n=0 then return []; fi;
if IsBound(PseudoBoundary[n][k]) then
return 1*PseudoBoundary[n][k]; fi;

if n=1 then
#########################
vts:=X!.boundaries[2][k]; 
ornt:=X!.orientation[2][k];
g:=gamma(k); 
fg:=PreImagesRepresentative(epi,g);
ng:=Position(FreeElts,fg);
if ng=fail then
   Add(Elts,g);
   Add(FreeElts,fg);
   ng:=Length(Elts);
fi;
PseudoBoundary[n][k]:= [[ornt[1]*vts[2],1], [ornt[2]*vts[3],ng]];

return 1*PseudoBoundary[n][k];
#########################
fi;

#########################
B:=1*X!.boundaries[n+1][k];
B:=1*B{[2..Length(B)]};
ornt:=X!.orientation[n+1][k];
bnd:=List([1..Length(B)],j->[ornt[j]*B[j],1]);
indx:=1*B{[2..Length(B)]};

while Length(indx)>0 do

bool:=false;
for ii in indx do
  bb:=1*Boundary(n-1,ii);
  bbb:=1*List(bb,x->AbsInt(x[1]));
  for j in Difference(B,indx) do
  xx:=Boundary(n-1,j);
     for kk in [1..Length(xx)] do
     pos:= Position(bbb, AbsInt(xx[kk][1]));
     if not pos=fail then
       #g:=(Elts[xx[kk][2]]*Elts[bb[pos][2]]^-1); I THINK THIS WAS THE SLIP!
g:=(Elts[bnd[Position(B,j)][2]]*Elts[xx[kk][2]]*Elts[bb[pos][2]]^-1);
       fg:=PreImagesRepresentative(epi,g);
       ng:=Position(FreeElts,fg);
          if ng=fail then
              Add(Elts,g);
              Add(FreeElts,fg);
              ng:=Length(Elts);
          fi;
       bnd[Position(B,ii)][2]:=ng; 
       indx:=Filtered(indx,a->not a=ii); bool:=true;
       break;
     fi;
     if bool then break; fi;
     od;
  if bool then break; fi;
  od;
  if bool then break; fi;
od;

od;
PseudoBoundary[n][k]:= bnd;
#########################

return 1*PseudoBoundary[n][k];
end;
##################################


U:=         Objectify(HapEquivariantCWComplex,
            rec(
            dimension:=X!.nrCells,
            boundary:=Boundary,
            elts:=Elts,
            group:=G,
            stabilizer:=Group(One(G)),
            baseSpace:=X,
            properties:=
            [["dimension",dim], 
            ]  ));

for n in [1..dim] do
for i in [1..U!.dimension(n)] do
U!.boundary(n,i);
od;od;

return U;

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(EquivariantCWComplexToRegularCWComplex,
function(U,H)
local boundaries, bnd, gbnd, ind, trans,pair2int, dimU, n, k,g,W;

if not IsHapEquivariantCWComplex(U)  then
    Print("This function applies only to G-CW-complexes.\n");
    return fail;
fi;

if not IsSubgroup(U!.group,H) then
    Print("The provided group is not a subgroup of the fundamental group of the G-CW-complex.\n");
    return fail;
fi;

##Apply(U!.elts,x->x^-1); #So now we have a right action!
trans:=RightCosets(U!.group,H);
ind:=Length(trans);
if not ind < infinity then
    Print("The provided subgroup is not of finite index.\n");
    return fail;
fi;
dimU:=EvaluateProperty(U,"dimension");


###############################
pair2int:=function(e,gH);
return (e-1)*ind + gH;
end;
###############################

boundaries:=[];
boundaries[1]:=List([1..ind*U!.dimension(0)],i->[1,0]);
for n in [1..dimU] do
boundaries[n+1]:=[];
for k in [1..U!.dimension(n)] do
bnd:=U!.boundary(n,k);

for g in trans do
gbnd:=List(bnd,x->[x[1],
Position(trans, g*U!.elts[x[2]])]);

gbnd:=List(gbnd,x->pair2int(AbsInt(x[1]),x[2]));
gbnd:=Concatenation([Length(gbnd)],gbnd);
Add(boundaries[n+1],gbnd);
od;
od;
od;
boundaries[dimU+2]:=[];

W:=HAPRegularCWComplex(boundaries);
W!.index:=ind;
return W;

end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallGlobalFunction(EquivariantCWComplexToRegularCWMap,
function(U,H)
local YH,Y, map, ind;

YH:=EquivariantCWComplexToRegularCWComplex(U,H);
Y:=U!.baseSpace;
ind:=YH!.index;

#######################
map:=function(n,k)
local m, a;
m:=k mod ind;
a:=Int(k/ind);
if m=0 then return a;
else return a+1; fi;
end;
#######################

return Objectify(HapRegularCWMap,
       rec(
           source:=YH,
           target:=Y,
           mapping:=map));


end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallGlobalFunction(ChainComplexOfUniversalCover,
function(arg)
local X, Y, dim, bool, nrCriticalCells,C,Boundary, critical,
         BASIS,BIJ,DEFORM,DEFORMrec, f, bnd, sn, def, def1, def2,
         mult,inv, BOUNDARY, BNDrec;

X:=ContractedComplex(arg[1]);
dim:=Dimension(X);
X!.allcocriticalcells:=dim;
Y:=UniversalCover(X);
bool:=true;

if Length(arg)>1 then
bool:=arg[2];
fi;

####################################################################
if  not bool then 
return Objectify(HapEquivariantChainComplex,
            rec(
            dimension:=Y!.dimension,
            boundary:=Y!.boundary,
            elts:=Y!.elts,
            group:=Y!.group,
            properties:=
            [["dimension",dim], ["characteristic",0], ["length",dim]
            ]  ));
fi;
####################################################################

#########################################################
mult:=function(g,h)
local m, pos;
m:=Y!.elts[g]*Y!.elts[h];
#pos:=Position(Y!.elts,m);  #CAN IMPROVE THIS
pos:=fail;
if pos=fail then pos:=Length(Y!.elts)+1; Add(Y!.elts,m); fi;
return pos;
end;
#########################################################

#########################################################
inv:=function(g)
local m, pos;
m:=Y!.elts[g]^-1;
# pos:=Position(Y!.elts,m); #CAN IMPROVE THIS
pos:=fail;
if pos=fail then pos:=Length(Y!.elts)+1; Add(Y!.elts,m); fi;
return pos;
end;
#########################################################



critical:=CriticalCells(X);
C:=ChainComplexOfRegularCWComplexWithVectorField(X);
BASIS:=C!.basis;
BIJ:=C!.bij;
dim:=Maximum(List(CriticalCells(X),x->x[1]));

DEFORMrec:=List([1..dim+1],i->[]);;
########################################################
DEFORM:=function(n,c)
local k,kk,kkk,sgnk,sgnn,g,r,x;


kk:=c[1];
k:=AbsInt(kk);
sgnk:=SignInt(kk);
g:=c[2];

if [n,k] in critical then return [c]; fi;

if n>0 then
if IsBound(X!.vectorField[n][k]) then return [];fi;
fi;

if IsBound(DEFORMrec[n+1][k]) then
r:=List(DEFORMrec[n+1][k],a->[a[1],mult(g,a[2])]);
if sgnk=1 then return r;
else
return List(r,a->[-a[1],a[2]]);
fi;
fi;

f:=X!.inverseVectorField[n+1][k];
bnd:=Y!.boundary(n+1,f);
kkk:=Position(List(bnd,a->AbsInt(a[1])),k);
kkk:=bnd[kkk][2];
kkk:=inv(kkk);
#Apply(bnd,a->[a[1],mult(g, mult(kkk,a[2])     )]);
Apply(bnd,a->[a[1], mult(kkk,a[2]) ]);
sn:=X!.orientation[n+2][f];

def:=[]; def1:=[];def2:=[];
for x in [1..Length(bnd)] do
if not AbsInt(bnd[x][1])=k then
Add(def1,[bnd[x][1],bnd[x][2]]);
else
sgnn:=sn[x];
break;
fi;
od;
cnt:=x+1;

for x in [cnt..Length(bnd)] do
Add(def2,[bnd[x][1],bnd[x][2]]);
od;

if sgnn=1 then   
def:=Concatenation(def1,def2);
def:=List(def,a->[-a[1],a[2]]);
else
def:=Concatenation(def2,def1);
fi;

def:=List(def,x->DEFORM(n,x));
def:=Concatenation(def);

def:=Filtered(def,x->Length(x)>0);

def:=AlgebraicReduction(def);

DEFORMrec[n+1][k]:=def;

def:=List(def,a->[a[1],mult(g,a[2])]);

if sgnk=1 then return def; fi;  
def:=List(def,a->[-a[1],a[2]]);
return def;


end;
########################################################

# outputs the number of critical n-cells in C_n
####################################################################
        nrCriticalCells:=function(n)
                return Length(Filtered(X!.criticalCells, x->x[1]=n));
        end;
####################################################################
# boundary of critical cells -- reindexed
####################################################################
Boundary:=function(n,k)
local bnd, BND, x,y,L;
bnd:=Y!.boundary(n,BASIS[n+1][k]);

BND:=[];
for x in bnd do

L:=DEFORM(n-1,x);

for y in L do
Add(BND,[SignInt(y[1])*BIJ[n][AbsInt(y[1])],y[2]]);
od;od;

return BND;
end;
####################################################################

BNDrec:=List([1..dim],i->[]);
####################################################################
BOUNDARY:=function(n,k)
local kk,sk;

if n<1 or n>dim then return []; fi; 
kk:=AbsInt(k); sk:=SignInt(k);

if not IsBound(BNDrec[n][kk]) then
BNDrec[n][kk]:=Boundary(n,kk);
fi;

if sk=1 then return BNDrec[n][kk];
else return List(BNDrec[n][kk],a->[-a[1],a[2]]);
fi;
end;
####################################################################

return Objectify(HapEquivariantChainComplex,
            rec(
            dimension:=nrCriticalCells,
            boundary:=BOUNDARY,
            elts:=Y!.elts,
            group:=Y!.group,
            properties:=
            [["dimension",dim],["characteristic",0], ["length",dim]
            ]  ));



end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(TensorWithIntegersOverSubgroup,
function(C,H)
local R;

if not IsBound(C!.homotopy) then
C!.homotopy:=fail; fi;

R:=ResolutionSubgroup(C,H);
return TensorWithIntegers(R);
end);
##########################################################
##########################################################

