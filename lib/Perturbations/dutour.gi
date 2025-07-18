HAP_GCOMPLEX_LIST:=0;
HAP_GCOMPLEX_SETUP:=[true];

################################################
################################################
InstallGlobalFunction(ContractibleGcomplex,
function(groupname)
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
        #InfGrps,
        x, lstMathieu, lstAlexanderSebastianSL, lstSebastianGL, lstAlexanderSebastianAlt,
        n,k,s,BI,SGN,tmp, LstEl , bool, name,dd;

n:=Position(groupname,'-');##This won't work for all file names
k:=Position(groupname,')');##So I'll re-introduce the number elsewhere
if IsInt(n) and IsInt(k) then dd:=groupname{[n+1..k-1]}; #dd:=-Int(dd);
else dd:=fail; fi;
groupname:=Filtered(groupname,x->not(x='(' or x=')' or x=',' or x='[' or x=']'));
name:=groupname;

if name="SL2Z" then return ContractibleSL2ZComplex(); fi;

groupname:=Concatenation("lib/Perturbations/Gcomplexes/",groupname);
bool:=ReadPackage("HAP",groupname);
#Print("\n\n",name,"\n\n");
#InfGrps:=["SL2O-5"];

if HAP_GCOMPLEX_SETUP[1] then 
TransMat:=function(x); return x^-1; end;
else
TransMat:=function(x); return x; end;
fi;


if bool = false then 
#if false then
Print("G-complexes are implemented for the following groups only:  \n\n");

lstMathieu:=" SL(2,Z) , SL(3,Z) , PGL(3,Z[i]) , PGL(3,Eisenstein_Integers) ,  PSL(4,Z) , Sp(4,Z) , \n\n PSL(4,Z)_b , PSL(4,Z)_c , PSL(4,Z)_d \n"; 

lstAlexanderSebastianAlt:=" SL(2,O-d)_a , d=2, 7, 11, 19  \n";

lstSebastianGL:=" GL(2,O-d) , d=1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 43  \n";

lstAlexanderSebastianSL:=" SL(2,O-d) , d=2, 3, 5, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 43, 67, 163 \n"; 
Print(lstMathieu,
"\n\n",
lstSebastianGL, 
"\n\n",
lstAlexanderSebastianSL,
"\n\n",
lstAlexanderSebastianAlt,
"\n\nwhere subscripts _a, _b,  etc denote alternative complexes for a given group and  O-d denotes the ring of integers of the number field Q(\sqrt(-d)).\n\n");
return fail;
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

Elts:=[Identity(C[1][1].TheMatrixStab)];   ##Modified 8 Feb 2012
StabilizerGroups:=[];
RotSubGroups:=[];
boundaryList:=[];


#######
for n in [1..lnth+1] do
boundaryList[n]:=[];
StabilizerGroups[n]:=[];
RotSubGroups[n]:=[];
  for k in [1..Dimension(n-1)] do
                             #BOTCH JOB!
  if IsFinite(C[n][k].TheMatrixStab) then 
  Append(Elts,Elements(C[n][k].TheMatrixStab));
  else                       
  Append(Elts,GeneratorsOfGroup(C[n][k].TheMatrixStab));
  fi;                        #BOTCH JOB OVER
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

G:=Group(One(Elts[1])*Elts);     #March 2025
G!.bianchiInteger:=dd;

####################
Boundary:=function(n,k);
if k>0 then
return boundaryList[n+1][k];
else
return NegateWord(boundaryList[n+1][-k]);
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
####################

if name="SL2Z" then G:=SL(2,Integers); fi;
SetName(G,"matrix group");

return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:=Elts,
            group:=G,
            stabilizer:=Stabilizer,
            action:=Action,
            properties:=
            [["length",Maximum(1000,lnth)],
             ["characteristic",0],
             ["type","resolution"],
             ["reduced",true]]  ));

end);
################################################
################################################


################################################
################################################
InstallGlobalFunction(RecalculateIncidenceNumbers,
function(R)
local
        NewBoundaryList,
        NewBoundary,
        S,N,pos,StandardWord,
        cnt,n,i,b,e,V,v,L,LL,x,xx;


if IsBound(R!.standardWord) then StandardWord:=R!.standardWord;
else StandardWord:=function(n,bnd); return bnd; end; fi;

NewBoundaryList:=[];
for n in [1..Length(R)] do
NewBoundaryList[n]:=[];
for i in [1..R!.dimension(n)] do
b:=StandardWord(n-1,1*R!.boundary(n,i));
b:=List(b,x->[AbsInt(x[1]),x[2]]);
Add(NewBoundaryList[n],b);
od;
od;


##############################
NewBoundary:=function(n,i);
if i>0 then
return NewBoundaryList[n][i];
else
return NegateWord(NewBoundaryList[n][-i]);
fi;
end;
##############################

########First Incidence Numbers###########
for b in NewBoundaryList[1] do
if b[1][1]*b[2][1]>0 then
b[1][1]:=-1*b[1][1];
fi;
od;
##########################################


R!.oldBoundary:=R!.boundary;
R!.boundary:=NewBoundary;

for N in [2..Length(NewBoundaryList)] do
######## N-th Incidence Numbers##########
cnt:=0;
for e in NewBoundaryList[N] do
  cnt:=cnt+1;
  #########################
  if Length(ResolutionBoundaryOfWord(R,N-1,e))>0 then
  #########################
    V:=List(e,x->[AbsInt(x[1]),x[2]]);
    v:=Random(V); b:=[v]; pos:=Position(V,v); Remove(V,pos);
    #b:=[V[1]]; V:=V{[2..Length(V)]}; 
    while Length(V)>0 do
      L:=ResolutionBoundaryOfWord(R,N-1,b); 
      L:=StandardWord(N-2,L);
      L:=AlgebraicReduction(L);
      LL:=List(L,x->[AbsInt(x[1]),x[2]]);
while L=[] do
 v:=Random(V); 
 Add(b,v); pos:=Position(V,v);Remove(V,pos);  
 L:=ResolutionBoundaryOfWord(R,N-1,b); 
 L:=StandardWord(N-2,L);
 L:=AlgebraicReduction(L);
 LL:=List(L,x->[AbsInt(x[1]),x[2]]);
od;

      for v in V do
        x:=ResolutionBoundaryOfWord(R,N-1,[v]);
        x:=StandardWord(N-2,x);
        xx:=List(x,a->[AbsInt(a[1]),a[2]]);  
        if Length(Intersection(xx,LL))>0 then
          if Length(Intersection(x,L))>0 then Add(b,[-v[1],v[2]]); 
             else Add(b,[v[1],v[2]]); 
          fi;
          pos:=Position(V,v);
          Remove(V,pos);
          break; 
        fi;
      od;
    od;
    NewBoundaryList[N][cnt]:=b;
  #########################
  fi;
  #########################
od;
##########################################
od;

#for N in [2..Length(NewBoundaryList)] do
#for e in NewBoundaryList[N] do
#   e:=StandardWord(N-1,e);
#   e:=StandardWord(N-2,ResolutionBoundaryOfWord(R,N-1,e));
#  if Length(e  )>0 then Print(e,"\n\n");
#  return false; fi;
#od;od;

end);
################################################
################################################	


################################################
################################################
InstallGlobalFunction(QuotientOfContractibleGcomplex,
function(C,S)
local
        Elts,G,Stabilizer,Action,D;

SetInfoLevel(InfoWarning,0);
D:=List(Elements(S),x->x);
Elts:=List(C!.elts, x->QuotientGroup(x,D));
G:=Group(Elts);
#Action:=function(a,b,c) return 1; end;

#####################
Action:=function(n,k,g)
local gg;
gg:=Position(C!.elts,Elts[g]!.element[1]);
if gg=fail then Add(C!.elts,Elts[g]!.element[1]);
gg:=Length(C!.elts);
fi;
return C!.action(n,k,gg);
end;
#####################

#####################
Stabilizer:=function(n,i);
return 
Group(List(Elements(C!.stabilizer(n,i)),x->QuotientGroup(x,D)));
end;
#####################

SetInfoLevel(InfoWarning,1);

return Objectify(HapNonFreeResolution,
            rec(
            dimension:=C!.dimension,
            boundary:=C!.boundary,
            homotopy:=fail,
            elts:=Elts,
            group:=G,
            stabilizer:=Stabilizer,
            action:=Action,
            properties:=
            [["length",EvaluateProperty(C,"length")],
             ["characteristic",0],
             ["type","resolution"],
             ["reduced",C!.dimension(0)=1]]  ));

end);
################################################
################################################


################################################
################################################
InstallOtherMethod(CanonicalRightCosetElement,
"for stabilizer in SL(2,Z) of point at infinite in upper-half plane",
[IsMatrixGroup,IsObject],
function(H,gg)
local g,e;
if not IsBound(H!.gname) then
TryNextMethod();
fi;
if H!.gname="sl2inf" then
#########################
if gg[2][2]<0 then g:=-gg;
else g:=gg; fi;
if g[1][2]=0 then return g; fi;
if g[1][2]>0 then e:=Int(g[2][2]/g[1][2]); 
return [[1,-e],[0,1]]*g; fi;
if g[1][2]<0 then e:=Int(g[2][2]/g[1][2]);
return [[1,e],[0,1]]*g; fi; 

#########################
return g;
fi;
TryNextMethod();
end);
################################################
################################################

###########################################################
###########################################################
InstallGlobalFunction(RecalculateIncidenceNumbers_NonFreeRes,
function(K)
local R, fn, PseudoBounds, boundary;

if K!.dimension(3)>0 then return fail; fi;
R:=FreeGResolution(K,2);
RecalculateIncidenceNumbers(R);

#####################
fn:=function(n,bnd)
local B, b, c;
B:=[];
c:=R!.dimension(n)-K!.dimension(n);

for b in bnd do
if AbsInt(b[1])>c then
Add(B,  [ SignInt(b[1])*(AbsInt(b[1])-c) , b[2]]);
fi;
od;

return B;
end;
#####################

PseudoBounds:=[];
PseudoBounds[1]:=List([1+R!.dimension(1)-K!.dimension(1)..R!.dimension(1)],k->R!.boundary(1,k));
PseudoBounds[2]:=List([1+R!.dimension(2)-K!.dimension(2)..R!.dimension(2)],k->R!.boundary(2,k));

Apply(PseudoBounds[1],bnd->fn(0,bnd));
Apply(PseudoBounds[2],bnd->fn(1,bnd));

#####################
boundary:=function(n,k);
if n<1 or n>2 then return []; fi;

if k>0 then return PseudoBounds[n][k];
else return PseudoBounds[n][-k]; fi;
end;
#####################

K!.boundary:=boundary;

end);
###########################################################
###########################################################

