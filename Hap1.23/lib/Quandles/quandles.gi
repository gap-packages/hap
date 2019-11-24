#####################################################################
#####################################################################
InstallMethod( ViewObj,
 "for HapConjQuandElt",
 [IsHapConjQuandElt],
 function(R)
 Print(R!.element, " ");
end);

#####################################################################
#####################################################################
InstallMethod( PrintObj,
 "for HapConjQuandElt",
 [IsHapConjQuandElt],
 function(R)
 Print(R!.element," ");
 end);

#####################################################################
#####################################################################
InstallGlobalFunction(NumberConnectedQuandles,
function(n)
local L;

if n>47 or n<1 then 
Print("The number of connected quandles of this order is not available.\n");
return fail; fi;

L:=[ 1, 0, 1, 1, 3, 2, 5, 3, 8, 1, 9, 10, 11, 0, 7, 9, 15, 12, 17, 10, 9, 0, 
  21, 42, 34, 0, 65, 13, 27, 24, 29, 17, 11, 0, 15, 73, 35, 0, 13, 33, 39, 
  26, 41, 9, 45, 0, 45 ];;

return L[n];
end);

#####################################################################
#####################################################################
InstallGlobalFunction(CosetsQuandle,
function(G,H,f)
local Trans, T, i, j;

Trans:=RightTransversal(G,H);
T:=List([1..Length(Trans)],i->[]);

for i in [1..Length(T)] do
for j in [1..Length(T)] do
T[i][j]:=PositionCanonical(Trans,f^-1*(Trans[i]*Trans[j]^-1)*f*Trans[j]);   
od;od;



T:=MagmaByMultiplicationTable(T);
if IsQuandle(T) then 
  if Size(T)<=1000 then return AsMagma(Elements(T));
  else return T;
  fi;
else
Print("The data does not yield a quandle.\n");
return fail; fi;
end);


#####################################################################
#####################################################################
InstallGlobalFunction(Cedric_ConjugateQuandleElement,
function(w,n);

return Objectify(HapConjQuandElt,
                    rec(
                    element:=w, power:=n) );
end);

#####################################################################
#####################################################################
InstallGlobalFunction(ConjugationQuandle,

#This will return the conjugation quandle of a finite group G

function(G,n)
local elts;
elts:=List(Elements(G),g->Cedric_ConjugateQuandleElement(g,n));

if Size(elts)<=1000 then
return AsMagma(elts);
else
return MagmaByGenerators(elts);
fi;

end);

#####################################################################
#####################################################################
InstallOtherMethod( \*,
    "composition in a conjugation group quandle",
    [ IsHapConjQuandElt, IsHapConjQuandElt],

function(x,y) local w,u;
u:=y!.element^-(y!.power);
w:=u*x!.element*u^-1;
return Cedric_ConjugateQuandleElement(w,x!.power);
      
end);

#####################################################################
#####################################################################
InstallOtherMethod( \<,
    "comparison in a conjugation group quandle",
    [ IsHapConjQuandElt, IsHapConjQuandElt],

function(x,y) ;

return x!.element<y!.element;

end);

#####################################################################
#####################################################################
InstallOtherMethod( \=,
    "equality in a conjugation group quandle",
    [ IsHapConjQuandElt, IsHapConjQuandElt],

function(x,y) ;

return x!.element=y!.element;

end);

#####################################################################
#####################################################################
InstallGlobalFunction(FirstQuandleAxiomIsSatisfied,
function(M)
local x;

for x in M do
if not x*x=x then return false; fi;
od;
return true;
end);

#####################################################################
#####################################################################
InstallGlobalFunction(SecondQuandleAxiomIsSatisfied,
function(M)
local x,y,z,secondAxSatis, nbInv;

secondAxSatis:=true;
for x in M do
for y in M do
nbInv:=0;
for z in M do
if z*y=x then nbInv:=nbInv+1; fi;
od;
secondAxSatis:=(nbInv=1) and secondAxSatis;
od;od;
return secondAxSatis;
end);

#####################################################################
#####################################################################
InstallGlobalFunction(ThirdQuandleAxiomIsSatisfied,
function(M)
local x,y,z;

for x in M do
for y in M do
for z in M do
if not (x*y)*z=(x*z)*(y*z) then return false; fi;
od; od; od;
return true;
end);

#####################################################################
#####################################################################
InstallMethod(IsQuandle,
"test to see if a magma is a quandle",
[IsMagma],
function(M)

return FirstQuandleAxiomIsSatisfied(M) and SecondQuandleAxiomIsSatisfied(M) and ThirdQuandleAxiomIsSatisfied(M);
end);

#####################################################################
#####################################################################
InstallGlobalFunction(Cedric_CheckThirdAxiomRow,
function(MultMat)
local i,j,k,l,rowA,rowB,rowC;

for i in [1..Length(MultMat)-1] do
for j in [i..Length(MultMat)] do

rowA:=MultMat[i];
rowB:=MultMat[j];
k:=rowB[i];
if k<=Length(MultMat) then
rowC:=MultMat[k];
for l in [1..Length(rowA)] do
if not rowB[rowA[l]]=rowC[rowB[l]] then return false; fi;
od; fi; od; od;
return true;
end);

#####################################################################
#####################################################################
InstallGlobalFunction(Cedric_Permute,
function(permu,permuMat,invPermuMat, Matrix)
local size,i,j,FinalMat;

size:=DimensionsMat(Matrix)[1];
FinalMat:=IdentityMat(size);
for i in [1..size] do
	for j in [1..size] do		
		FinalMat[i][j]:=permu[Matrix[i][j]];			
	od;
od;
FinalMat:=invPermuMat*FinalMat*permuMat;
return FinalMat;
end);

#####################################################################
#####################################################################
InstallGlobalFunction(Cedric_Quandle1,
function(FiltList,PermL,permMatList,iPermMatList, listMinT)
local r1,t;

for r1 in FiltList[1] do
t:=TransposedMat([r1]);
if ThirdQuandleAxiomIsSatisfied(MagmaByMultiplicationTable(t)) then Add(listMinT, t); fi;
od;
end);

#####################################################################
#####################################################################
InstallGlobalFunction(Cedric_Quandle2,
function(FiltList,PermL,permMatList,iPermMatList, listMinT)
local r1,r2,t;

for r1 in FiltList[1] do
for r2 in FiltList[2] do
t:=TransposedMat([r1,r2]);
if ThirdQuandleAxiomIsSatisfied(MagmaByMultiplicationTable(t)) then Add(listMinT, t); fi;
od;od;
end);

#####################################################################
#####################################################################
InstallGlobalFunction(Cedric_Quandle3,
function(FiltList,PermL,permMatList,iPermMatList, listMinT)
local r1,r2,r3,t;

for r1 in FiltList[1] do
for r2 in FiltList[2] do
if Cedric_CheckThirdAxiomRow([r1,r2]) then
for r3 in FiltList[3] do
if Cedric_CheckThirdAxiomRow([r1,r2,r3]) then
t:=TransposedMat([r1,r2,r3]);
if ThirdQuandleAxiomIsSatisfied(MagmaByMultiplicationTable(t)) then Add(listMinT, t); fi;
fi;od;fi;od;od;
end);

#####################################################################
#####################################################################
InstallGlobalFunction(Cedric_Quandle4,
function(FiltList,PermL,permMatList,iPermMatList, listMinT)
local r1,r2,r3,r4,t,i,PT,listPT,minT;

for r1 in FiltList[1] do
for r2 in FiltList[2] do
if Cedric_CheckThirdAxiomRow([r1,r2]) then
for r3 in FiltList[3] do
if Cedric_CheckThirdAxiomRow([r1,r2,r3]) then
for r4 in FiltList[4] do
if Cedric_CheckThirdAxiomRow([r1,r2,r3,r4]) then
t:=TransposedMat([r1,r2,r3,r4]);
if ThirdQuandleAxiomIsSatisfied(MagmaByMultiplicationTable(t)) then 
minT:=t;
for i in [1..Length(permMatList)] do
PT:=Cedric_Permute(PermL[i],permMatList[i],iPermMatList[i],t);
if PT<minT then minT:=PT; fi;
od;
AddSet(listMinT,minT);
fi;
fi;od;fi;od;fi;od;od;
end);

#####################################################################
#####################################################################
InstallGlobalFunction(Cedric_Quandle5,
function(FiltList,PermL,permMatList,iPermMatList,listMinT)
local r1,r2,r3,r4,r5,t,i,PT,listPT,minT;

for r1 in FiltList[1] do
for r2 in FiltList[2] do
if Cedric_CheckThirdAxiomRow([r1,r2]) then
for r3 in FiltList[3] do
if Cedric_CheckThirdAxiomRow([r1,r2,r3]) then
for r4 in FiltList[4] do
if Cedric_CheckThirdAxiomRow([r1,r2,r3,r4]) then
for r5 in FiltList[5] do
if Cedric_CheckThirdAxiomRow([r1,r2,r3,r4,r5]) then
t:=TransposedMat([r1,r2,r3,r4,r5]);
if ThirdQuandleAxiomIsSatisfied(MagmaByMultiplicationTable(t)) then 
minT:=t;
for i in [1..Length(permMatList)] do
PT:=Cedric_Permute(PermL[i],permMatList[i],iPermMatList[i],t);
if PT<minT then minT:=PT; fi;
if minT in listMinT then break; fi;
od;
AddSet(listMinT,minT);
fi;
fi;od;fi;od;fi;od;fi;od;od;
end);

#####################################################################
#####################################################################
InstallGlobalFunction(Cedric_Quandle6,
function(FiltList,PermL,permMatList,iPermMatList,listMinT)
local r1,r2,r3,r4,r5,r6,t,i,PT,listPT,minT,n;

for r1 in FiltList[1] do
for r2 in FiltList[2] do
if Cedric_CheckThirdAxiomRow([r1,r2]) then
for r3 in FiltList[3] do
if Cedric_CheckThirdAxiomRow([r1,r2,r3]) then
for r4 in FiltList[4] do
if Cedric_CheckThirdAxiomRow([r1,r2,r3,r4]) then
for r5 in FiltList[5] do
if Cedric_CheckThirdAxiomRow([r1,r2,r3,r4,r5]) then
for r6 in FiltList[6] do
if Cedric_CheckThirdAxiomRow([r1,r2,r3,r4,r5,r6]) then
t:=TransposedMat([r1,r2,r3,r4,r5,r6]);
if ThirdQuandleAxiomIsSatisfied(MagmaByMultiplicationTable(t)) then 
minT:=t;
for i in [1..Length(permMatList)] do
PT:=Cedric_Permute(PermL[i],permMatList[i],iPermMatList[i],t);
if PT<minT then minT:=PT; fi;
if minT in listMinT then break; fi;
od;
AddSet(listMinT,minT);
fi;
fi;od;fi;od;fi;od;fi;od;fi;od;od;
end);

#####################################################################
#####################################################################

Cedric_XYXYQuandles:=[];
HAPRIGXXX:=[];

#####################################################################
#####################################################################
InstallMethod(Quandles,"for an integer leq to 6",[IsInt],
function(n)
local P,C,peMatList,i,iPeMatList,ListQ;

if n>6 then return fail; fi;
if IsBound(Cedric_XYXYQuandles[n]) then return Cedric_XYXYQuandles[n]; fi;

P:=PermutationsList([1..n]);
C:=List([1..n],i->Filtered(P,p->p[i]=i));
peMatList:=List([1..Length(P)],i->PermutationMat(PermList(P[i]),n));
iPeMatList:=List([1..Length(P)],i->Inverse(peMatList[i]));

ListQ:=[];

if n = 1 then Cedric_Quandle1(C,P,peMatList,iPeMatList,ListQ); fi;
if n = 2 then Cedric_Quandle2(C,P,peMatList,iPeMatList,ListQ); fi;
if n = 3 then Cedric_Quandle3(C,P,peMatList,iPeMatList,ListQ); fi;
if n = 4 then Cedric_Quandle4(C,P,peMatList,iPeMatList,ListQ); fi;
if n = 5 then Cedric_Quandle5(C,P,peMatList,iPeMatList,ListQ); fi;
if n = 6 then Cedric_Quandle6(C,P,peMatList,iPeMatList,ListQ); fi;

Apply(ListQ,MagmaByMultiplicationTable);

Cedric_XYXYQuandles[n]:=ListQ; MakeImmutable("Cedric_XYXYQuandles[n]");
return ShallowCopy(ListQ);
end);

#####################################################################
#####################################################################
InstallMethod(Quandle,"for two integers, the first one leq to 6",[IsInt,IsInt],
function(n,k);

if (n=1 and k>1) or (n=2 and k>1) or (n=3 and k>5) or (n=4 and k>7) or (n=5 and k>22) or (n=6 and k>73) or n>6 then return fail; fi;
return Quandles(n)[k];
end);

#####################################################################
#####################################################################
InstallGlobalFunction(IsQuandleEnvelope,
function(Q,G,e,st);

return (IsTransitive(G,Q)) and (e in Elements(Q)) and (st in Center(Stabilizer(G,e))) and (NormalClosure(G,Group([st]))=G);
end);

#####################################################################
#####################################################################
InstallGlobalFunction(QuandleQuandleEnvelope,
function(Q,G,e,stig)
local x,y,g,yChapSt,multTab;

if not IsQuandleEnvelope(Q,G,e,stig) then return fail; fi;

multTab:=List(Q,x->[]);

for y in Q do
for g in G do
if e^g=y then yChapSt:=stig^g; break; fi;
od;
for x in Q do
multTab[x][y]:=x^yChapSt;
od;
od;
return MagmaByMultiplicationTable(multTab);
end);

#####################################################################
#####################################################################
InstallMethod(IsLatinQuandle,"for a magma",[IsMagma],
function(Q)
local x,y,a,bool;

if not IsQuandle(Q) then TryNextMethod(); fi;

for x in Q do
for y in Q do
bool:=false;
for a in Q do
if x*a=y then bool:=true; break;fi;
od;
if not bool then return false; fi;
od;od;
return true;
end);

#####################################################################
#####################################################################
InstallMethod(IsHomogeneousQuandle,"for a quandle magma",[IsMagma],
function(Q)
local A;

if not IsQuandle(Q) then TryNextMethod(); fi;
if IsConnected(Q) then return true; fi;
A:=AutomorphismGroupQuandle(Q);
if Size(Random(Q)^A)=Size(Q) then return true;
else return false; fi;

end);


#####################################################################
#####################################################################
InstallMethod(IsConnected,"for a magma",[IsMagma],
function(Q)
local i,j,q,L1,L2,X;

if not IsQuandle(Q) then TryNextMethod(); fi;

for i in Q do
L1:=[i];
X:=[i];
while not Size(L1)=0 do
L2:=[];
for q in Q do
for j in L1 do
if not j*q in X then Add(L2, j*q); fi;
od;od;
L1:=ShallowCopy(L2);
X:=Union(X,L1);
od;
if Size(X) < Size(Q) then return false; fi;
od;
return true;
end);

#####################################################################
#####################################################################

Cedric_XYXYConnQuan:=[[MagmaByMultiplicationTable([[1]])]];
for n in [14,22,26,34,38,46] do Cedric_XYXYConnQuan[n]:=[]; MakeImmutable("Cedric_XYXYConnQuan[n]"); od;
Unbind(n);

#####################################################################
#####################################################################
InstallGlobalFunction(ConnectedQuandles,
function(arg)
local n,Stab,ListMultTabConnQuan,Q,i,G,der,Norm,listKSI,st,bool,ksi;

n:=arg[1];
if IsBound(Cedric_XYXYConnQuan[n]) then return Cedric_XYXYConnQuan[n]; fi;
if n>30 and Length(arg)=1 then
Print("Transitive groups of degree >30 are not available in the standard GAP distribution. Use the command  ConnectedQuandles(n, \"rig\"); if n<48.\n");
return fail;
fi;


if Length(arg)>1 then
if n<48 then
ReadPackage("HAP", "lib/Quandles/rig.gi");
if n>30 then
Print("This function uses data from the Rig package for racks, quandles and Nichols algebras due to L. Vendramin.\n");
fi;
if n>30 then
Cedric_XYXYConnQuan[n]:=List(HAPRIGXXX[n],T->AsMagma(Elements(MagmaByMultiplicationTable(TransposedMat(T)))));; 
else
Cedric_XYXYConnQuan[n]:=List(HAPRIGXXX[n],T->AsMagma(Elements(MagmaByMultiplicationTable(T))));;
fi;
MakeImmutable("Cedric_XYXYConnQuan[n]");
return Cedric_XYXYConnQuan[n];
else
Print("The connected quandles of order >47 are not available.\n");
return fail;
fi;
fi;

Stab:=Stabilizer(SymmetricGroup(n),1);
ListMultTabConnQuan:=[];
Q:=[1..n];

for G in AllTransitiveGroups(DegreeAction,n,DerivedSubgroup,IsTransitive,function(g) return g/DerivedSubgroup(g);end,IsCyclic) do 

#Norm:=Normalizer(Stab,G);

listKSI:=[];
for st in Center(Stabilizer(G,1)) do
if Order(NormalClosure(G,Group([st])))=Order(G) then
bool:=true;
for ksi in listKSI do

#
#if IsConjugate(Norm,st,ksi) then bool:=false; break; fi;
## This approach is too slow!!


od;
if bool then
AddSet(ListMultTabConnQuan,MultiplicationTable(QuandleQuandleEnvelope(Q,G,1,st))); Add(listKSI,st); fi;
fi;
od;

od;

Apply(ListMultTabConnQuan,MagmaByMultiplicationTable);

ListMultTabConnQuan:=   #This approach is much faster
QuandleIsomorphismRepresentatives(ListMultTabConnQuan); 

Cedric_XYXYConnQuan[n]:=ListMultTabConnQuan; MakeImmutable("Cedric_XYXYConnQuan[n]");
return ShallowCopy(ListMultTabConnQuan);
end);

#####################################################################
#####################################################################
InstallGlobalFunction(ConnectedQuandle,
function(arg)
local CQ, n, k;

n:=arg[1];
k:=arg[2];
if Length(arg)=2 then
CQ:=ConnectedQuandles(n);
else
CQ:=ConnectedQuandles(n,"rig");
fi;
if CQ=fail then
#Print("The connected quandles of order ",n," are not available. The command ConnectedQuandle(n,k,\"rig\") works for n<48.\n"); 
return fail; fi;
if k<=Length(CQ) then return CQ[k]; fi;
return fail;
end);

#####################################################################
#####################################################################
InstallGlobalFunction(IdQuandle,
function(Q) local size,multTab,pmultTab,P,peMatList,iPeMatList,i,MTQ;

size:=Size(Q);
if size>6 or not IsQuandle(Q) then return [size,fail]; fi;
if Q in Quandles(size) then return [size,Position(Quandles(size),Q)]; fi;
if size=1 then return [1,1]; fi;

multTab:=MultiplicationTable(Q);
P:=PermutationsList([1..size]);;
peMatList:=List([1..Length(P)],i->PermutationMat(PermList(P[i]),size));;
iPeMatList:=List([1..Length(P)],i->Inverse(peMatList[i]));;
MTQ:=ShallowCopy(Quandles(size)); Apply(MTQ,MultiplicationTable);
for i in [1..Length(peMatList)] do
pmultTab:=Cedric_Permute(P[i],peMatList[i],iPeMatList[i],multTab);
if pmultTab in MTQ then return [size,Position(MTQ,pmultTab)]; fi;
od;

return [size,fail];
end);

#####################################################################
#####################################################################
InstallGlobalFunction(IdConnectedQuandle,
function(Q) local size, invariant, C, D, invQ, i;

size:=Size(Q);
if not IsConnected(Q) then return [size,fail]; fi;
if Q in ConnectedQuandles(size) then return [size,Position(ConnectedQuandles(size),Q)]; fi;
if size=1 then return [1,1]; fi;

#################
invariant:=function(Q)
local R;
R:=RightMultiplicationGroupOfQuandleAsPerm(Q);
if HasIdGroup(R) then
return IdGroup(R);
else return AbelianInvariants(R); fi;
end;
#################

invQ:=invariant(Q);

C:=ConnectedQuandles(size);
D:=Filtered([1..Length(C)],i->invariant(C[i])=invQ);

for i in D do
if Length(QuandleIsomorphismRepresentatives([Q,C[i]]))=1 then
return [size,i]; fi;;
od;

end);

#####################################################################
#####################################################################
InstallMethod(Projection, "For semi-direct products of groups", [IsGroup,IsInt],
function(S,n)
local pr,ShomG,GhomS,NhomS,SmapN;

if not IsBound(S!.SemidirectProductInfo) then TryNextMethod(); fi;
if not n in [1,2] then TryNextMethod(); fi;

if n=2 then return Projection(S); fi;

#so n=1
ShomG:=Projection(S);
GhomS:=Embedding(S,1);
NhomS:=Embedding(S,2);

######################
pr:=function(s)
local x ;

x:=Image(GhomS,Image(ShomG,s));
x:=s*x^-1;
return PreImagesRepresentative(NhomS,x);
end;
######################

SmapN:=MappingByFunction(S,Image(Projection(S)),pr);
return SmapN;
end);

#####################################################################
#####################################################################
InstallMethod(RightMultiplicationGroupOfQuandleAsPerm,"for a magma", [IsMagma],
function(M)
local listRx,x,multTab,R;

if not IsQuandle(M) then TryNextMethod(); fi;
multTab:=TransposedMat(MultiplicationTable(M));
listRx:=[];
for x in multTab do
AddSet(listRx,PermList(x));
od;

R:=Group(listRx);
if Size(R)>1 then
return Group(SmallGeneratingSet(R));
else return R; fi;
end);

#####################################################################
#####################################################################
InstallMethod(RightMultiplicationGroupOfQuandle,"for a magma",[IsMagma],
function(M)
local RM;

if not IsQuandle(M) then TryNextMethod(); fi;
RM:=RightMultiplicationGroupOfQuandleAsPerm(M);
RM:=List(GeneratorsOfGroup(RM),a-> MagmaHomomorphismByFunctionNC(M,M,function(q) return MagmaElement(M,Position(Elements(M),q)^a); end)); 
return Group(RM);
end);

#####################################################################
#####################################################################
InstallMethod(InnerAutomorphismGroupQuandle,"for a magma",[IsMagma],
function(M)
return RightMultiplicationGroupOfQuandle(M);
end);

#####################################################################
#####################################################################
InstallMethod(InnerAutomorphismGroupQuandleAsPerm,"for a magma",[IsMagma],
function(M)
return RightMultiplicationGroupOfQuandleAsPerm(M);
end);





#####################################################################
#####################################################################
InstallGlobalFunction(Cedric_FromAutGeReToAutQe,
function(Ralpha,RightMultGrpOfQ,Q)
local L,g,A,elts;

elts:=Elements(RightMultGrpOfQ);
A:=List(elts,g->1^g);
L:=List([1..Size(Q)],function(x) g:=elts[Position(A,x)]; return 1^(g^Ralpha); end);
return Inverse(PermList(L));
end);

#####################################################################
#####################################################################
InstallMethod(AutomorphismGroupQuandleAsPerm,"for a magma",[IsMagma],
function(Q)
local G,Re,Ge,AutG,elAutGeRe,AutGeRe,semiDP,AutGeReHomSDP,GHomSDP,L,p,sour,ImL,x,pos,tau,ImageTau,P1,P2,action,AutQ;

if not  IsQuandle(Q) then TryNextMethod(); fi;

if not IsConnected(Q) then
return AutomorphismGroupQuandleAsPerm_nonconnected(Q); fi;

G:=RightMultiplicationGroupOfQuandleAsPerm(Q);;
Re:=PermList(TransposedMat(MultiplicationTable(Q))[1]);;
Ge:=Stabilizer(G,1);;
AutG:=AutomorphismGroup(G);;
elAutGeRe:=Filtered(Elements(AutG),phi->Re^phi=Re and OnSets(Elements(Ge),phi)=Ge);;
AutGeRe:=AsGroup(elAutGeRe);;

semiDP:=SemidirectProduct(AutGeRe,G);;

AutGeReHomSDP:=Embedding(semiDP,1);;
GHomSDP:=Embedding(semiDP,2);;

L:=List(Elements(Ge),a->[a,InnerAutomorphism(Ge,a)^-1]);

p:=Maximum(Elements(Source(AutGeReHomSDP)));
while IsMapping(p) do p:=Maximum(Elements(Source(p))); od;

sour:=List(Elements(Source(AutGeReHomSDP)),s->p^s);
ImL:=[];;
for x in L do
pos:=Position(sour,p^x[2]);
Add(ImL,Image(GHomSDP,x[1])*Elements(Image(AutGeReHomSDP))[pos]);
od;

tau:=NaturalHomomorphismByNormalSubgroup(semiDP,Core(semiDP,Group(ImL)));

ImageTau:=Image(tau);
P1:=Projection(semiDP,1);
P2:=Projection(semiDP,2);

########
action:=function(a,q) local atilde,alpha,phi;

atilde:=PreImagesRepresentative(tau,a);	# in GxAutGeRe
alpha:=ImageElm(P1,atilde);		# in G=RightMultiplicationGroup(Q)
phi:=ImageElm(P2,atilde);		# in AutGeRe

return (q^(Cedric_FromAutGeReToAutQe(phi,G,Q)))^alpha;
end;
########

AutQ:=List(GeneratorsOfGroup(ImageTau),function(a) L:=List([1..Size(Q)],i->action(a,i)); return PermList(L); end);

#AutQ:=List(ImageTau,function(a) L:=List([1..Size(Q)],i->action(a,i)); return PermList(L); end);

if Size(Group(AutQ))>1 then
return Group(SmallGeneratingSet(Group(AutQ)));
else
return Group(AutQ); fi;
end);

#####################################################################
#####################################################################
InstallMethod(AutomorphismGroupQuandle,"for a magma",[IsMagma],
function(Q)
local Au;

if not IsQuandle(Q) then TryNextMethod(); fi;
#if IsConnected(Q) then TryNextMethod(); fi;

Au:=AutomorphismGroupQuandleAsPerm(Q);
Au:=GeneratorsOfGroup(Au); 
Au:=List(Au,a-> MagmaHomomorphismByFunctionNC(Q,Q,function(q) return MagmaElement(Q,Position(Elements(Q),q)^a); end)); 

return Group(Au);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(AutomorphismGroupQuandleAsPerm_nonconnected,
function(Q)
local S, T, Au, s, bool, x, y;

S:=SymmetricGroup(Size(Q));
T:=MultiplicationTable(Q);   #Changed S to Q!!!!
Au:=[];
for s in S do
bool:=true;
for x in [1..Size(Q)] do
for y in [1..Size(Q)] do
if not T[x][y]^s= T[x^s][y^s]
then bool:=false; break; fi;
od;
if not bool then break; fi;
od;
if bool then Add(Au,s); fi;
od;

if Size(Group(Au))>1 then
Au:=Group(SmallGeneratingSet(Group(Au)));
else
Au:=Group(Au);
fi;
return Au;

end);

#####################################################################
#####################################################################
InstallGlobalFunction(AdjointGroupOfQuandle,
function(Q)
local rel,F,G,Gen;

F:=FreeGroup(Size(Q));
Gen:=GeneratorsOfGroup(F);

rel:=List([1..Size(Q)],function(k)
local L,i,j;
L:=[];
for i in [1..Size(Q)] do
for j in [1..Size(Q)] do
if MultiplicationTable(Q)[i][j]=k then Add(L,Gen[j]^-1*Gen[i]*Gen[j]*Gen[k]^-1); fi;
od; od;
return L;
end);;

return F/Flat(rel);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallMethod(PathComponents,
"Path components of a quandle",
[IsMagma],
function(Q) 
local act, R;

if not IsQuandle(Q) then
Print("Magma is not a quandle.\n");
return fail; fi;

R:=RightMultiplicationGroupOfQuandle(Q);
act:=function(q,r) return Image(r,q); end;
return Orbits(R,Q,act);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallMethod(ConnectedComponentsQuandle,
"components of a quandle as quandles",
[IsMagma],
function(Q)
local act, C;

if not IsQuandle(Q) then
Print("Magma is not a quandle.\n");
return fail; fi;

C:=PathComponents(Q);
return List(C,AsMagma);
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod(FundamentalGroup,
"fundamental group of a connected quandle",
[IsMagma],

function(Q)
local A, Inn, AhomInn, S,gensA,gensAI, elts, Z,gensZ,H,R,RR,T,AhomZ,
      D,ZhomA,i,j,a,q;

if not IsConnected(Q) then
Print("The quandle is not connected.\n");
return fail;
fi;

q:=Elements(Q)[1];
Inn:=RightMultiplicationGroupOfQuandle(Q);
Z:=FreeGroup(1);
gensZ:=GeneratorsOfGroup(Z);
A:=AdjointGroupOfQuandle(Q);
gensA:=GeneratorsOfGroup(A); #gensA[i] corresponds to Elements(Q)[i]
gensAI:=List(Elements(Q),q->MagmaHomomorphismByFunctionNC(Q,Q,x->x*q));
AhomInn:=GroupHomomorphismByImagesNC(A,Inn,gensA,gensAI);

AhomZ:=GroupHomomorphismByImagesNC(A,Z,gensA,
       List([1..Length(gensA)],i->gensZ[1]));

ZhomA:=GroupHomomorphismByImagesNC(Z,A,gensZ,[gensA[1]]);

H:=Subgroup(A,[gensA[1]]);   #H=ZZ and A=R x| H with R finite.
S:=RightTransversal(A,H);
R:=List(S, x->Image(ZhomA,Image(AhomZ,x))^-1*x);
RR:=Filtered([1..Length(R)],r->Image(Image(AhomInn,R[r]),q)=q);

R:=R{RR};
T:=[];
for i in [1..Length(R)] do
T[i]:=[];
for j in [1..Length(R)] do
a:=PositionCanonical(S,R[i]*R[j]);
T[i][j]:=Position(RR,a);
od;od;

D:= GroupByMultiplicationTable(T);
D:=SmallGeneratingSet(D);
return Group(D);
end);

##########################################################
##########################################################
InstallGlobalFunction(DerivedGroupOfQuandle,
function(Q)
local A, S,gensA, elts, Z,gensZ,H,R,RR,T,AhomZ,              
      D,ZhomA,i,j,a;

if not IsConnected(Q) then
Print("The quandle is not connected.\n");
return fail;
fi;

Z:=FreeGroup(1);
gensZ:=GeneratorsOfGroup(Z);
A:=AdjointGroupOfQuandle(Q);
gensA:=GeneratorsOfGroup(A); #gensA[i] corresponds to Elements(Q)[i]

AhomZ:=GroupHomomorphismByImagesNC(A,Z,gensA,
       List([1..Length(gensA)],i->gensZ[1]));

ZhomA:=GroupHomomorphismByImagesNC(Z,A,gensZ,[gensA[1]]);

H:=Subgroup(A,[gensA[1]]);   #H=ZZ and A=R x| H with R finite.
S:=RightTransversal(A,H);
R:=List(S, x->Image(ZhomA,Image(AhomZ,x))^-1*x);

T:=[];
for i in [1..Length(R)] do
T[i]:=[];
for j in [1..Length(R)] do
T[i][j]:=PositionCanonical(S,R[i]*R[j]);
od;od;

D:= GroupByMultiplicationTable(T);
D:=SmallGeneratingSet(D);
return Group(D);

end);
##########################################################
##########################################################



