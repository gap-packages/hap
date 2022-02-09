


#####################################################
#####################################################
InstallOtherMethod(NonabelianTensorProduct,
"Nonabelian tensor product of two crossed modules",
[IsHapCrossedModule,IsHapCrossedModule],
function(MG,NG)
local M, N, G, F, MhomG, NhomG, gensM, gensN, gensMf, gensNf, gensF, 
      gensF1, gensF2, MhomMfp, NhomNfp, MfhomF, NfhomF, Mfp, Nfp,
      MhomF, NhomF, FhomG, relsF,r,m,x,y,z,zz,v,w, FhomT, T,actM, actN,
      TN, TM, ShomTT, TT, tensor, delta, gensTT, gens, TThomG, action,
      S, ThomS, X, Y ;;

if not (IsHapCrossedModule(MG) and IsHapCrossedModule(NG)) 
then Print("Arguments must be crossed modules"); return fail; fi;
if not Target(MG!.map)=Target(NG!.map) then
Print("Arguments must have common codomain"); return fail; fi;
if not 
(IsNilpotent(Source(CatOneGroupByCrossedModule(MG)!.sourceMap))
and
IsNilpotent(Source(CatOneGroupByCrossedModule(MG)!.sourceMap)) )
then
Print("WARNING: The largest nilpotent quotient of the tensor product will be returned.\n\n");
fi;


MhomG:=MG!.map;
NhomG:=NG!.map;
actM:=MG!.action;
actN:=NG!.action;
M:=Source(MhomG);
N:=Source(NhomG);
G:=Target(MhomG);

MhomMfp:=IsomorphismFpGroup(M);#
NhomNfp:=IsomorphismFpGroup(N);#

Mfp:=Range(MhomMfp);
Nfp:=Range(NhomNfp);

gensM:=List(GeneratorsOfGroup(Mfp),x->PreImagesRepresentative(MhomMfp,x));#
gensN:=List(GeneratorsOfGroup(Nfp),x->PreImagesRepresentative(NhomNfp,x));#
gensMf:=FreeGeneratorsOfFpGroup(Mfp);
gensNf:=FreeGeneratorsOfFpGroup(Nfp);
F:=FreeGroup(Length(gensM)+Length(gensN));
gensF:=GeneratorsOfGroup(F);
gensF1:=gensF{[1..Length(gensM)]};
gensF2:=gensF{Length(gensM)+[1..Length(gensN)]};
FhomG:=GroupHomomorphismByImagesNC(F,G,gensF,
Concatenation(List(gensM,x->Image(MhomG,x)),
List(gensN,x->Image(NhomG,x))));

MfhomF:=GroupHomomorphismByImagesNC(Group(gensMf),F,gensMf,gensF1);
NfhomF:=GroupHomomorphismByImagesNC(Group(gensNf),F,gensNf,gensF2);
MhomF:=GroupHomomorphismByImagesNC(M,F,gensM,gensF1);
NhomF:=GroupHomomorphismByImagesNC(N,F,gensN,gensF2);



relsF:=[];
for r in RelatorsOfFpGroup(Mfp) do
Add(relsF,Image(MfhomF,r));
od;
for r in RelatorsOfFpGroup(Nfp) do
Add(relsF,Image(NfhomF,r));
od;

X:=Kernel(MhomG);
Y:=MinimalGeneratingSet(X);
X:=RightTransversal(M,X);
X:=Concatenation(X,Y);

for z in X do
zz:=Image(MhomG,z);
for x in gensM do
for y in gensN do
v:=Comm(Image(MhomF,x),Image(NhomF,y))^Image(MhomF,z) ;
w:=Comm(Image(NhomF,actN(zz,y)),Image(MhomF,actM(zz,x)) );
Append(relsF,[v*w]);
od;od;od;

X:=Kernel(NhomG);
Y:=MinimalGeneratingSet(X);
X:=RightTransversal(N,X);
X:=Concatenation(X,Y);

for z in X do
zz:=Image(NhomG,z);
for x in gensM do
for y in gensN do
v:=Comm(Image(MhomF,x),Image(NhomF,y))^Image(NhomF,z) ;
w:=Comm(Image(NhomF,actN(zz,y)),Image(MhomF,actM(zz,x)) );
Append(relsF,[v*w]);
od;od;od;


T:=F/relsF;
FhomT:=GroupHomomorphismByImagesNC(F,T,gensF,GeneratorsOfGroup(T));;

ThomS:=IsomorphismSimplifiedFpGroup(T);
S:=Image(ThomS);

ShomTT:=NqEpimorphismNilpotentQuotient(S);

TT:=Image(ShomTT);

TM:=List(gensF1,x->Image(ShomTT,Image(ThomS,Image(FhomT,x))));;
TM:=Group(TM);
TM:=NormalClosure(TT,TM);
TN:=List(gensF2,x->Image(ShomTT,Image(ThomS,Image(FhomT,x))));;
TN:=Group(TN);
TN:=NormalClosure(TT,TN);
tensor:= Intersection(TM,TN); 

gensTT:=GeneratorsOfGroup(TT);
gens:=List(gensTT,x->PreImagesRepresentative(ShomTT,x));
gens:=List(gens,x->PreImagesRepresentative(ThomS,x));
gens:=List(gens,x->PreImagesRepresentative(FhomT,x));
gens:=List(gens,x->Image(FhomG,x));
TThomG:=GroupHomomorphismByImagesNC(TT,G,gensTT,gens);

delta:=GroupHomomorphismByFunction(tensor,G,x->Image(TThomG,x));
########################
action:=function(g,t)
local gg;
gg:=PreImagesRepresentative(TThomG,g);
return gg^-1*t*gg;
end;
########################


return Objectify(HapCrossedModule,
       rec(map:=delta,action:=action));
end);
#####################################################
#####################################################




#####################################################
#####################################################
InstallOtherMethod(NonabelianExteriorProduct,
"Nonabelian exterior product of two crossed modules",
[IsHapCrossedModule,IsHapCrossedModule],
function(MG,NG)
local M, N, G, F, MhomG, NhomG, gensM, gensN, gensMf, gensNf, gensF,
      gensF1, gensF2, MhomMfp, NhomNfp, MfhomF, NfhomF, Mfp, Nfp,
      MhomF, NhomF, FhomG, relsF,r,m,x,y,z,zz,v,w, FhomT, T,actM, actN,
      TN, TM, ShomTT, TT, tensor, delta, gensTT, gens, TThomG, action,
      S, ThomS, X, Y, imM, imN, MN,PreMN,KerM, KerN ;;

if not (IsHapCrossedModule(MG) and IsHapCrossedModule(NG))
then Print("Arguments must be crossed modules"); return fail; fi;
if not Target(MG!.map)=Target(NG!.map) then
Print("Arguments must have common codomain"); return fail; fi;



MhomG:=MG!.map;
NhomG:=NG!.map;
imM:=Image(MhomG);
imN:=Image(NhomG);
KerM:=Kernel(MhomG);
KerN:=Kernel(NhomG);
KerM:=MinimalGeneratingSet(KerM);
KerN:=MinimalGeneratingSet(KerN);
MN:=Intersection(imM,imN);
PreMN:=List(MN,x->[PreImagesRepresentative(MhomG,x),PreImagesRepresentative(NhomG,x)]);
actM:=MG!.action;
actN:=NG!.action;
M:=Source(MhomG);
N:=Source(NhomG);
G:=Target(MhomG);
#gensM:=MinimalGeneratingSet(M);
#gensN:=MinimalGeneratingSet(N);
#MhomMfp:=IsomorphismFpGroupByGenerators(M,gensM);
#NhomNfp:=IsomorphismFpGroupByGenerators(N,gensN);

MhomMfp:=IsomorphismFpGroup(M);#
NhomNfp:=IsomorphismFpGroup(N);#

Mfp:=Range(MhomMfp);
Nfp:=Range(NhomNfp);

gensM:=List(GeneratorsOfGroup(Mfp),x->PreImagesRepresentative(MhomMfp,x));#
gensN:=List(GeneratorsOfGroup(Nfp),x->PreImagesRepresentative(NhomNfp,x));#
gensMf:=FreeGeneratorsOfFpGroup(Mfp);
gensNf:=FreeGeneratorsOfFpGroup(Nfp);
F:=FreeGroup(Length(gensM)+Length(gensN));
gensF:=GeneratorsOfGroup(F);
gensF1:=gensF{[1..Length(gensM)]};
gensF2:=gensF{Length(gensM)+[1..Length(gensN)]};
FhomG:=GroupHomomorphismByImagesNC(F,G,gensF,
Concatenation(List(gensM,x->Image(MhomG,x)),
List(gensN,x->Image(NhomG,x))));

MfhomF:=GroupHomomorphismByImagesNC(Group(gensMf),F,gensMf,gensF1);
NfhomF:=GroupHomomorphismByImagesNC(Group(gensNf),F,gensNf,gensF2);
MhomF:=GroupHomomorphismByImagesNC(M,F,gensM,gensF1);
NhomF:=GroupHomomorphismByImagesNC(N,F,gensN,gensF2);


relsF:=[];
for r in RelatorsOfFpGroup(Mfp) do
Add(relsF,Image(MfhomF,r));
od;
for r in RelatorsOfFpGroup(Nfp) do
Add(relsF,Image(NfhomF,r));
od;

X:=Kernel(MhomG);
Y:=MinimalGeneratingSet(X);
X:=RightTransversal(M,X);
X:=Concatenation(X,Y);

for z in X do
zz:=Image(MhomG,z);
for x in gensM do
for y in gensN do
v:=Comm(Image(MhomF,x),Image(NhomF,y))^Image(MhomF,z) ;
w:=Comm(Image(NhomF,actN(zz,y)),Image(MhomF,actM(zz,x)) );
Append(relsF,[v*w]);
od;od;od;

X:=Kernel(NhomG);
Y:=MinimalGeneratingSet(X);
X:=RightTransversal(N,X);
X:=Concatenation(X,Y);

for z in X do
zz:=Image(NhomG,z);
for x in gensM do
for y in gensN do
v:=Comm(Image(MhomF,x),Image(NhomF,y))^Image(NhomF,z) ;
w:=Comm(Image(NhomF,actN(zz,y)),Image(MhomF,actM(zz,x)) );
Append(relsF,[v*w]);
od;od;od;

for x in PreMN do
v:=Comm(Image(MhomF,x[1]),Image(NhomF,x[2]));
Add(relsF,v);
for y in KerN do
v:=Comm(Image(MhomF,x[1]),Image(NhomF,y));
Add(relsF,v);
od;
for y in KerM do
v:=Comm(Image(MhomF,y),Image(NhomF,x[2]));
Add(relsF,v);
od;


od;


T:=F/relsF;
FhomT:=GroupHomomorphismByImagesNC(F,T,gensF,GeneratorsOfGroup(T));;

ThomS:=IsomorphismSimplifiedFpGroup(T);
S:=Image(ThomS);

ShomTT:=NqEpimorphismNilpotentQuotient(S);
TT:=Image(ShomTT);


TM:=List(gensF1,x->Image(ShomTT,Image(ThomS,Image(FhomT,x))));;
TM:=Group(TM);
TM:=NormalClosure(TT,TM);
TN:=List(gensF2,x->Image(ShomTT,Image(ThomS,Image(FhomT,x))));;
TN:=Group(TN);
TN:=NormalClosure(TT,TN);
tensor:= Intersection(TM,TN);

gensTT:=GeneratorsOfGroup(TT);
gens:=List(gensTT,x->PreImagesRepresentative(ShomTT,x));
gens:=List(gens,x->PreImagesRepresentative(ThomS,x));
gens:=List(gens,x->PreImagesRepresentative(FhomT,x));
gens:=List(gens,x->Image(FhomG,x));
TThomG:=GroupHomomorphismByImagesNC(TT,G,gensTT,gens);

delta:=GroupHomomorphismByFunction(tensor,G,x->Image(TThomG,x));
########################
action:=function(g,t)
local gg;
gg:=PreImagesRepresentative(TThomG,g);
return gg^-1*t*gg;
end;
########################


return Objectify(HapCrossedModule,
       rec(map:=delta,action:=action));
end);
#####################################################
#####################################################

