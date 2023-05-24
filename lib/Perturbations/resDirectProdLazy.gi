#(C) Graham Ellis, 2005-2006

#################################################################
#################################################################
InstallGlobalFunction(ResolutionFiniteCyclicGroup,
function(MM,N)
local G, GG, M, iso, R, Dimension, Boundary, Homotopy, fn, posfn, pcp, Elts,g;

if IsInt(MM) then M:=MM;
G:=AbelianPcpGroup([M]);
else
if not IsCyclic(MM) then 
Print("The group must be finite cyclic.\n");
return fail; fi;
M:=Order(MM);
G:=AbelianPcpGroup([M]);
iso:=GroupHomomorphismByImagesNC(G,MM,MinimalGeneratingSet(G),MinimalGeneratingSet(MM));
fi;
pcp:=Pcp(G);
g:=MinimalGeneratingSet(G);
g:=g[1];

##########################
##########################
Dimension:=function(n);
if n>=0  then return 1;
else return 0;
fi;
end;
##########################
##########################

##########################
##########################
Boundary:=function(n,k)
if n<=0 or AbsInt(k)>1 then return []; fi;

if IsOddInt(n) then 
if k>0 then return [[1,2], [-1,1]]; 
else return [[-1,2], [1,1]]; fi;
fi;

if IsEvenInt(n) then 
if k>0 then return 
List([1..M],i->[1,i]); 
else return List([1..M],i->[-1,i]); fi;
fi;
end;
##########################
##########################

##########################
##########################
fn:=function(k);
return g^(k-1);
end;
Elts:=LazyList(fn,[["length",M]]);

posfn:=function(x)
local e;
e:=ExponentsByPcp(pcp,x)[1];
if e>0 then return e+1; else return 1; fi;
end;
Elts!.posfun:=posfn;
##########################
##########################

if IsGroup(MM) then 
##########################
##########################
fn:=function(k);
return Image(iso,g^(k-1));
end;
Elts:=LazyList(fn,[["length",M]]);

posfn:=function(x)
local e,y;
y:=PreImagesRepresentative(iso,x);
e:=ExponentsByPcp(pcp,y)[1];
if e>0 then return e+1; else return 1; fi;
end;
Elts!.posfun:=posfn;
##########################
##########################
fi;

##########################
##########################
Homotopy:=function(n,x)
local k;
k:=x[2];
if IsEvenInt(n) then
if x[1]>0 then return List([1..k-1],i->[1,i]);
else return List([1..k-1],i->[-1,i]);fi;
else
if not k=M then return [];
else
if x[1]>0 then
return [[1,1]];
else return [[-1,1]]; fi;
fi;
fi;

end;
##########################
##########################

if IsInt(MM) then GG:=G;
else GG:=Image(iso); fi;

return Objectify(HapResolution,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                homotopy:=Homotopy,
                elts:=Elts,
                group:=GG,
                properties:=
                   [["length",N],
                    ["reduced",true],
                    ["type","resolution"],
                    ["characteristic",0]  ]));

end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(ResolutionInfiniteCyclicGroup,
function(arg)
local N,G, Dimension, Boundary, Homotopy, fn, posfn, pcp, Elts,g;

if Length(arg)=1 then N:=arg[1]; 
G:=AbelianPcpGroup([0]);
pcp:=Pcp(G);
g:=GeneratorsOfGroup(G);
g:=g[1];
fi;

if Length(arg)=2 then N:=arg[2]; 
G:=arg[1];
pcp:=Pcp(G);
g:=GeneratorsOfGroup(G);
g:=g[1];
fi;

##########################
##########################
Dimension:=function(n);
if n=0 or n=1 then return 1;
else return 0;
fi;
end;
##########################
##########################

##########################
##########################
Boundary:=function(n,k);
if n=1 and k=1 then return [[1,2], [-1,1]];
fi;
return [];
end;
##########################
##########################

##########################
##########################
fn:=function(k) local kk;
if IsOddInt(k) then kk:=(1-k)/2;
return g^kk;
else kk:=k/2;
return g^kk;
fi;
end;
Elts:=LazyList(fn);

posfn:=function(x)
local e;
e:=ExponentsByPcp(pcp,x)[1];
if e<=0 then return 1-2*e; fi;
 return 2*e; ;
end;
Elts!.posfun:=posfn;
##########################
##########################

##########################
##########################
Homotopy:=function(n,x)
local k,kk;
if not n=0 then return []; fi;

k:=x[2];
if IsEvenInt(k) then
kk:=k/2; kk:=2*[0..kk-1]; kk[1]:=1;
return List(kk,i->[1,i]);
else
kk:=(1-k)/2; kk:=2*[2..1-kk]; kk:=kk-1;
return List(kk,i->[-1,i]);
fi;
end;
##########################
##########################

return Objectify(HapResolution,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                homotopy:=Homotopy,
                elts:=Elts,
                group:=G,
                properties:=
                   [["length",N],
                    ["reduced",true],
                    ["type","resolution"],
                    ["characteristic",0]  ]));

end);
#################################################################
#################################################################



#####################################################################
InstallGlobalFunction(ResolutionDirectProductLazy,
function(arg)
local
	R,S,T,
	G,H,E,K,GhomE,HhomE,EhomG,EhomH,EltsE,
	pcpEmod,Emod,Ehom, homE,
 	DimensionR,BoundaryR,HomotopyR,
        DimensionS,BoundaryS,HomotopyS,
	Lngth,Dimension,Boundary,Homotopy,
	PseudoBoundary,
	DimPQ,
	Int2Pair, Pair2Int,
	Charact,
	AddWrds,
	Int2Vector,
	Vector2Int,
	Elts2Int,
	HomotopyRec,
	HomotopyGradedGen,
	HomotopyOfWord,
	FinalHomotopy,
	HorizontalBoundaryGen,
	HorizontalBoundaryWord,
	F,FhomE,
	gensE, gensE1, gensE2, gensG, gensH, gens2, gens1,
	Boole, 
	i,j,k,g,h,fn,posfn;

R:=arg[1];
S:=arg[2];
G:=R!.group; 
H:=S!.group; 

####################### DIRECT PRODUCT OF GROUPS ###########
E:=DirectProduct(G,H);
if Size(G)=infinity or Size(H)=infinity then SetSize(E,infinity);fi;
GhomE:=Embedding(E,1);
HhomE:=Embedding(E,2);
EhomG:=Projection(E,1);
EhomH:=Projection(E,2);

#if Parent(G)=Parent(H) and not G=H  then
if One(G) in H and not G=H then

gensG:=MinimalGeneratingSet(G);
gens1:=List(gensG,x->One(G));
gensH:=MinimalGeneratingSet(H);
gens2:=List(gensH,x->One(H));
gensE:=Concatenation(gensG,gensH);
gensE2:=Concatenation(gensG,gens2);
gensE1:=Concatenation(gens1,gensH);
E:=Group(gensE);

GhomE:=GroupHomomorphismByImages(G,E,gensG,gensG);
HhomE:=GroupHomomorphismByImages(H,E,gensH,gensH);
EhomG:=GroupHomomorphismByImages(E,G,gensE,gensE2);
EhomH:=GroupHomomorphismByImages(E,H,gensE,gensE1);
fi;

################ DIRECT PRODUCT OF GROUPS CONSTRUCTED #########


if not (IsFinite(G) or IsFinite(H)) then
#########################
fn:=function(n)
local d,r;
d:=1;
while (d+1)*d/2<=n do
d:=d+1;
od;
d:=d-1;
r:=n-((d+1)*d/2);
if r=0 then return Image(GhomE,R!.elts[1])*Image(HhomE,S!.elts[d]); fi;
return Image(GhomE,R!.elts[d+2-r])*Image(HhomE,S!.elts[r]);
end;
#EltsE:=LazyList(fn,[["length",Order(G)*Order(H)]]);
EltsE:=LazyList(fn,[["length",infinity]]);
#EltsE:=LazyList(fn);


posfn:=function(x)
local g,h, posg,posh,d;
g:=Image(EhomG,x);
h:=Image(EhomH,x);
posg:=Position(R!.elts,g);
posh:=Position(S!.elts,h);
d:=posg+posh-2;
return ((d+1)*d/2) + posh ;
end;
EltsE!.posfun:=posfn;  #This doesn't seem to affect the runtime.
#########################
fi;

if IsFinite(H) then
#########################
fn:=function(n)
local d,r;
d:=Int(n/Order(H));
r:=n mod Order(H);
if r=0 then return Image(GhomE,R!.elts[d])*Image(HhomE,S!.elts[Order(H)]); fi;
return Image(GhomE,R!.elts[d+1])*Image(HhomE,S!.elts[r]);
end;
EltsE:=LazyList(fn, [["length",Size(H)*Size(G)]]);

posfn:=function(x)
local g,h, posg,posh,d;
g:=Image(EhomG,x);
h:=Image(EhomH,x);
posg:=Position(R!.elts,g);
posh:=Position(S!.elts,h);
return (posg -1)*Order(H)+posh ;
end;
EltsE!.posfun:=posfn;  #This doesn't seem to affect the runtime.
#########################
fi;

if IsFinite(G) and not IsFinite(H) then
#########################
fn:=function(n)
local d,r;
d:=Int(n/Order(G));
r:=n mod Order(G);
if r=0 then return Image(GhomE,R!.elts[Order(G)])*Image(HhomE,S!.elts[d]); fi;
return Image(GhomE,R!.elts[r])*Image(HhomE,S!.elts[d+1]);
end;
EltsE:=LazyList(fn);

posfn:=function(x)
local g,h, posg,posh,d;
g:=Image(EhomG,x);
h:=Image(EhomH,x);
posg:=Position(R!.elts,g);
posh:=Position(S!.elts,h);
return (posh -1)*Order(G)+posg ;
end;
EltsE!.posfun:=posfn;  #This doesn't seem to affect the runtime.
#########################
fi;



PseudoBoundary:=[];
#########################
DimensionR:=function(n);
if n<0 or n>Length(R) then return 0;
else return R!.dimension(n);
fi;
end;; 
#########################
#########################
DimensionS:=function(n);
if n<0 or n>Length(S) then return 0;
else return S!.dimension(n);
fi;
end;;
#########################

BoundaryS:= S!.boundary;
	   
BoundaryR:=R!.boundary;  
HomotopyR:=R!.homotopy;
HomotopyS:=S!.homotopy;  

#################DETERMINE VARIOUS PROPERTIES########################
Lngth:=Minimum(EvaluateProperty(R,"length"),EvaluateProperty(S,"length"));
#Lngth:=Length(R)+Length(S);

for i in [3..Length(arg)] do
if IsInt(arg[i]) then Lngth:=arg[i]; fi;
od;

if EvaluateProperty(R,"characteristic")=0
and EvaluateProperty(S,"characteristic")=0
then Charact:=EvaluateProperty(R,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")=0
and EvaluateProperty(S,"characteristic")>0
then Charact:=EvaluateProperty(S,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")>0
and EvaluateProperty(S,"characteristic")=0
then Charact:=EvaluateProperty(R,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")>0
and EvaluateProperty(S,"characteristic")>0
then Charact:=Product(Intersection([
DivisorsInt(EvaluateProperty(R,"characteristic")),
DivisorsInt(EvaluateProperty(S,"characteristic"))
]));
fi;

if Charact=0 then AddWrds:=AddFreeWords; else
        AddWrds:=function(v,w);
        return AddFreeWordsModP(v,w,Charact);
        end;
fi;
####################PROPERTIES DETERMINED############################


#####################################################################
Dimension:=function(i)
local D,j;
if i<0 then return 0; fi;
if i=0 then return 1; else
D:=0;

for j in [0..i] do
D:=D+DimensionR(j)*DimensionS(i-j);
od;

return D; fi;
end;
#####################################################################

for i in [1..Lngth] do
PseudoBoundary[i]:=[1..Dimension(i)];
od;

#####################################################################
DimPQ:=function(p,q)
local D,j;

if (p<0) or (q<0) then return 0; fi;

D:=0;
for j in [0..q] do
D:=D+DimensionR(p+q-j)*DimensionS(j);
od;

return D;
end;
#####################################################################

#####################################################################
Int2Pair:=function(i,p,q)       #Assume that x<=DimR(p)*DimS(q).
local s,r,x;
                                #The idea is that the generator f_i in F
				#corresponds to a tensor (e_r x e_s)
x:=AbsInt(i)-DimPQ(p+1,q-1);     #with e_r in R_p, e_s in S_q. If we
s:= x mod DimensionS(q);                #input i we get output [r,s].
r:=(x-s)/DimensionS(q);

if s=0 then return [SignInt(i)*r,DimensionS(q)];
else return [SignInt(i)*(r+1),s]; fi;

end;
#####################################################################

#####################################################################
Pair2Int:=function(x,p,q)
local y;                        #Pair2Int is the inverse of Int2Pair.

y:=[AbsInt(x[1]),AbsInt(x[2])];
return SignInt(x[1])*SignInt(x[2])*((y[1]-1)*DimensionS(q)+y[2]+DimPQ(p+1,q-1));end;
#####################################################################

#####################################################################
Int2Vector:=function(k,j)
local tmp,p,q;

p:=k;q:=0;
while j>=DimPQ(p,q)+1 do
p:=p-1;q:=q+1;
od;				#p,q are now computed from k,j

tmp:=Int2Pair(j,p,q);

return [p,q,tmp[1],tmp[2]];
end;
#####################################################################

#####################################################################
Vector2Int:=function(p,q,r,s);
return Pair2Int([r,s],p,q);
end;
#####################################################################

T:=0;
#####################################################################
Elts2Int:=function(x)
local pos;

pos:=Position(EltsE,x);
if IsPosInt(pos) then return pos;
else Print(x,"\n\n");
	Append(EltsE,[x]);
	if not IsInt(T) then
	Append(T!.elts,[x]); fi;
	return Length(EltsE); 
fi;
end;
#####################################################################

#####################################################################
Boundary:=function(k,jj)
local j, p,q,r,s,tmp, horizontal, vertical;

if k<1 then return []; fi;
j:=AbsInt(jj);
	#################IF BOUNDARY NOT ALREADY COMPUTED############
if IsInt(PseudoBoundary[k][j]) then
tmp:=Int2Vector(k,j);
p:=tmp[1]; q:=tmp[2]; r:=tmp[3]; s:=tmp[4];

horizontal:=ShallowCopy(BoundaryR(p,r));
Apply(horizontal,x->[x[1],Elts2Int(   ImagesRepresentative(GhomE,R!.elts[x[2]])   )  ]);
Apply(horizontal,x->[Vector2Int(p-1,q,x[1],s),x[2]]);

vertical:=ShallowCopy(BoundaryS(q,s));
Apply(vertical,x->[x[1],Elts2Int(   ImagesRepresentative(HhomE,S!.elts[x[2]])  )    ]);
Apply(vertical,x->[Vector2Int(p,q-1,r,x[1]),x[2]]);
if IsOddInt(p) then
vertical:=NegateWord(vertical);
fi;


PseudoBoundary[k][j]:= Concatenation(horizontal, vertical);
fi;
	################IF ENDS######################################

if SignInt(jj)=1 then
return PseudoBoundary[k][j];
else
return NegateWord(PseudoBoundary[k][j]);
fi;
end;
#####################################################################

#####################################################################
HorizontalBoundaryGen:=function(n,y)
local a,i, p,q,r,s, tmp,horizontal;

a:=AbsInt(y[1]);
tmp:=Int2Vector(n,a);

p:=tmp[1]; q:=tmp[2]; r:=tmp[3]; s:=tmp[4];

horizontal:=StructuralCopy(BoundaryR(p,r));

Apply(horizontal,x->[x[1],Elts2Int( EltsE[y[2]]*ImagesRepresentative(GhomE,R!.elts[x[2]]) )]);


Apply(horizontal,x->[Vector2Int(p-1,q,x[1],s),x[2]]);

return horizontal;

end;
#####################################################################

#####################################################################
HorizontalBoundaryWord:=function(n,w)
local x, bnd;

bnd:=[];
for x in w do
Append(bnd,HorizontalBoundaryGen(n,x));
od;
return bnd;

end;
#####################################################################

#####################################################################
HomotopyGradedGen:=function(g,p,q,r,s,bool)    #Assume EltsE[g] exists!
local aa,hty, hty1, Eg, Eg1, Eg2, g1, g2;	#bool=true for vertical homotopy

#This function seems to work! But I should really check the maths again!!

Eg:=EltsE[g];
Eg1:=ImagesRepresentative(EhomG,Eg);
				    
Eg2:=ImagesRepresentative(EhomH,Eg);
				    
g2:=Position(S!.elts,Eg2); 
g1:=Position(R!.elts,Eg1); 
Eg1:=ImagesRepresentative(GhomE,Eg1);
Eg2:=ImagesRepresentative(HhomE,Eg2);


hty:=HomotopyS(q,[s,g2]);
Apply(hty,x->[ Vector2Int(p,q+1,r,x[1]), ImagesRepresentative(HhomE,S!.elts[x[2]])]); 
Apply(hty,x->[ x[1], Elts2Int(Eg1*x[2])]);
if IsOddInt(p) then
hty:=NegateWord(hty); fi;

if (p=0 and q>0) or bool then return hty; fi;

if p>0 then
hty1:=HomotopyOfWord(p+q,StructuralCopy(HorizontalBoundaryWord(p+q+1,hty)),false);
Append(hty, NegateWord(hty1));
fi;

if q>0 then return hty; fi;


hty1:=HomotopyR(p,[r,g1]);
Apply(hty1,x->[ Vector2Int(p+1,q,x[1],s), ImagesRepresentative(GhomE,R!.elts[x[2]])]);
Apply(hty1,x->[ x[1], Elts2Int(x[2])]); #Here

Append(hty,hty1);

hty1:=HomotopyOfWord(p+q,StructuralCopy(HorizontalBoundaryWord(p+q+1,hty1)),true);

Append(hty,NegateWord(hty1));

hty1:=HomotopyOfWord(p+q,StructuralCopy(HorizontalBoundaryWord(p+q+1,hty1)),false);


Append(hty,hty1); 	#I think this perturbation term is always zero and
			#thus not necessary.

return hty;
end;
#####################################################################

#####################################################################
Homotopy:=function(n,x,bool)
local vec,a;


a:=AbsInt(x[1]);
vec:=Int2Vector(n,a);
if SignInt(x[1])=1 then 
return HomotopyGradedGen(x[2],vec[1],vec[2],vec[3],vec[4],bool);
else
return NegateWord(HomotopyGradedGen(x[2],vec[1],vec[2],vec[3],vec[4],bool));
fi;
end;
#####################################################################

#####################################################################
HomotopyOfWord:=function(n,w,bool)
local x, hty;

hty:=[];
for x in w do
Append(hty,Homotopy(n,x,bool));
od;
return hty;

end;
#####################################################################

HomotopyRec:=[];
for i in [1..Lngth] do
HomotopyRec[i]:=[];
for j in [1..Dimension(i-1)] do
HomotopyRec[i][j]:=[];
od;od;

#####################################################################
FinalHomotopy:=function(n,x)
local a;

a:=AbsInt(x[1]);
if not IsBound(HomotopyRec[n+1][a][x[2]]) then
HomotopyRec[n+1][a][x[2]]:=Homotopy(n,[a,x[2]],false);
fi;

if SignInt(x[1])=1 then
return StructuralCopy(HomotopyRec[n+1][a][x[2]]);
else
return NegateWord(StructuralCopy(HomotopyRec[n+1][a][x[2]]));
fi;

end;
#####################################################################

if HomotopyR=fail or HomotopyS=fail then
FinalHomotopy:=fail;
fi;

#for i in [1..Lngth] do           #I'm guessing there is no need for this
#for j in [1..Dimension(i)] do    #given that we are "lazy computing".
#g:=Boundary(i,j);
#od;
#od;


Boole:=false;
if EvaluateProperty(R,"reduced")=true
and  EvaluateProperty(S,"reduced")=true
then Boole:=true;
fi;


T:= rec(
            dimension:=Dimension,
	    boundary:=Boundary,
	    homotopy:=FinalHomotopy,
	    elts:=EltsE,
	    group:=E,
	    firstProjection:=EhomG,
	    secondProjection:=EhomH,
            Int2Vector:=Int2Vector,
            Vector2Int:=Vector2Int,
	    properties:=
	    [["type","resolution"],
	    ["length",Lngth],
	    ["characteristic",Charact],
	    ["reduced",true]]
	    );


return Objectify(HapResolution,T);
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(HAP_ResolutionAbelianGroupFromInvariants,
function(LL,n)
local L,G,R,S,pi,L1, L2, i, ResolutionAbelianGroup;
ResolutionAbelianGroup:=HAP_ResolutionAbelianGroupFromInvariants;

L:=SortedList(LL);
L1:=Filtered(L,x->x=0);
L2:=Filtered(L,x->x>0);
L2:=AbelianGroup(L2);
L2:=AbelianInvariants(L2);
L2:=AbelianInvariantsToTorsionCoefficients(L2);
L:=Concatenation(L1,L2);

L:=Reversed(L);
if not IsList(L) then 
Print("First argument must be a list of non-negative integers.\n");
return fail;
fi;

if Length(L)=0 then G:=AbelianPcpGroup([1]); return ResolutionFiniteGroup(G,n); fi;
if Length(L)=1 then 
if L[1]=0 then return ResolutionInfiniteCyclicGroup(n); fi;
pi:=[2..L[1]];Add(pi,1);pi:=PermList(pi);G:=Group(pi);
return ResolutionFiniteGroup(G,n);
fi;

R:=ResolutionAbelianGroup(L{[1]},n);
S:=ResolutionAbelianGroup(L{[2..Length(L)]},n);
return ResolutionDirectProductLazy(R,S);


end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ResolutionAbelianGroup,
function(A,n)
local mgens,L1,L2,R, S, T ;

if IsList(A) then return HAP_ResolutionAbelianGroupFromInvariants(A,n); fi;

if not IsAbelian(A) then Print("The first variable is not an abelian group\n");
return fail;
fi;


mgens:=MinimalGeneratingSet(A);
L1:=Filtered(mgens,x->Order(x)=infinity);
L2:=Filtered(mgens,x->Order(x)<infinity);
#if Length(L2)>0 then L2:=MinimalGeneratingSet(Group(L2)); fi;
if Length(L2)>0 then L2:=TorsionGeneratorsAbelianGroup(Group(L2)); fi;

mgens:=Concatenation(L1,L2);

if Length(mgens)=0 then return ResolutionFiniteGroup(Group(One(A)),n); fi;

if Length(mgens)=1 then
if Order(mgens[1])<infinity then 
     return ResolutionFiniteCyclicGroup(Group(mgens[1]),n); 
else
   return ResolutionInfiniteCyclicGroup(A,n);
fi;
fi;

R:=ResolutionAbelianGroup(Group(mgens{[1]}),n);
S:=ResolutionAbelianGroup(Group(mgens{[2..Length(mgens)]}),n);
SetParent(R!.group,A);
SetParent(S!.group,A);
T:=ResolutionDirectProductLazy(R,S);

return T;

end);
#####################################################################
#####################################################################

