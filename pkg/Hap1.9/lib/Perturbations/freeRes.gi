#(C) 2008 Graham Ellis


################################################################
InstallGlobalFunction(FreeResolution,
function(P,N)
local 
	Dimension, DimensionRecord, DimRecs, 
	Boundary,
	BoundaryP,
	Pair2Quad, Pair2QuadRec,
      	Quad2Pair,Quad2PairRec,
	HtpyGen, HtpyWord,
	StabGrps, 
	StabResls, 
	Action,
	EltsG, G, Mult,
	DelGen, DelWord, DelGenRec,
	PseudoBoundary,FinalBoundary,
	L,i,k,n,q,r,s,t ;

N:=Minimum(Length(P),N);
G:=P!.group;
EltsG:=P!.elts;
BoundaryP:=P!.boundary;

##############################################
if IsBound(P!.action) then 
Action:=P!.action;
else 
Action:=function(k,j,g) return 1; end;
fi;
##############################################

################################################################
Mult:=function(g,h);
return Position(EltsG,EltsG[g]*EltsG[h]);
end;
################################################################

StabGrps:= List([0..Length(P)],n->
           List([1..P!.dimension(n)], k->P!.stabilizer(n,k))); 

StabResls:=[];
i:=N;
##################################
for L in StabGrps do
Add(StabResls,List(L,	
g->ExtendScalars(ResolutionFiniteGroup(g,i),G,EltsG))
);
i:=Maximum(0,AbsInt(i-1));
od;
#################################

DimRecs:=List([0..N],i->[]);

###################################################################
Dimension:=function(k)
local dim,i,R;
dim:=0;
for i in [0..k] do
DimRecs[k+1][i+1]:=[];
for R in StabResls[i+1] do
dim:=dim+R!.dimension(k-i);
Add(DimRecs[k+1][i+1],dim);
od;
od;
return dim;
end;

DimensionRecord:=List([0..N],Dimension);

Dimension:=function(k);
return DimensionRecord[k+1];
end;
###################################################################

###################################################################
Quad2PairRec:=[];
for q in [0..N] do
Quad2PairRec[q+1]:=[];
for r in [1..Length(StabGrps[q+1])] do
Quad2PairRec[q+1][r]:=[];
for s in [0..N-q] do
Quad2PairRec[q+1][r][s+1]:=[];
od;od;od;
###################################################################


###################################################################
Pair2Quad:=function(k,n)
local qq,q,r,s,t;
#The n-th generator in degree k of our final resolution is actually the 
#t-th generator in degree s of the resolution of the r-th stabilizer group
#of the q-th chain module of the non-free resolution. We need the
#function f(k,n)=[q,r,s,t] .

for qq in [0..N] do
if n <= DimRecs[k+1][qq+1][Length(DimRecs[k+1][qq+1])] then q:=qq; break; fi;
od;

r:=PositionProperty(DimRecs[k+1][q+1],x->(n<=x));

s:=k-q;

if r-1>0 then
t:=n-DimRecs[k+1][q+1][r-1];
else
if q>=1 then
t:=n-DimRecs[k+1][q][Length(  DimRecs[k+1][q] )];;
else t:=n;
fi;
fi;
 
Quad2PairRec[q+1][r][s+1][t]:=[k,n];
return [q,r,s,t];

end;

Pair2QuadRec:=List([1..N+1],i->[]);
for k in [0..N] do
for n in [1..Dimension(k)] do
Pair2QuadRec[k+1][n]:=Pair2Quad(k,n);
od;
od;

##############
Pair2Quad:=function(k,n)
local a;
if n>0 then
return StructuralCopy(Pair2QuadRec[k+1][n]);
else
a:=StructuralCopy(Pair2QuadRec[k+1][-n]);
a[4]:=-a[4];
return a;
fi;
end;
##############

##############
Quad2Pair:=function(q,r,s,t)
local a,pr,pt;
if r>0 then pr:=r;pt:=t;
else
pr:=-r;pt:=-t;
fi;

if pt>0 then
return StructuralCopy(Quad2PairRec[q+1][pr][s+1][pt]);
else
a:=StructuralCopy(Quad2PairRec[q+1][pr][s+1][-pt]);
a[2]:=-a[2];
return a;
fi;
end;
##############
###################################################################

###################################################################
HtpyGen:=function(q,s,r,t,g)
local y,pr,pt;
#The applies the "vertical homotopy" to the free group generator [r,t,g]
#in "dimension" [q,s]. The output is an "r-word" in "dimension" [q,s+1].

if r>0 then pr:=r;pt:=t;
else
pr:=-r;pt:=-t;
fi;

y:=StabResls[q+1][pr]!.homotopy(s,[pt,g]);
return List(y,x->[pr,x[1],x[2]]);
end;
###################################################################

###################################################################
HtpyWord:=function(q,s,w)
local h,x;
#This applies the "vertical homotopy" to the r-word w in "dimension" 
#[q,s]. The output is an r-word in "dimension" [q,s+1].

h:=[];
for x in w do
Append(h,HtpyGen(q,s,x[1],x[2],x[3]));
od;

return h;
end;
###################################################################

DelGenRec:=[];
for k in [1..N+1] do
DelGenRec[k]:=[];
for q in [1..N+1] do
DelGenRec[k][q]:=[];
for s in [1..N+1] do
DelGenRec[k][q][s]:=[];
for r in [1..P!.dimension(q-1)] do
DelGenRec[k][q][s][r]:=[];
od;
od;
od;
od;
###################################################################
DelGen:=function(k,q,s,r,t)
local y,pr,pt,i;
#For k=0,1,2 ... this is the equivariant homomorphism
#Del_k:A_{q,s} ---> A_{q-k,s+k-1} applied to a free r-generator [r,t]
#in dimension [q,s].

if r>0 then pr:=r;pt:=t;
else
pr:=-r;pt:=-t;
fi;

##############
if IsBound(DelGenRec[k+1][q+1][s+1][pr][AbsInt(pt)]) then

if pt>0 then return DelGenRec[k+1][q+1][s+1][pr][pt];
else
return List(DelGenRec[k+1][q+1][s+1][pr][-pt], a->[a[1],-a[2],a[3]]);
fi;

fi;
##############

if k=0 then
if s=0 then return [];
else
y:=List(StabResls[q+1][pr]!.boundary(s,pt),x->[Action(q,r,x[2])*x[1],x[2]]);
if pt>0 then
DelGenRec[k+1][q+1][s+1][pr][pt]:= List(y,x->[pr,x[1],x[2]]);
return DelGenRec[k+1][q+1][s+1][pr][pt];
else DelGenRec[k+1][q+1][s+1][pr][-pt]:=List(y,x->[pr,-x[1],x[2]]);
return List(y,x->[pr,x[1],x[2]]);
fi;
fi;
fi;

if k=1 then
if s=0 then  
if q=0 then return [];
fi;
y:=BoundaryP(q,pr);
if pt>0 then 
DelGenRec[k+1][q+1][s+1][pr][pt]:= List(y,x->[x[1],1,x[2]]);
return DelGenRec[k+1][q+1][s+1][pr][pt];
else
DelGenRec[k+1][q+1][s+1][pr][-pt]:= List(y,x->[x[1],1,x[2]]);
return List(y,x->[x[1],-1,x[2]]);
fi;

else

if pt>0 then
DelGenRec[k+1][q+1][s+1][pr][pt]:=
HtpyWord(q-1,s-1,DelWord(1,q,s-1,DelGen(0,q,s,pr,-pt))) ;
return DelGenRec[k+1][q+1][s+1][pr][pt];
else
DelGenRec[k+1][q+1][s+1][pr][-pt]:=
HtpyWord(q-1,s-1,DelWord(1,q,s-1,DelGen(0,q,s,pr,pt))) ;
return 
List(DelGenRec[k+1][q+1][s+1][pr][-pt], a->[a[1],-a[2],a[3]]);
fi;

fi;
fi;

y:=[];
for i in [1..k] do
Append(y,
HtpyWord(q-k,s+k-2,DelWord(i,q-k+i,s+k-i-1,DelGen(k-i,q,s,pr,-pt)))
);
od;


if pt>0 then
DelGenRec[k+1][q+1][s+1][pr][pt]:=y;
else
DelGenRec[k+1][q+1][s+1][pr][-pt]:=List(y,a->[a[1],-a[2],a[3]]);
fi;

return y;
end;
###################################################################

###################################################################
DelWord:=function(k,q,s,w)
local y,x;
#For k=0,1,2 ... this is the equivariant homomorphism
#Del_k:A_{q,s} ---> A_{q-k,s+k-1} applied to an r-word [[r,t,g],...]
#in dimension [q,s].

y:=[];
for x in w do
Append(y,List(DelGen(k,q,s,x[1],x[2]),
a->[a[1],Action(q-k,a[1],x[3])*a[2],Mult(x[3],a[3])])); #THIS CAN'T BE RIGHT!
od;
return y;

end;
###################################################################

###################################################################
Boundary:=function(k,n)
local q,s,r,t,x,y,z,i;

y:=Pair2Quad(k,n); q:=y[1];s:=y[3];r:=y[2];t:=y[4];

y:=[];

for i in [0..k] do
if q>=i then
z:=DelGen(i,q,s,r,t);
Append(y,
List(z,x->[Quad2Pair(q-i,x[1],s+i-1,x[2])[2],x[3]])  );
else break;
fi;
od;

return AlgebraicReduction(y);
end;
###################################################################

PseudoBoundary:=[];
for n in [1..N+1] do
PseudoBoundary[n]:=[];
od;
#######################################
FinalBoundary:=function(n,k)
local  pk;

pk:=AbsInt(k);
if  not IsBound(PseudoBoundary[n+1][pk]) then
PseudoBoundary[n+1][pk]:= Boundary(n,pk);
fi;

if k>0 then return PseudoBoundary[n+1][k];
else return NegateWord(PseudoBoundary[n+1][pk]); fi;
end;
#######################################



return Objectify(HapResolution,
                rec(
                dimension:=Dimension,
                boundary:=FinalBoundary,
                homotopy:=fail,
                elts:=P!.elts,
                group:=P!.group,
                properties:=
                   [["length",N],
                    ["reduced",true],
                    ["type","resolution"],
                    ["characteristic",EvaluateProperty(P,"characteristic")]  ]));


end);
################################################################
################################################################



################################################################
################################################################
InstallGlobalFunction(ExtendScalars,
function(arg)
# Here H=S!.group is a subgroup of G, and EltsG is a list of the 
# elements of G.
local 
	S,G,EltsG, 	
	R,T,
	H,EltsH,
	HhomG,GmapTH,THmapG,
	BoundaryS,Boundary,
	HomotopyS,Homotopy,
	PseudoBoundary,n,k,PseudoHomotopy,FinalHomotopy;

S:=arg[1];
G:=arg[2];
EltsG:=arg[3];
H:=S!.group;
EltsH:=S!.elts;

BoundaryS:=S!.boundary;
HomotopyS:=S!.homotopy;

#######################################
HhomG:=function(i);
return Position(EltsG,EltsH[i]);
end;
#######################################

PseudoBoundary:=[];
for n in [1..Length(S)+1] do
PseudoBoundary[n]:=[];
od;
#######################################
Boundary:=function(n,k)
local  pk;

pk:=AbsInt(k);
if  not IsBound(PseudoBoundary[n+1][pk]) then
PseudoBoundary[n+1][pk]:= List(BoundaryS(n,pk),x->[x[1],HhomG(x[2])]);
fi;

if k>0 then return PseudoBoundary[n+1][k];
else return NegateWord(PseudoBoundary[n+1][pk]); fi;
end;
#######################################

#######################################
GmapTH:=function(g)    #ht=g^-1 ==> g=t^-1 h^-1
local t,h,gg;
gg:=EltsG[g]^-1;
t:=CanonicalRightCosetElement(H,gg)^-1;
h:=(gg*t)^-1;

return [Position(EltsG,t),Position(EltsH,h)];
end;
#######################################

#######################################
THmapG:=function(t,h);
return Position(EltsG,EltsG[t]*EltsG[HhomG(h)]);
end;
#######################################

#######################################
Homotopy:=function(n,p)
local ht,h,t,g,k,htpy;

k:=p[1];
g:=p[2];
ht:=GmapTH(g);
h:=ht[2];t:=ht[1];
htpy:=StructuralCopy(HomotopyS(n,[k,h]));
Apply( htpy,x->[x[1],THmapG(t,x[2])] );
return htpy;
end;
#######################################

PseudoHomotopy:=[];
for n in [1..Length(S)] do
PseudoHomotopy[n]:=[];
for k in [1..S!.dimension(n-1)] do
PseudoHomotopy[n][k]:=[];
od;
od;
#######################################
FinalHomotopy:=function(n,p)
local t,g,pt;


t:=p[1];g:=p[2];pt:=AbsInt(t);

if  not IsBound(PseudoHomotopy[n+1][pt][g]) then
PseudoHomotopy[n+1][pt][g]:= Homotopy(n,[pt,g]);
fi;

if t>0 then return PseudoHomotopy[n+1][pt][g];
else return NegateWord(PseudoHomotopy[n+1][pt][g]); fi;
end;
#######################################

return Objectify(HapResolution,
                rec(
                dimension:=S!.dimension,
                boundary:=Boundary,
                homotopy:=FinalHomotopy,
                elts:=EltsG,
                group:=G,
                properties:=S!.properties
                ));

end);
################################################################
################################################################


################################################################
################################################################
InstallGlobalFunction(InduceScalars,
function(S,hom)
local
	G,Q,N,R,StabilizerSubgroup,
	Boundary,BoundaryS,
	QmapG,  
	EltsG, EltsQ,
	PseudoBoundary, FinalBoundary,n;

G:=Source(hom);
N:=Kernel(hom);
EltsG:=Elements(G);
EltsQ:=S!.elts;
BoundaryS:=S!.boundary;

###############################################################
StabilizerSubgroup:=function(k,n);
return N;
end;
################################################################

#################################################################
QmapG:=function(q);
return Position(EltsG,PreImagesRepresentative(hom,EltsQ[q]));
end;
#################################################################


#################################################################
Boundary:=function(k,n);
return List(BoundaryS(k,n),x->[x[1],QmapG(x[2])]);
end;
#################################################################

PseudoBoundary:=[];
for n in [1..Length(S)+1] do
PseudoBoundary[n]:=[];
od;
#######################################
FinalBoundary:=function(n,k)
local  pk;

pk:=AbsInt(k);
if  not IsBound(PseudoBoundary[n+1][pk]) then
PseudoBoundary[n+1][pk]:= Boundary(n,pk);
fi;

if k>0 then return PseudoBoundary[n+1][k];
else return NegateWord(PseudoBoundary[n+1][pk]); fi;
end;
#######################################

return Objectify(HapNonFreeResolution,
           rec(
            dimension:=S!.dimension,
            boundary:=FinalBoundary,
            homotopy:=fail,
            elts:=EltsG,
            group:=G,
            stabilizer:=StabilizerSubgroup,
            properties:=
             [["type","resolution"],
              ["length",EvaluateProperty(S,"length")],
              ["characteristic", EvaluateProperty(S,"characteristic")] ]));

end);
################################################################
################################################################

################################################################
################################################################
InstallGlobalFunction(CoxeterComplex,
function(arg)
local 
	D,N,R, A, W, WP, EltsWP, WPev,
	AhomW, WhomWP, AhomWP, 
	ResGens, 
	StabilizerSubgroup,
	Action,
	i, n, k,x;

D:=arg[1];
###########################
if not CoxeterDiagramIsSpherical(D) then
Print("This function is only implemented for finite Coxeter groups.\n");
return fail;
fi;
###########################



if Length(arg)>1 then N:=arg[2];
else N:=1000; fi;

R:=ResolutionArtinGroup(D,N); #I guess no one will ever try 
				 #more than 1000 generators!
A:=R!.group;
N:=Minimum(N,Length(GeneratorsOfGroup(A)));

for n in [1..N] do
for k in [1..R!.dimension(n)] do
i:=R!.boundary(n,k);
od; od;

ResGens:=R!.resGens;
W:=CoxeterDiagramFpCoxeterGroup(D);
W:=W[1]/W[2];
AhomW:=GroupHomomorphismByImagesNC(A,W,GeneratorsOfGroup(A),GeneratorsOfGroup(W));
WhomWP:=IsomorphismPermGroup(W);
WP:=Image(WhomWP);
WPev:=EvenSubgroup(WP);
AhomWP:=GroupHomomorphismByFunction(A,WP,x->Image(WhomWP,Image(AhomW,x)));

EltsWP:=List(R!.elts,x->Image(AhomWP,x));
EltsWP:=Concatenation(
EltsWP,
Filtered(Elements(WP),x->not x in EltsWP));


###############################################################
StabilizerSubgroup:=function(n,k)
local G;

G:=List(ResGens[n+1][k], x->GeneratorsOfGroup(WP)[x]);
if Length(G)>0 then return Group(G);fi;
return Group(());
end;
###############################################################

###############################################################
# This describes how the group WP acts on the orientation.
Action:=function(k,j,g);
if 
EltsWP[g] in WPev then return 1;
else return -1; fi;
end;
###############################################################

return Objectify(HapNonFreeResolution,
           rec(
            dimension:=R!.dimension,
            boundary:=R!.boundary,
            homotopy:=fail,
            elts:=EltsWP,
            group:=WP,
            stabilizer:=StabilizerSubgroup,
	    action:=Action,
            properties:=
             [["type","resolution"],
              ["length",N],
              ["characteristic", 0] ]));

end);
################################################################
################################################################

################################################################
################################################################
InstallGlobalFunction(ResolutionCoxeterGroup,
function(D,n)
local P,R;

###########################
if not CoxeterDiagramIsSpherical(D) then
Print("This function is only implemented for finite Coxeter groups.\n");
return fail;
fi;
###########################

P:=CoxeterComplex(D,n);
return FreeResolution(P,n);

end);
################################################################
################################################################
