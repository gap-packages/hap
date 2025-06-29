#(C) 2008 Graham Ellis
#RT:=0;
################################################################
InstallGlobalFunction(FreeGResolution,
function(arg)
local 
	P,N,prime,
	Dimension, DimensionRecord, DimRecs, FiltDimRecs,
        BinGp,
	Boundary,
	BoundaryP,
	Pair2Quad, Pair2QuadRec,
      	Quad2Pair,Quad2PairRec,
	HtpyGen, HtpyWord,
	StabGrps, 
	StabResls, 
        ResolutionFG,
	Action,
        AlgRed,  
	EltsG, G, Mult, MultRecord,
	DelGen, DelWord, DelGenRec,
	PseudoBoundary,FinalBoundary,
        FilteredLength, FilteredDimension, FilteredDimensionRecord,
	L,i,k,n,q,r,s,t,bool;

SetInfoLevel(InfoWarning,0);

P:=arg[1];
N:=arg[2];
if Length(arg)>2 then prime:=Gcd(arg[3],EvaluateProperty(P,"characteristic"));
else prime:=EvaluateProperty(P,"characteristic"); fi;

N:=Minimum(EvaluateProperty(P,"length"),N);
G:=P!.group;
bool:=not IsComponentObjectRep(One(G));
bool:=IsHapSL2Subgroup(G) or IsHapSL2OSubgroup(G) or IsBound(G!.bianchiInteger);
EltsG:=P!.elts;
BoundaryP:=P!.boundary;

BinGp:=ContractibleGcomplex("SL(2,O-2)");
BinGp:=BinGp!.stabilizer(0,4);;
BinGp:=Image(RegularActionHomomorphism(BinGp));
#BinGp:=Group(ReduceGenerators(GeneratorsOfGroup(BinGp),BinGp));

#############################
ResolutionFG:=function(G,n)
local x, tmp, iso,iso1,iso2,iso3,res,Q, fn;

##Added April 2017
if prime>0 and Order(G)>1 and IsPGroup(G) then return ResolutionPrimePowerGroup(G,n); fi;
##
##Added Jan 2012
if IsBound(P!.resolutions) and HasName(G) then
x:=Position(P!.resolutions[2], Name(G));
if not x=fail then return P!.resolutions[1][x]; fi;
fi;
##


###
if IsBianchiAbelianGroup(G) and Order(G)=infinity then 
#if IsAbelian(G) and Order(G)=infinity then
res:=ResolutionAbelianBianchiSubgroup(G,n);
return res;
fi;
if IsAbelian(G) and not bool then   #NOT SURE WHY I NEED TO DO THIS
#if false then
res:=ResolutionAbelianGroup_alt(G,n);

return res;
fi;
###

iso:=RegularActionHomomorphism(G);
Q:=Image(iso);

if Order(Image(iso))=24 then
if IdGroup(Image(iso))=[24,3] then 
iso1:=IsomorphismGroups(Q,BinGp); 
res:=ResolutionFiniteGroup(BinGp,n);
res!.group:=G;
res!.elts:=List(res!.elts,x->
PreImagesRepresentative(iso,PreImagesRepresentative(iso1,x)));
return res;
fi;
fi;

#res:=ResolutionFiniteGroup(Q,n);
res:=ResolutionGenericGroup(Q,n);
res!.group:=G;
res!.elts:=List(res!.elts,x->PreImagesRepresentative(iso,x));
return res;
###

end;
#############################

if prime>0 then 
##############################################
AlgRed:= function(ww)
local w,x,v,pos,u;

w:=StructuralCopy(ww);

v:=Collected(w);
for x in v do
if x[1][1]<0 then x[1][1]:=-x[1][1]; x[2]:=-x[2] mod prime; fi;
if x[1][2]<0 then x[1][2]:=-x[1][2]; x[2]:=-x[2] mod prime; fi;
x[2]:=x[2] mod prime;
od;

u:=[];
for x in v do
Append(u,MultiplyWord(x[2],[x[1]]));
od;

v:=Collected(u);
for x in v do
x[2]:=x[2] mod prime;
od;

u:=[];
for x in v do
Append(u,MultiplyWord(x[2],[x[1]]));
od;

return u;

end;
##############################################

else

##############################################
AlgRed:= function(ww)
local x,i,v,k,u,w, vv, vvv, vvvv, pos, pos2;

w:=ww;
        for x in w do
        if x[2]<0 then x[1]:=-x[1];x[2]:=-x[2];fi;
        od;

#v:=Filtered(w,x->x[1]>0);
#        for x in w do
#        if x[1]<0 then
##########################
#        k:=Position(v,[-x[1],x[2],x[3]]);  
#        if (k=fail) then Add(v,x); 
#        else
#        Unbind(v[k]);
#        fi;
##########################
#        fi;
#        od;
#        v:=Filtered(v,x->IsBound(x));
#        return v;
         v:=Collected(w);
         vv:=List(v,x->x[1]);
         vvv:=List(vv,x->[AbsInt(x[1]),x[2],x[3]]);
         vvv:=SSortedList(vvv);
         vvvv:=[];
         for x in vvv do
             pos:=Position(vv,x);
             pos2:=Position(vv,[-x[1],x[2],x[3]]);
             k:=0;
             if not pos=fail then k:=k+v[pos][2]; fi;
             if not pos2=fail then k:=k-v[pos2][2]; fi;
             if not k=0 then Append(vvvv,MultiplyWord(k,[ [x[1],x[2], x[3]] ]  )); fi;

         od;
         return vvvv;

end;
##############################################
fi;



##############################################
if IsBound(P!.action) and not prime=2 then 
Action:=P!.action;
else 
Action:=function(k,j,g) return 1; end;
fi;
##############################################

MultRecord:=[];
################################################################
Mult:=function(g,h)
local pos;
if not IsBound(MultRecord[g]) then MultRecord[g]:=[]; fi;
if not IsBound(MultRecord[g][h]) then
    pos:= Position(EltsG,EltsG[g]*EltsG[h]);
    if pos=fail then Add(EltsG,EltsG[g]*EltsG[h]);  
    MultRecord[g][h]:= Length(EltsG);
    else MultRecord[g][h]:= pos; 
    fi;
fi;
return MultRecord[g][h];
end;
################################################################

StabGrps:= List([0..Length(P)],n->
           List([1..P!.dimension(n)], k->P!.stabilizer(n,k))); 

StabResls:=[];
i:=N;
#if true then
#################################
for L in StabGrps do
Add(StabResls,List(L,	
g->ExtendScalars(ResolutionFG(g,i),G,EltsG))
);
i:=Maximum(0,AbsInt(i-1));
od;
#################################
#else
##################################
#for L in StabGrps do
#Add(StabResls,List(L,
##g->ExtendScalars(ResolutionFiniteGroup(g,i,false,prime),G,EltsG)));
#g->ExtendScalars(ResolutionGenericGroup(g,i,false,prime),G,EltsG)));
#i:=Maximum(0,AbsInt(i-1));
#od;
##################################
#fi;

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
#This applies the "vertical homotopy" to the free group generator [r,t,g]
#in "dimension" [q,s]. The output is an "r-word" in "dimension" [q,s+1].

if r>0 then pr:=r;pt:=t;
else
pr:=-r;pt:=-t;
fi;

y:=StructuralCopy(StabResls[q+1][pr]!.homotopy(s,[pt,g]));
Apply(y,x->[pr,x[1],x[2]]);
return y;
end;
###################################################################

###################################################################
HtpyWord:=function(q,s,w)
local h,z,x,y;
#This applies the "vertical homotopy" to the r-word w in "dimension" 
#[q,s]. The output is an r-word in "dimension" [q,s+1].

h:=[];
for y in w do
x:=[Action(q,y[1],y[3])*y[1],y[2],y[3]];
z:=HtpyGen(q,s,x[1],x[2],x[3]);
z:=List(z,a->[Action(q,a[1],a[3])*a[1],a[2],a[3]]);
Append(h,z);
od;

return AlgRed(h);
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
DelGenRec[k+1][q+1][s+1][pr][pt]:= AlgRed(List(y,x->[pr,x[1],x[2]]));
return DelGenRec[k+1][q+1][s+1][pr][pt];
else DelGenRec[k+1][q+1][s+1][pr][-pt]:=AlgRed(List(y,x->[pr,-x[1],x[2]]));
return AlgRed(List(y,x->[pr,x[1],x[2]]));   
fi;
fi;
fi;

if k=1 then
if s=0 then  
if q=0 then return [];
fi;
y:=BoundaryP(q,pr);
if pt>0 then 
DelGenRec[k+1][q+1][s+1][pr][pt]:= AlgRed(List(y,x->[x[1],1,x[2]]));
return DelGenRec[k+1][q+1][s+1][pr][pt];
else
DelGenRec[k+1][q+1][s+1][pr][-pt]:= AlgRed(List(y,x->[x[1],1,x[2]]));
return List(y,x->[x[1],-1,x[2]]);
fi;

else

if pt>0 then
DelGenRec[k+1][q+1][s+1][pr][pt]:=
AlgRed(HtpyWord(q-1,s-1,DelWord(1,q,s-1,DelGen(0,q,s,pr,-pt)))) ;
return DelGenRec[k+1][q+1][s+1][pr][pt];
else
DelGenRec[k+1][q+1][s+1][pr][-pt]:=
AlgRed(HtpyWord(q-1,s-1,DelWord(1,q,s-1,DelGen(0,q,s,pr,pt)))) ;

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
y:=AlgRed(y);

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
a->[a[1],a[2],Mult(x[3],a[3])]));
od;

return y;  #Added Jan 2013. Speeds up the calculation in some(!!) examples.
return AlgRed(y);

end;
###################################################################

###################################################################
Boundary:=function(k,n)
local q,s,r,t,x,y,z,i;
y:=Pair2Quad(k,n); q:=y[1];s:=y[3];r:=y[2];t:=y[4];

y:=[];

for i in [0..k] do
#for i in [0..1] do
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
if k>0 then  return PseudoBoundary[n+1][k];
else  return NegateWord(PseudoBoundary[n+1][pk]); fi;
end;
#######################################



################spectral sequence requirements##################


FiltDimRecs:=[];
for k in [0..N] do
FiltDimRecs[k+1]:=[];
for i in [1..Dimension(k)] do
FiltDimRecs[k+1][i]:=Pair2Quad(k,i)[1];
od;
od;

FilteredLength:=Maximum(Flat(FiltDimRecs));

##################################################
FilteredDimension:=function(r,k);

return Length(Filtered(FiltDimRecs[k+1],x->x<=r));

end;
##################################################

SetInfoLevel(InfoWarning,1);

return Objectify(HapResolution,
                rec(
                dimension:=Dimension,
                filteredDimension:=FilteredDimension,
                boundary:=FinalBoundary,
                homotopy:=fail,
                elts:=P!.elts,
                group:=P!.group,
                pseudoBoundary:=PseudoBoundary,
                pair2Quad:=Pair2Quad,
                quad2Pair:=Quad2Pair,
                properties:=
                   [["length",N],
                    ["filtration_length",FilteredLength],
                    ["initial_inclusion",true],
                    ["reduced",true],
                    ["type","resolution"],
                    ["characteristic",prime]  ]));


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
	HhomGrec,GmapTHrec,THmapGrec,
	BoundaryS,Boundary,
	HomotopyS,Homotopy,
	PseudoBoundary,n,k,PseudoHomotopy,FinalHomotopy,
	PosG;

S:=arg[1];
G:=arg[2];
EltsG:=arg[3];

PosG:=Position;

H:=S!.group;
EltsH:=S!.elts;

BoundaryS:=S!.boundary;
HomotopyS:=S!.homotopy;

HhomGrec:=[];
#######################################
HhomG:=function(i)
local pos;

if IsBound(HhomGrec[i]) then return HhomGrec[i]; fi;

pos:= PosG(EltsG,EltsH[i]);
if pos=fail then Add(EltsG,EltsH[i]);   HhomGrec[i]:=Length(EltsG); 
else HhomGrec[i]:= pos; fi;

return HhomGrec[i];
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

GmapTHrec:=[];
#######################################
GmapTH:=function(g)    #ht=g^-1 ==> g=t^-1 h^-1
local t,h,gg,pos1,pos2;

if IsBound(GmapTHrec[g]) then return GmapTHrec[g]; fi;

gg:=EltsG[g]^-1;
#t:=CanonicalRightCosetElement(H,gg)^-1;
t:=CanonicalRightCountableCosetElement(H,gg)^-1;
#t:=CanonicalRightCosetElement(H,gg)^-1;

h:=(gg*t)^-1;

pos1:=PosG(EltsG,t);
if pos1=fail then Add(EltsG,t); pos1:=Length(EltsG);fi;
pos2:=Position(EltsH,h);
if pos2=fail then Add(EltsH,h); pos2:=Length(EltsH);fi;

GmapTHrec[g]:= [pos1,pos2];

return GmapTHrec[g];
end;
#######################################

THmapGrec:=[];
#######################################
THmapG:=function(t,h)
local pos, g;

if not IsBound(THmapGrec[t]) then THmapGrec[t]:=[]; fi;

if IsBound( THmapGrec[t][h] ) then return THmapGrec[t][h]; fi;

g:=EltsG[t]*EltsG[HhomG(h)];

pos:= PosG(EltsG,g);


if pos=fail then Add(EltsG, g);

 THmapGrec[t][h]:= Length(EltsG);
else  THmapGrec[t][h]:= pos; fi;

return  THmapGrec[t][h];
end;
#######################################

#######################################
Homotopy:=function(n,p)
local ht,h,t,g,k,htpy;

k:=p[1];
g:=p[2];
ht:=GmapTH(g);
h:=ht[2];t:=ht[1];
htpy:=(HomotopyS(n,[k,h]));
return List( htpy,x->[x[1],THmapG(t,x[2])] );
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
QmapG:=function(q)
local pos;
pos:= Position(EltsG,PreImagesRepresentative(hom,EltsQ[q]));
if pos = fail then Add(EltsG,PreImagesRepresentative(hom,EltsQ[q]));
return Length(EltsG);
else return pos; fi;
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
        Dimension,
	StabilizerSubgroup,
	Action,
        StandardWord, 
	i, n, k,x;

D:=arg[1];
if Length(arg)>1 then N:=arg[2];
else N:=1000; fi;

###########################
if not CoxeterDiagramIsSpherical(D) then
#Print("This function is only implemented for finite Coxeter groups.\n");
#return fail;
return CoxeterComplex_alt(D,N);
fi;
###########################




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


#####################################################################
Dimension:=function(n);

if n=0 then return 1; fi;
if n>Length(R) then return 0;
else return R!.dimension(n); fi;

end;
#####################################################################

#####################################################################
StandardWord:=function(k,bnd)
local w;
w:=
List(bnd,x->[x[1],
Position(EltsWP,  CanonicalRightCosetElement(StabilizerSubgroup(k,AbsInt(x[1])), EltsWP[x[2]]^-1 )^-1)
]);
return AlgebraicReduction(w);
end;
#####################################################################



return Objectify(HapNonFreeResolution,
           rec(
            dimension:=Dimension,
            boundary:=R!.boundary,
            homotopy:=fail,
            elts:=EltsWP,
            group:=WP,
            stabilizer:=StabilizerSubgroup,
	    action:=Action,
            standardWord:=StandardWord,
            properties:=
             [["type","resolution"],
              ["length",N],
              ["characteristic", 0] ]));

end);
################################################################
################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(CoxeterComplex_alt,
function(D,K)
local
        Dimension,
        Boundary,
        Contraction,    #not yet used
        EltsG,
        Vertices,
        G, gensG,
        W, Wgens, GhomW,
        ResGens,
        BoundaryCoeff,
        PseudoBoundary,
        BoundaryRecord,
        Action, StabilizerSubgroup,
        m, n, S, SD, LN;

Vertices:=CoxeterDiagramVertices(D);
G:=CoxeterDiagramMatCoxeterGroup(D);
gensG:=GeneratorsOfGroup(G);
LN:=Length(gensG[1]);
EltsG:=[];

ResGens:=[];
ResGens[1]:=[[]];
for n in [1..K] do
ResGens[n+1]:=[];
for S in Combinations(Vertices,n) do
SD:=CoxeterSubDiagram(D,S);
if CoxeterDiagramIsSpherical(SD) then AddSet(ResGens[n+1],S); fi;
od;
od;

#####################################################################
Dimension:=function(n);

if n=0 then return 1; fi;
if n>=Length(ResGens) then return 0;
else return Length(ResGens[n+1]); fi;

end;
#####################################################################

BoundaryRecord:=[];
for n in [1..K] do
BoundaryRecord[n]:=[];
for m in [1..Dimension(n)] do
BoundaryRecord[n][m]:=true;
od;
od;


#####################################################################
BoundaryCoeff:=function(S,T)    #S is a set of vertices generating a
                                #finite Coxeter group WS. T is a
                                #subset of S, and WT is the corresponding
                                #subgroup of WS.
local   SD, WS, gensWS,
        WT, gensWT,
        Trans,
        WShomG, Ggens,
        x,y;

SD:=CoxeterSubDiagram(D,S);
WS:=CoxeterDiagramFpCoxeterGroup(SD);
WS:=WS[1]/WS[2];
Ggens:=List(S,x->gensG[Position(Vertices,x)]);
gensWS:=GeneratorsOfGroup(WS);
WShomG:=GroupHomomorphismByImagesNC(WS,G,gensWS,Ggens);
gensWT:=List(T,x->gensWS[Position(S,x)]);
if Length(T)>0 then WT:=Group(gensWT);
else WT:=Group(Identity(WS)); fi;

Trans:=List(Elements(RightTransversal(WS,WT)),x->x^-1);

for x in Trans do
y:=Image(WShomG,x);
if not y in EltsG then Append(EltsG,[y]); fi;
od;

return List(Trans,x->Image(WShomG,x));
end;
#####################################################################

#####################################################################
PseudoBoundary:=function(S)     #S is a subset of vertices with finite
                                #Coxeter group WS.
local T, bndry, a;

bndry:=[];
for T in Combinations(S,Length(S)-1) do
a:=Difference(S,T)[1];
Append(bndry,[  [T,BoundaryCoeff(S,T),Position(S,a)]  ]);
od;

return bndry;
end;
#####################################################################

#####################################################################
Boundary:=function(n,kk)
local B, B1, FreeGWord, x, y, k;

#n:=AbsoluteValue(m);
if n<1 then return 0; fi;

k:=AbsoluteValue(kk);

if not BoundaryRecord[n][k]=true then
if kk>0 then return BoundaryRecord[n][k];
else return NegateWord(BoundaryRecord[n][k]);fi;
fi;

B:=PseudoBoundary(ResGens[n+1][k]);
#B1:=List(B,x->[Position(ResGens[n],x[1]),
#        List(x[2],y->(-1)^(Length(y)+x[3])*Position(EltsG,y))  ]);
B1:=List(B,x->[Position(ResGens[n],x[1]),
        List(x[2],y->Determinant(y)*(-1)^(x[3])*Position(EltsG,y))  ]);

FreeGWord:=[];
for x in B1 do
for y in x[2] do
Append(FreeGWord,[ [SignInt(y)*x[1],AbsoluteValue(y)] ]);
od;
od;
BoundaryRecord[n][k]:=FreeGWord;
if kk>0 then return FreeGWord;
else return NegateWord(FreeGWord); fi;
end;
#####################################################################

###############################################################
StabilizerSubgroup:=function(n,k)
local G,S;

S:=ResGens[n+1][k];
S:=List(S,i->gensG[i]  ) ;
if Length(S)=0 then return Group(IdentityMat(LN)); fi;
return Group(S);

end;
###############################################################

###############################################################
# This describes how the group WP acts on the orientation.
Action:=function(n,k,g);
if n=0 then return 1; fi;
return Determinant( EltsG[g]);
end;
###############################################################



return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:=EltsG,
            group:=G,
            stabilizer:=StabilizerSubgroup,
            action:=Action,
            properties:=
            [["length",n],
             ["characteristic",0],
             ["type","resolution"],
             ["reduced",true]]  ));

end);
#####################################################################
#####################################################################

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
return FreeGResolution(P,n);

end);
################################################################
################################################################



################################################################
################################################################
InstallGlobalFunction(TwistedResolution,
function(R,Action)
local N, n,k,g,Boundary, BoundaryRec,Homotopy,HomotopyRec;

N:=Length(R);

##################################################
Boundary:=function(n,k)
local bnd;

bnd:=StructuralCopy(R!.boundary(n,k));
Apply(bnd,x->[Action(x[2])*x[1],x[2]]);

return bnd;
end;
##################################################
BoundaryRec:=[];
for n in [1..N] do
BoundaryRec[n]:=[];
for k in [1..R!.dimension(n)] do
BoundaryRec[n][k]:=StructuralCopy(Boundary(n,k));
od;
od;
##################################################
Boundary:=function(n,k);
if k> 0 then
return BoundaryRec[n][k];
else return
NegateWord(BoundaryRec[n][AbsInt(k)]);
fi;
end;
##################################################


##################################################
Homotopy:=function(n,x)
local htpy;

htpy:=[Action(x[2])*x[1],x[2]];
htpy:=StructuralCopy(R!.homotopy(n,htpy));
Apply(htpy,y->[Action(y[2])*y[1],y[2]]);
return htpy;

end;
##################################################

HomotopyRec:=[];
for n in [0..N-2] do
HomotopyRec[n+1]:=[];
for k in [1..R!.dimension(n)] do
HomotopyRec[n+1][k]:=[];
for g in [1..Length(R!.elts)] do
HomotopyRec[n+1][k][g]:=Homotopy(n,[k,g]);
od;od;od;
##################################################
Homotopy:=function(n,x);

if x[1]>0 then return
HomotopyRec[n+1][x[1]][x[2]];
else return
NegateWord(HomotopyRec[n+1][AbsInt(x[1])][x[2]]);
fi;

end;
##################################################

return Objectify(HapNonFreeResolution,
           rec(
            dimension:=R!.dimension,
            boundary:=Boundary,
            homotopy:=Homotopy,
            elts:=R!.elts,
            group:=R!.group,
            properties:=
             [["type","resolution"],
              ["length",N],
              ["characteristic", 0] ]));


end);
################################################################
################################################################





