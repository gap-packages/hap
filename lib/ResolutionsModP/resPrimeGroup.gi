#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(ResolutionPrimePowerGroup,
function(arg)
local
	G,n,
	eltsG,
	gensG,
	Dimension,
	Boundary,
	PartialHomotopy,
	Homotopy,
	PseudoBoundary,
	PseudoBoundaryAsVec,
	WordToVectorList,
	VectorListToWord,
	prime, pp,
	BoundaryMatrix,
	MT,
	GactMat,
	ZGbasisOfKernel,
	one,
	InverseFlat,
	ComplementaryBasis,
	zero,
	pcgens,
	BoundaryMatrices,
	BoundaryMatrices2,
	Echelonize,
	EchelonMatrices,
	SolutionMatBoundaryMatrices,
	Toggle,
	SMBM,
	SMBMhom,
	ImageFBasis,
	HomotopyRec,
	InitHomotopyRec,
	Toggle2,
	pos,A,g,h,i,x,tmp, htpy,
#################################
AbsInt,                         #
SignInt;                        #
                                #
AbsInt:=AbsInt_HAP;             #
SignInt:=SignInt_HAP;           #
#################################


G:=arg[1];
n:=arg[2];
tmp:=SSortedList(Factors(Order(G)));

##A HACK April 2017
if n=0 then x:=ResolutionFiniteGroup(G,0);
i:=Position(x!.properties,["characteristic", 0 ]);
x!.properties[i][2]:=tmp[1];
return x;
fi;
##

if Length(tmp)>1 and Length(arg)<3 then 
Print("The third input variable must be a prime when this function is applied to non prime-power nilpotent groups. \n");
return fail;
fi;

if Length(arg)>2 then prime:=arg[3];
 if not IsNilpotentGroup(arg[1]) then Print("The group is not nilpotent.\n");
 return fail; fi;
else prime:=tmp[1]; fi;

pp:=Order(G);
one:=Identity(GaloisField(prime));
zero:=0*one;

gensG:=ReduceGenerators(GeneratorsOfGroup(G),G);
eltsG:=Elements(G);
if Order(eltsG[1])>1 then
eltsG:=List(eltsG,x->x);
pos:=Position(eltsG,One(G));
eltsG[pos]:=eltsG[1];
eltsG[1]:=One(G);
fi;
MT:=MultiplicationTable(eltsG);

pcgens:=Pcgs(SylowSubgroup(G,prime));
pcgens:=ReduceGenerators(pcgens,G);
pcgens:=List(pcgens,x->Position(eltsG,x));

PseudoBoundary:=[];
PseudoBoundaryAsVec:=[];
for i in [1..n] do
PseudoBoundary[i]:=[];
PseudoBoundaryAsVec[i]:=[];
od;
ImageFBasis:=[1..n];
#####################################################################
Dimension:=function(i);
if i<0 then return 0; fi;
if i=0 then return 1; fi;
return Length(PseudoBoundaryAsVec[i]);
end;
#####################################################################

#####################################################################
Boundary:=function(i,j);

if i<=0 then return []; fi;
if j>0 then
return PseudoBoundary[i][j]; 
else return
 NegateWord(PseudoBoundary[i][-j]);
 fi;

end;
#####################################################################

#####################################################################
WordToVectorList:=function(w,k)	#w is a word in R_k. 
local v,x,a;			#v is a list of vectors mod p.
#v:=List([1..Dimension(k)],x->List([1..pp],y->0) );
#v:=ListWithIdenticalEntries(Dimension(k),ListWithIdenticalEntries(pp,0) );
#How many more times will I make this mistake?

#v:=List([1..Dimension(k)],x-> ListWithIdenticalEntries(pp,0)      );
v:=NullMat(Dimension(k),pp);
for x in w do
a:=AbsoluteValue(x[1]);
v[a][x[2]]:=v[a][x[2]] + SignInt(x[1]);
od;

return v mod prime;
end;
#####################################################################

for x in gensG do
Add(PseudoBoundary[1], [[-1,1],[1,Position(eltsG,x)]]  );
Add(PseudoBoundaryAsVec[1], Flat(WordToVectorList
( [[-1,1],[1,Position(eltsG,x)]] ,0))*one);
od;

#####################################################################
VectorListToWord:=function(v)
local w, i, x;

w:=[];
for x in [1..Length(v)] do
for i in [1..Length(v[x])] do
if not v[x][i]=0 then 
Append(w, MultiplyWord(v[x][i],[ [x,i]   ]));
fi;
od;
od;

return w;
end;
#####################################################################


#####################################################################
GactMat:=function(g,tB)
local k,q,h,C,ppp;

ppp:=pp;

C:=[];

k:=0;
for q in [0..(-1+Length(tB)/ppp)] do

for h in [1..ppp] do
C[k+MT[g][h]]:=tB[k+h];
od;
k:=k+ppp;
od;

ConvertToMatrixRepNC(C);

return C;
end;
#####################################################################



#####################################################################
BoundaryMatrix:=function(k)	#Returns the matrix of d_k:R_k->R_k-1
local B,M,r, b, i, g,j,gB;		
				#M is actually the transpose of the matrix!
B:=TransposedMat(PseudoBoundaryAsVec[k]);

M:=[];

for g in [1..pp] do
gB:=TransposedMat(GactMat(g,B));

for i in [0..Dimension(k)-1] do
M[i*pp+g]:=gB[i+1];
od;
od;


return M;
end;
#####################################################################


#####################################################################
InverseFlat:=function(v)
local w,x,cnt,i;


w:=[];
cnt:=0;
while cnt< Length(v) do
x:=[];
for i in [cnt+1..cnt+pp] do
Add(x,v[i]);
cnt:=cnt+1;
od;
Add(w,x);
od;

return w;
end;
#####################################################################

#####################################################################
ComplementaryBasis:=function(arg)
local B, NS, BC, heads,ln, i, v,zeroheads;


B:=arg[1];
if Length(arg)>1 then
NS:=arg[2];
fi;

ConvertToMatrixRepNC(B,prime); 
heads:=SemiEchelonMat(B).heads;
ln:=Length(B[1]);
BC:=[];

zeroheads:=Filtered([1..ln],i->heads[i]=0);

if Length(arg)=1 then
	for i in zeroheads do
	v:=ListWithIdenticalEntries(ln,zero);
	v[i]:=one;
	ConvertToVectorRep(v,prime);
	Add(BC,v);
	od;
else
	for i in zeroheads do
	v:=NS.vectors[NS.heads[i]];
	ConvertToVectorRep(v,prime);
	Add(BC,v);
	od;
fi;



return BC;
end;
#####################################################################

#####################################################################
ZGbasisOfKernel:=function(k)		#The workhorse!
local  	tB,i, v, g, h, b,bb, ln, B, B1, B2,NS, 
	Bcomp, BndMat, ReducedB1,rrb, rrb1;

BndMat:=BoundaryMatrix(k); 


ConvertToMatrixRepNC(BndMat,prime); 
NS:=SemiEchelonMat(NullspaceMat(BndMat));

tB:=TransposedMat(NS.vectors);
Bcomp:=ComplementaryBasis(NS.vectors);


for g in pcgens do     	 
Append(Bcomp,SemiEchelonMat(TransposedMat(tB-GactMat(g,tB))).vectors);
od;							



B1:=ComplementaryBasis(Bcomp,NS);
NS:=Length(NS.vectors);

if not IsPrimePowerInt(Order(G)) then
	ReducedB1:=[];
	rrb:=0;
	rrb1:=0;
		for v in B1 do
		if rrb1=NS then break; fi;
		rrb1:=FpGModule
		(Concatenation(ReducedB1,[v]),G,prime)!.dimension;
		if rrb1>rrb then
		rrb:=rrb1;
		ReducedB1:=Concatenation(ReducedB1,[v]);
		fi;
		od;


B1:=StructuralCopy(ReducedB1);
		for v in B1 do
		ReducedB1:=Filtered(ReducedB1,x->not x=v);
		rrb:=FpGModule
                (ReducedB1,G,prime)!.dimension;
		if rrb<NS then Add(ReducedB1,v); fi;
		od;

B1:=ReducedB1;

fi;

Apply(B1,v->List(v,x->IntFFE(x)));
return B1;
end;
#####################################################################

for i in [2..n] do


for x in ZGbasisOfKernel(i-1) do
Add(PseudoBoundary[i], VectorListToWord(InverseFlat(x))   );
Add(PseudoBoundaryAsVec[i], x*one   );
od;

od;


#####################################################################

#####################################################################

#####################################################################

BoundaryMatrices:=[];

#####################################################################
Echelonize:=function()
local i, Mt, T;

#BoundaryMatrices:=[];
BoundaryMatrices2:=[];

for i in [1..n] do
BoundaryMatrices2[i]:=one*BoundaryMatrix(i);
ConvertToMatrixRep(BoundaryMatrices2[i],prime);
BoundaryMatrices[i]:=TransposedMat(BoundaryMatrices2[i]);
od;

EchelonMatrices:=[];

for i in [1..n] do
#We want to solve XM=W, so work on (Mt)(Xt)=(Wt)
# where
ConvertToMatrixRepNC(BoundaryMatrices[i],prime);
A:=TransposedMat(BoundaryMatrices[i]);
ConvertToMatrixRepNC(A);
#T:=SemiEchelonMatTransformation(TransposedMat(BoundaryMatrices[i]));
T:=SemiEchelonMatTransformation(A);
ConvertToMatrixRepNC(T.vectors);
ConvertToMatrixRepNC(T.coeffs);
ConvertToMatrixRepNC(T.relations);
EchelonMatrices[i]:=[T,Length(A),Length(A[1])];
od;

end;
#####################################################################

#####################################################################
InitHomotopyRec:=function()
local i,j,k;

HomotopyRec:=[];
for i in [1..n] do
HomotopyRec[i]:=[];
for j in [1..Dimension(i-1)] do
HomotopyRec[i][j]:=[];
for k in [1..pp] do
HomotopyRec[i][j][k]:=0;
od;
od;
od;

end;
#####################################################################


#####################################################################
SolutionMatBoundaryMatrices:=function(m,vec)
local i,ncols,sem, vno, z,x, row, sol;

ncols := Length(vec);
z := zero;
sol := ListWithIdenticalEntries(EchelonMatrices[m][2],z);
ConvertToVectorRepNC(sol);
sem := EchelonMatrices[m][1];
for i in [1..ncols] do
vno := sem.heads[i];
if vno <> 0 then
x := vec[i];
if x <> z then
AddRowVector(vec, sem.vectors[vno], -x);
 AddRowVector(sol, sem.coeffs[vno], x);
       fi;
       fi;
    od;
       if IsZero(vec) then

      return sol;
        fi;

end;
#####################################################################


Toggle:=true;
Toggle2:=true;
#####################################################################
SMBM:=function(m,w);

if Toggle then Echelonize(); Toggle:=false; fi;

return SolutionMatBoundaryMatrices(m,w);
end;
#####################################################################

#####################################################################
SMBMhom:=function(m,w) #Need SMBH as a vector space homomorphism
local coeffs, answer,i;

return SMBM(m,w); #I think SMBH is alread a group homomorphism!!!!!!

if IsInt(ImageFBasis[m]) then
ImageFBasis[m]:=(SemiEchelonMat(BoundaryMatrices2[m]).vectors);
ConvertToMatrixRep(ImageFBasis[m]);
fi;

coeffs:=SolutionMat(ImageFBasis[m],w);
#answer:=List([1..Dimension(m)*pp], a->zero);
answer:=ListWithIdenticalEntries(Dimension(m)*pp,zero);
ConvertToVectorRep(answer);
for i in [1..Length(coeffs)] do
if not IsZero(coeffs[i]) then
answer:=answer + coeffs[i]*SMBM(m,ImageFBasis[m][i]);
fi;
od;

return answer;
end;
#####################################################################

#####################################################################
Homotopy:=function(k,w)         #assume w is in kernel d_n
local u,v,s,Ab;

Ab:=AbsInt(w[1]);
if Toggle2 then InitHomotopyRec(); Toggle2:=false; fi;

if not IsInt(HomotopyRec[k+1][Ab][w[2]]) then
if SignInt(w[1]) > 0 then

return 1*HomotopyRec[k+1][Ab][w[2]]; 
else

return 1*NegateWord(HomotopyRec[k+1][Ab][w[2]]);
fi;
fi;


if Toggle then Echelonize(); Toggle:=false; fi;

v:=Flat(WordToVectorList([w],k))*one;

ConvertToVectorRep(v,prime);
if k=0 then 

s:=Sum(v);
u:=ListWithIdenticalEntries(Length(v),zero);
u[1]:=s;
ConvertToVectorRep(u);
v:=v-u;
v:=SMBM(k+1,v);

else
v:=v-SMBM(k,
 v*(BoundaryMatrices2[k]));

v:=SMBM(k+1,v);

fi;
Apply(v,i->IntFFE(i));

v:=VectorListToWord(InverseFlat(v)); 

if SignInt(w[1])>0 then
HomotopyRec[k+1][w[1]][w[2]]:=v;
else
HomotopyRec[k+1][-w[1]][w[2]]:=NegateWord(v);
fi;

return 1*v;
end;
#####################################################################



return Objectify(HapResolution,
	        rec(
		dimension:=Dimension,
		boundary:=Boundary ,
		boundaryMatrices:=BoundaryMatrices,
		homotopy:=Homotopy,
		elts:=eltsG,
		group:=G,
		properties:=
			[["length",n],
			 ["reduced",true],
			 ["type","resolution"],
			 ["characteristic",prime],
			 ["isMinimal",true]],
		solutionMatBoundaryMatrices:=
		  SMBM,
		actMat:=GactMat
		  ));
end);
#####################################################################
#####################################################################



#####################################################################
#####################################################################
InstallGlobalFunction(RankHomologyPGroup,
function(arg)
local
	G, R, N, Ranks, expand, pol, i;

if IsGroup(arg[1]) and Length(arg)=3 then
G:=arg[1];
pol:=arg[2];
N:=arg[3];
Ranks:=[];
Ranks[1]:=Length(AbelianInvariants(G));
expand:=ExpansionOfRationalFunction(pol,N);
for i in [2..N] do
Ranks[i]:= expand[i+1]- Ranks[i-1];
od;
return Ranks[N] ;
fi;


##############################################
if IsHapResolution(arg[1]) then
if EvaluateProperty(arg[1],"isMinimal")=true then 
R:=arg[1]; G:=R!.group; N:=arg[2];
fi;
fi;
##############################################

##############################################
if IsGroup(arg[1]) and Length(arg)=2 then
G:=arg[1]; N:=arg[2];
R:=ResolutionPrimePowerGroup(G,N);
fi;
##############################################

if not IsBound(R) then
 return fail; 
fi;

Ranks:=[];
Ranks[1]:=Length(AbelianInvariants(G));
for i in [2..N] do
Ranks[i]:=R!.dimension(i) - Ranks[i-1];
od;
return Ranks[N];

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(NumberGeneratorsOfGroupHomology,
function(arg)
local
	G,p,n,N,pol,
	Ranks,
	expand,
	i;

G:=arg[1];
p:=arg[2];
n:=arg[3];
if Length(arg)>3 then N:=arg[4];
else N:=n; fi;

pol:=PoincareSeriesPrimePart(G,p,n);

Ranks:=[];
Ranks[1]:=Length(AbelianInvariants(SylowSubgroup(G/DerivedSubgroup(G),p)));
expand:=ExpansionOfRationalFunction(pol,N);

for i in [2..N] do
Ranks[i]:= expand[i+1]- Ranks[i-1];
od;

return Ranks;

end);
#####################################################################
#####################################################################
