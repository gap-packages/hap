#(C) Graham Ellis, 2005-2006


#####################################################################
#####################################################################
InstallGlobalFunction(ResolutionFpGModule,
function(arg) 
local
	MDL,G,n,
	eltsG,
	gensG,
	Dimension,
	CorrectedDimension,
	Boundary,
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
	Mgens,
	one,
	InverseFlat,
	ComplementaryBasis,
	zero,
	pcgens,
	BoundaryMatrices,
	Echelonize,
	EchelonMatrices,
	SolutionMatBoundaryMatrices,
	Toggle,
	SMBM,
	t,g,h,i,x,tmp;


MDL:=arg[1];
G:=MDL!.group;

if not IsPGroup(G) then
Print("Group must be of prime-power order.\n");
return fail;
fi;

n:=arg[2]+1;
tmp:=SSortedList(Factors(Order(G)));


prime:=MDL!.characteristic;

pp:=Order(G);
one:=Identity(GaloisField(prime));
zero:=0*one;

gensG:=ReduceGenerators(GeneratorsOfGroup(G),G);
eltsG:=Elements(G);
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

#####################################################################
Dimension:=function(i);
if i<0 then return 0; fi;
return Length(PseudoBoundaryAsVec[i]);
end;
#####################################################################

#####################################################################
CorrectedDimension:=function(t)
local i;
i:=t+1;
if i<0 then return 0; fi;
return Length(PseudoBoundaryAsVec[i]);
end;
#####################################################################

#####################################################################
Boundary:=function(t,j)
local i;
i:=t+1;
if i<0 then return []; fi;
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

v:=ListWithIdenticalEntries(CorrectedDimension(k),ListWithIdenticalEntries(pp,0) );


for x in w do
a:=AbsoluteValue(x[1]);
v[a][x[2]]:=v[a][x[2]] + SignInt(x[1]);
od;

return v mod prime;
end;
#####################################################################


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

Mgens:=GeneratorsOfFpGModule(MDL);
for x in Mgens do
tmp:=[];
for i in [1..Length(x)/pp] do
tmp[i]:=x{[1+(i-1)*pp..i*pp]};
tmp[i]:=List(tmp[i],j->IntFFE(j));
od;
Add(PseudoBoundary[1], VectorListToWord(tmp));
Add(PseudoBoundaryAsVec[1],x);
od;



#####################################################################
GactMat:=function(g,tB)
local k,q,h,C;

C:=[];

k:=0;
for q in [0..(-1+Length(tB)/pp)] do

for h in [1..pp] do
C[k+MT[g][h]]:=tB[k+h];
od;
k:=k+pp;
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
	Bcomp, BndMat, ReducedB1, rrb, rrb1;

BndMat:=BoundaryMatrix(k); 


ConvertToMatrixRepNC(BndMat,prime); 
NS:=NullspaceMat(BndMat);
#if Length(NS)>0 then     #Need to handle case =0 at some stage!!
NS:=SemiEchelonMat(NS);
#fi;

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

B1:=List(B1,v->List(v,x->IntFFE(x)));
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


#####################################################################
Echelonize:=function()
local i,  Mt, T;

BoundaryMatrices:=[];
for i in [1..n] do
BoundaryMatrices[i]:=TransposedMat(BoundaryMatrix(i)*one);
od;

EchelonMatrices:=[];

for i in [1..n] do
#We want to solve XM=W, so work on (Mt)(Xt)=(Wt)
# where
ConvertToMatrixRepNC(BoundaryMatrices[i],prime);
T:=SemiEchelonMatTransformation(BoundaryMatrices[i]);
EchelonMatrices[i]:=[T.coeffs,Reversed(T.heads),T.vectors];
od;

end;
#####################################################################


#####################################################################
SolutionMatBoundaryMatrices:=function(m,w)
local h,u,v,i,cnt,pos,col,row,diff;

ConvertToVectorRep(w);
v:=EchelonMatrices[m][1]*w;
h:=StructuralCopy(EchelonMatrices[m][2]);
u:=ListWithIdenticalEntries(Length(h),0*one);

while not Sum(h)=0 do

col:=PositionProperty(h,x->not x=0);
row:=h[col];
h[col]:=0;
diff:=EchelonMatrices[m][3][row]*u;
pos:=Length(u)+1-col;
u[pos]:=v[row]-diff;
ConvertToVectorRep(u);
od;


return u;
end;
#####################################################################


Toggle:=true;
#####################################################################
SMBM:=function(m,w);

if Toggle then Echelonize(); Toggle:=false; fi;

return SolutionMatBoundaryMatrices(m,w);
end;
#####################################################################

#####################################################################
Homotopy:=function(k,w)         
local u,v,s,uu,i;

if Toggle then Echelonize(); Toggle:=false; fi;

v:=Flat(WordToVectorList([w],k))*one;
ConvertToVectorRep(v,prime);


if k=0 then 


u:=SignInt(w[1])*Mgens[AbsInt(w[1])];
u:=MDL!.action(eltsG[w[2]],[u]);
u:=u[1];
u:=FpGModuleSection(MDL,u);
u:=WordToVectorList(u,k);
u:=Flat(u);

v:=v-u;
v:=SolutionMatBoundaryMatrices(k+1,v);

else

v:=v-SolutionMatBoundaryMatrices(k,BoundaryMatrices[k]*v);
v:=SolutionMatBoundaryMatrices(k+1,v);

fi;
Apply(v,i->IntFFE(i));

if not v=fail then v:=VectorListToWord(InverseFlat(v)); fi;
return v;
end;
#####################################################################



return Objectify(HapResolution,
	        rec(
		dimension:=CorrectedDimension,
		boundary:=Boundary,
		homotopy:=Homotopy,
		elts:=eltsG,
		group:=G,
		properties:=
			[["length",n-1], 
			 ["reduced",true],
			 ["type","resolution"],
			 ["characteristic",prime],
			 ["isMinimal",true]],
		solutionMatBoundaryMatrices:=
		  SMBM));
end);
#####################################################################
#####################################################################



