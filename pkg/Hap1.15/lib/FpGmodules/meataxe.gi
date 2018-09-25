#####################################################################
#####################################################################
InstallGlobalFunction(GeneratorsOfMtxModule,
function(M) # It is assumed that G acts on column vectors on the 
	    # left (g,v)-->g.v in all FG-modules. 
local
	G,
        dim,
        P,
        prime,
        Mat,
        one,
        gensP,
        Rad,
        x,v,
        gensM,
        RedgensM,
        DimOfMod,
        d1,d2;

#####################################################################
DimOfMod:=function(A)
local B,v,x;

if Length(A)=0 then return 0;fi;

B:=[];
for x in G do
Append(B,A*x);
od;

ConvertToMatrixRep(B);
return Length(SemiEchelonMatDestructive(B).vectors);
end;
#####################################################################


G:=Group(List(M.generators,x->TransposedMat(x)));
# This should handle to left-right module conversion.

dim:=M.dimension;
prime:=Characteristic(M.field);
one:=One(GF(prime));
Mat:=IdentityMat(dim)*one;
ConvertToMatrixRep(Mat); #The ROWS of MAT are a basis for the module.

P:=SylowSubgroup(G,prime);
gensP:=GeneratorsOfGroup(P);
gensP:=ReduceGenerators(gensP,P);
for x in gensP do
ConvertToMatrixRep(x);
od;

Rad:=[];
for x in gensP do
Append(Rad,Mat-Mat*x);
od;

if Length(Rad)=0 then Add(Rad,0*Mat[1]); fi;
Rad:=SemiEchelonMat(Rad);

gensM:=Mat{Filtered([1..dim],i->Rad.heads[i]=0)};
Rad:=0;

##################
if not IsPrimePowerInt(Order(G)) then

RedgensM:=[];
d1:=0;
for v in gensM do
        Add(RedgensM,v);
        d2:=DimOfMod(RedgensM);
        if d2=d1 then
        RedgensM:=Filtered(RedgensM,x->not x=v);
        else
        if d2=dim then break; fi;
        fi;
od;

gensM:=StructuralCopy(RedgensM);
for v in gensM do
RedgensM:=Filtered(RedgensM,x->not x=v);
if DimOfMod(RedgensM)<dim then Add(RedgensM,v); fi;
od;

gensM:=RedgensM;
fi;
##############

return gensM;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(DesuspensionMtxModule,
function(M)
local 	gensM, Mat, G, GG,hom,  pp, elts, elts1,elts2,x, v, MT,Action, GactMat;

gensM:=GeneratorsOfMtxModule(M);

G:=Group(M.generators);
if "gens" in NamesOfComponents(M) then
GG:=Group(M.gens); else GG:=G; M.gens:=M.generators; fi;
hom:=GroupHomomorphismByImages(GG,G,M.gens,M.generators);

Mat:=[];
elts1:=Elements(GG);;
elts:=List(elts1,x->(Image(hom,x)));
elts2:=List(elts,x->TransposedMat(x));

pp:=Order(G);

for v in gensM do
for x in elts2 do
Add(Mat,v*x);
# This is really x*Transpose(v) !!
od;
od;


ConvertToMatrixRep(Mat);
Mat:=NullspaceMat(Mat);


############################
############################Should make this a global function####
MT:=MultiplicationTable(elts1);

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
Action:=function(g,B);
return TransposedMat(GactMat(
Position(elts1,g),
TransposedMat(B)));
end;
#####################################################################
############################

return Objectify(HapFPGModule,
	rec 
	(matrix:=Mat,
	group:=GG,
	action:=Action,
	dimension:=Length(Mat),
	ambientDimension:=Length(gensM)*pp,
	characteristic:=Characteristic(M.field)
	));

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(FpG_to_MtxModule,
function(M) 
local
	MDL,
	G,
	gensG,
	MatgensG,
	prime,
	F,
	B,
	x;

prime:=M!.characteristic;
F:=GF(prime);

G:=M!.group;
gensG:=GeneratorsOfGroup(G);
MatgensG:=[];

for x in gensG do
B:=MutableCopyMat(M!.action(x,M!.matrix));
Add(MatgensG, 
SolutionsMatDestructive(MutableCopyMat(M!.matrix),B));
od;

MatgensG:=List(MatgensG,x->TransposedMat(x));

MDL:=GModuleByMats(MatgensG,F);
MDL.gens:=gensG;

return MDL;
end);
#####################################################################
#####################################################################


