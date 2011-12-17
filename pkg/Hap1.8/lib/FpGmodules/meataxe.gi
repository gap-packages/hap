#####################################################################
#####################################################################
InstallGlobalFunction(GeneratorsOfMtxModule,
function(M)
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


G:=Group(M.generators);
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
Rad:=SemiEchelonMatDestructive(Rad);
gensM:=Mat{Filtered([1..dim],i->not i in Rad.heads)};
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
local 	DM, gensM, Mat, G, pp, elts, x, v, MT,Action, GactMat;

gensM:=GeneratorsOfMtxModule(M);
Mat:=[];
G:=Group(M.generators);
elts:=Elements(G);
pp:=Order(G);

for v in gensM do
for x in elts do
Add(Mat,v*TransposedMat(x));
#Why on earth does GAP multiply on the right?!!
od;
od;

ConvertToMatrixRep(Mat);
Mat:=NullspaceMat(Mat);


############################
############################Should make this a global function####
MT:=MultiplicationTable(elts);

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
Position(elts,g),
TransposedMat(B)));
end;
#####################################################################
############################

return Objectify(HapFPGModule,
	rec 
	(matrix:=Mat,
	group:=G,
	action:=Action,
	dimension:=Length(Mat),
	ambientDimension:=Length(gensM)*pp,
	characteristic:=Characteristic(M.field)
	));

end);
#####################################################################
#####################################################################
