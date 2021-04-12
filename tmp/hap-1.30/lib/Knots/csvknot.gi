
############################################################
InstallGlobalFunction(ReadCSVfileAsPureCubicalKnot, 
function(arg)
local  FILE, CNT, i,j,atoms, cnt, file, scl1, scl2, scl3, f, v,x,
       AA, A, x1, x2,y1,y2,z1,z2, L, M , MM, B;

if IsString(arg[1]) then
FILE:=[arg[1]];
else FILE:=arg[1];
fi;

A:=[];
CNT:=0;
x1:=0;x2:=0;y1:=0;y2:=0;z1:=0;z2:=0;

for file in FILE do
CNT:=CNT+1;

if  Length(arg)=2 then
scl1:=arg[2][1];
scl2:=arg[2][2];
scl3:=arg[2][3];
else
scl1:=13;
scl2:=13;
scl3:=13;
fi;


AA:=ReadCSV(file);
A[CNT]:=[];
for x in AA do
v:=[];
v[1]:=EvalString(x.A);
v[2]:=EvalString(x.B);
v[3]:=EvalString(x.C);
v[1]:=(scl1*v[1]);
v[2]:=(scl2*v[2]);
v[3]:=(scl3*v[3]);

Add(A[CNT],v);
od;
Add(A[CNT],StructuralCopy(A[CNT][1]));

atoms:=Length(A[CNT]);

if atoms=0 then return fail; fi;


Print("Reading chain containing ",atoms," atoms.\n");

for j in [1..13] do
B:=[];
for i in [1..Length(A[CNT])-1] do
B[2*i-1]:=A[CNT][i];
v:=(A[CNT][i]+A[CNT][i+1])/2;
B[2*i]:=v;
od;
B[2*i+1]:=A[CNT][i+1];
A[CNT]:=StructuralCopy(B);
B:=[];
od;

A[CNT]:=List(A[CNT],a->[Int(a[1]),Int(a[2]),Int(a[3])]);

L:=List(A[CNT],a->a[1]);
x1:=Minimum(x1,Minimum(L));
x2:=Maximum(x2,Maximum(L));
L:=List(A[CNT],a->a[2]);
y1:=Minimum(y1,Minimum(L));
y2:=Maximum(y2,Maximum(L));
L:=List(A[CNT],a->a[3]);
z1:=Minimum(z1,Minimum(L));
z2:=Maximum(z2,Maximum(L));

od;

x1:=x1-1;
y1:=y1-1;
z1:=z1-1;

M:=NullMat(y2-y1+1,z2-z1+1);
M:=List([1..x2-x1],i->StructuralCopy(M));

for i in [1..Length(A)] do
for f in A[i] do
M[f[1]-x1][f[2]-y1][f[3]-z1]:=1;
od;
od;


M:=FrameArray(M);
M:=FrameArray(M);

M:=PureCubicalComplex(M);
ContractPureCubicalComplex(M);

for i in [1..10] do
MM:=ContractedComplex(M);
MM:=RegularCWComplex(MM);;
#if Homology(MM,0)=[0] and Homology(MM,1)=[0] then return M; fi;
#if Homology(MM,0)=[0] and Length(Homology(MM,1))=Length(FILE) then return M; fi;
if Length(Homology(MM,0))=Length(FILE) and Length(Homology(MM,1))=Length(FILE) then return M; fi;
Print("thickening ...\n");
M:=ThickenedPureCubicalComplex(M);
M:=PureCubicalComplex(FrameArray(M!.binaryArray));
od;

return fail;


end);
#######################################################

############################################################
############################################################
InstallGlobalFunction(DisplayCSVknotFile,  
function(arg)
local  FILE, File, i,j,  tmpdir, file, cnt,
 scl,  f, x,  v, AA, A,  M , B, BLUE, AppendTo, PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

if IsString(arg[1]) then
FILE:=[arg[1]];
else
FILE:=arg[1];
fi;


scl:=100; 

A:=[];
cnt:=0;

for File in FILE do
cnt:=cnt+1;
AA:=ReadCSV(File);
A[cnt]:=[];

for x in AA do
v:=[];
v[1]:=EvalString(x.A);
v[2]:=EvalString(x.B);
v[3]:=EvalString(x.C);
v[1]:=(scl*v[1]);
v[2]:=(scl*v[2]);
v[3]:=(scl*v[3]);

Add(A[cnt],v);
od;

Add(A[cnt],A[cnt][1]);

A[cnt]:=List(A[cnt],a->[Int(a[1]),Int(a[2]),Int(a[3])]);

A[cnt]:=List(A[cnt],a->Concatenation("(",String(a[1]),",",String(a[2]),",",String(a[3]),")"));

od;



if Length(arg)=1 then
BLUE:=List([1..Length(FILE)],i->[]);
else BLUE:=arg[2];
fi;

tmpdir := DirectoryTemporary();;
file:=Filename( tmpdir , "tmp.asy" );

PrintTo(file, "import three;\n\n");
AppendTo(file, "size(700);\n\n");
AppendTo(file, "currentprojection=orthographic(5,15,-13);\n\n");
AppendTo(file, "defaultpen(1.0);\n\n");

for cnt in [1..Length(FILE)] do

AppendTo(file,"path3 g=", A[cnt][1]);
for i in [2..Length(A[cnt])] do
AppendTo(file, "..",A[cnt][i]);
od;
AppendTo(file,";\n\n");

AppendTo(file,"draw(g);\n");
AppendTo(file,"dot(g,red+opacity(0.8));\n");

for i in BLUE[cnt] do
AppendTo(file,"dot(" , A[cnt][i] , ", blue);\n");
od;

od;

Exec( Concatenation( ASY_PATH, "-V ", file) );

RemoveFile(file);
file:=Filename(tmpdir,"");
RemoveFile(file);

end);
#############################################################
#############################################################

