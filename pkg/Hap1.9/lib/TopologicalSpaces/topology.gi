#C 2007 Graham Ellis 

#####################################################################
#####################################################################
InstallGlobalFunction(MatrixToTopologicalManifold,
function(AA,threshold)	
local
	rows,cols,
	i,j,A;
A:=StructuralCopy(AA);
rows:=Length(A);
cols:=Length(A[1]);

if IsInt(threshold) then
for i in [1..rows] do
for j in [1..cols] do
if A[i][j] <threshold then A[i][j]:=1; else
A[i][j]:=0;fi;
od;od;
fi;

if IsList(threshold) then
for i in [1..rows] do
for j in [1..cols] do
if A[i][j] >threshold[1] and A[i][j] <threshold[2] then A[i][j]:=1; else
A[i][j]:=0;fi;
od;od;
fi;




return Objectify(HapTopologicalManifold,
	   rec(
	   BinaryList:=A,
	   properties:=[
	   ["dimension",2],
	   ["arraySize",[rows,cols]],
	   ["genericCWtype","cubical"],
	   ["CW",false]]
	   ));
end);
#####################################################################
#####################################################################

HAPAAA:=0;
#################################################################
#################################################################
InstallGlobalFunction(ReadImageAsTopologicalManifold,
function(file,threshold)
local i,j,prog,A;

prog:=Concatenation(GAP_ROOT_PATHS[1],"pkg/Hap1.9/lib/TopologicalSpaces/prog");

i:=Concatenation("convert ",file," /tmp/im.txt");
Exec(i);
i:=Concatenation("perl ",prog," /tmp/im.txt >/tmp/im.g");
Exec(i);

Read("/tmp/im.g");
Exec("rm /tmp/im.g");
Exec("rm /tmp/im.txt");
A:=StructuralCopy(HAPAAA);
HAPAAA:=0;

return MatrixToTopologicalManifold(A,threshold);
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(WriteTopologicalManifoldAsImage,
function(T,file,ext)
local
        A,i,j,rows,cols,colour,filetxt;

if not EvaluateProperty(T,"dimension")=2 then
Print("There is no method for viewing a topological manifold of dimension ",
EvaluateProperty(T,"dimension"),".\n"); return fail; fi;

A:=T!.BinaryList;
rows:=Length(A);;
cols:=Length(A[1]);

filetxt:=Concatenation(file,".txt");

PrintTo(filetxt,"# ImageMagick pixel enumeration: ",
Length(A),",",Length(A[1]),",255,RGB\n");

for i in [1..rows] do
for j in [1..cols] do
if A[i][j]=0 then colour:="(255,255,255)"; else colour:="(0,0,100)";fi;
AppendTo(filetxt,i,",",j,": ",colour,"\n");
od;
od;

i:=Concatenation("convert ",filetxt," ",file,".",ext);
Exec(i);
i:=Concatenation("rm ",filetxt);
Exec(i);

end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ViewTopologicalManifold,
function(arg)
local i,A,viewer,T;

T:=arg[1];
if not EvaluateProperty(T,"dimension")=2 then
Print("There is no method for viewing a topological maqnifold of dimension ",
EvaluateProperty(T,"dimension"),".\n"); return fail; fi;

if Length(arg)>1 then viewer:=arg[2];
else viewer:="mozilla";fi;
WriteTopologicalManifoldAsImage(T,"HAPtmpImage","png");
Exec(Concatenation(viewer," ","HAPtmpImage.png"));
Sleep(2);
Exec(Concatenation("rm  ","HAPtmpImage.png"));

end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ConcatenatedTopologicalManifold,
function(L)
local
	S,
	dim,i,
	genericCWtype,
	row,col,arraySize;

dim:=EvaluateProperty(L[1],"dimension");
genericCWtype:=EvaluateProperty(L[1],"genericCWtype");
arraySize:=EvaluateProperty(L[1],"arraySize");
for i in [2..Length(L)] do
if not 
arraySize=EvaluateProperty(L[i],"arraySize")
and genericCWtype=EvaluateProperty(L[i],"genericCWtype")
then
Print("The concatenated manifolds must all have similar dimensions and CW-structures. \n");
return fail; fi;
od;

return
Objectify(HapTopologicalManifold,
           rec(
	        BinaryList:=List(L,T->StructuralCopy(T!.BinaryList)),
	        properties:=[
	        ["dimension",dim+1],
	        ["genericCWtype",genericCWtype],
		["arraySize",Concatenation([Length(L)],arraySize)],
	        ["CW",false]]
	       ));

end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ZerothBettiNumber,
function(T)
local 
	ChoosePoint,
	ArraySize,ArraySizePO,
	Tarray,
	RawData,
	SmallCircle,SmallCircumference,
	Neighbours,
	Dim,
	point,colour,
	neibs,newneibs,
	lenpoint,
	i,tmp;


tmp:= EvaluateProperty(T,"pathComponents");
if not tmp=fail then return tmp; fi;

Dim:=EvaluateProperty(T,"dimension");
ArraySize:=EvaluateProperty(T,"arraySize");
ArraySizePO:=ArraySize+List([1..Length(ArraySize)],i->1);
SmallCircle:=Cartesian(List([1..Dim],i->[-1,0,1]));
SmallCircle:=Filtered(SmallCircle,x->not IsZero(x));
SmallCircumference:=Length(SmallCircle);
RawData:=T!.BinaryList;

Tarray:=Cartesian(List(ArraySize,x->[1..x]));


#################################################
ChoosePoint:=function(RawData,Tarray)
local p,Tpoint,i;

for p in Tarray do

if EvalT(RawData,ArraySizePO,p)=1 then return p; fi;
od;
return fail;
end;
#################################################

#################################################
Neighbours:=function(RawData,ArraySizePO,SmallCircle,SmallCircumference,point)
local nbs;

nbs:=List([1..SmallCircumference],i->point)+SmallCircle;
return Filtered(nbs,x->EvalT(RawData,ArraySizePO,x)=1);
end;
#################################################

point:=ChoosePoint(RawData,Tarray);


if not point=fail then
lenpoint:=Length(point)-1;
fi;
colour:=1;

while not point=fail do
	colour:=colour+1;
	SetT(RawData,ArraySizePO,point,lenpoint,colour);

neibs:=Neighbours(RawData,ArraySizePO,SmallCircle,SmallCircumference,point);

	while Length(neibs)>0 do
	for point in neibs do
	SetT(RawData,ArraySizePO,point,lenpoint,colour);
	od;

	newneibs:=Concatenation(List(neibs,x->
Neighbours(RawData,ArraySizePO,SmallCircle,SmallCircumference,x)
	));
        neibs:=SSortedList(newneibs);
	od;
point:=ChoosePoint(RawData,Tarray);

od;
Add(T!.properties,["pathComponents",colour-1]);

return colour-1;
end);
################################################################
################################################################

################################################################
################################################################
InstallOtherMethod(Bettinumbers,
"betti numbers of a topological Manifold",
[IsHapTopologicalManifold,IsInt],

function(T,n);
if n=0 then return ZerothBettiNumber(T); fi;
if n=1 then return BettiNumbersOfMatrix(
ComplementTopologicalManifold(ComplementTopologicalManifold(T))!.BinaryList)[2]; fi;
if n>EvaluateProperty(T,"dimension") then return 0;fi;

Print("Theis function is not yet defined for n= ",n,"\n");
return fail;
end);
################################################################
################################################################

################################################################
################################################################
InstallOtherMethod(Bettinumbers,
"betti numbers of a topological Manifold",
[IsHapTopologicalManifold],

function(T);
if EvaluateProperty(T,"dimension")=2 then
return BettiNumbersOfMatrix(
ComplementTopologicalManifold(ComplementTopologicalManifold(T))!.BinaryList); fi;

Print("This function is not yet defined for manifolds if dimension >2.","\n");
return fail;
end);
################################################################
################################################################


################################################################
################################################################
InstallGlobalFunction(PathComponent,
function(T,n)
local A,PC,lp,RawData, N, Dim, Tarray,ArraySize, ArraySizePO, point;

N:=ZerothBettiNumber(T);
N:=n+1;
Dim:=EvaluateProperty(T,"dimension");
ArraySize:=EvaluateProperty(T,"arraySize");
ArraySizePO:=ArraySize+List([1..Length(ArraySize)],i->1);
RawData:=StructuralCopy(T!.BinaryList);
Tarray:=Cartesian(List(ArraySize,x->[1..x]));

A:=0*T!.BinaryList;
lp:=Dim-1;
for point in Tarray do
if EvalT(RawData,ArraySizePO,point)=N
then SetT(A,ArraySizePO,point,lp,1);
else SetT(A,ArraySizePO,point,lp,0);
fi;
od;

return 
Objectify(HapTopologicalManifold,
           rec(
                BinaryList:=A,
                properties:=[
                ["dimension",Dim],
                ["genericCWtype",EvaluateProperty(T,"genericCWtype")],
                ["arraySize",ArraySize],
                ["CW",EvaluateProperty(T,"CW")]]
                 ));

end);
################################################################
################################################################

######################################################################
######################################################################
InstallGlobalFunction(ThickenedTopologicalManifold,
function(T)
local 
        ArraySize,ArraySizePO,
        Tarray,
        RawData,
        SmallCircle,SmallCircumference,
        Neighbours,
        Dim,
	A,lp,point,p;

Dim:=EvaluateProperty(T,"dimension");
lp:=Dim-1;
ArraySize:=EvaluateProperty(T,"arraySize");
ArraySizePO:=ArraySize+List([1..Length(ArraySize)],i->1);
SmallCircle:=Cartesian(List([1..Dim],i->[-1,0,1]));
SmallCircle:=Filtered(SmallCircle,x->not IsZero(x));
SmallCircumference:=Length(SmallCircle);
RawData:=1*(T!.BinaryList);
Tarray:=Cartesian(List(ArraySize,x->[1..x]));


#################################################
Neighbours:=function(SmallCircle,SmallCircumference,point)
local nbs;

return List([1..SmallCircumference],i->point)+SmallCircle;
end;
#################################################

for point in Tarray do
if not EvalT(T!.BinaryList,ArraySizePO,point)=0 then
SetT(RawData,ArraySizePO,point,lp,1);
for p in Neighbours(SmallCircle,SmallCircumference,point) do
SetT(RawData,ArraySizePO,p,lp,1); 
od;
fi;
od;

return
Objectify(HapTopologicalManifold,
           rec(
               BinaryList:=RawData,
               properties:=[
               ["dimension",Dim],
               ["genericCWtype",EvaluateProperty(T,"genericCWtype")],
               ["arraySize",ArraySize],
               ["CW",EvaluateProperty(T,"CW")]]
               ));
end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(ComplementTopologicalManifold,
function(T)
local
        ArraySize,ArraySizePO,
        Tarray,
        RawData,
        Dim,
        lp,point;

Dim:=EvaluateProperty(T,"dimension");
lp:=Dim-1;
ArraySize:=EvaluateProperty(T,"arraySize");
ArraySizePO:=ArraySize+List([1..Length(ArraySize)],i->1);
RawData:=1*(T!.BinaryList);
Tarray:=Cartesian(List(ArraySize,x->[1..x]));

for point in Tarray do
if EvalT(RawData,ArraySizePO,point)=0 then
SetT(RawData,ArraySizePO,point,lp,1);
else SetT(RawData,ArraySizePO,point,lp,0); fi;
od;

return
Objectify(HapTopologicalManifold,
           rec(
	       BinaryList:=RawData,
	       properties:=[
	       ["dimension",Dim],
	       ["genericCWtype",EvaluateProperty(T,"genericCWtype")],
	       ["arraySize",ArraySize],
	       ["CW",EvaluateProperty(T,"CW")]]
	       ));

end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(SetT,
function(RawData,ArraySizePO,p,lp,value)
local i,Tvalue;

if 0 in p then return 0;fi;
if 0 in (p-ArraySizePO) then return 0;fi;

Tvalue:=RawData;
for i in [1..lp] do
Tvalue:=Tvalue[p[i]];
od;
Tvalue[p[Length(p)]]:=value;

end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(EvalT,
function(RawData,ArraySizePO,p)
local i,Tvalue;

if 0 in p then return 0;fi;
if 0 in (p-ArraySizePO) then return 0;fi;
Tvalue:=RawData;
for i in p do
Tvalue:=Tvalue[i];
od;
return Tvalue;

end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(BoundaryTopologicalManifold,
function(T)
local
	C,
	Dim,
	lp,
	ArraySize,
	ArraySizePO,
	Tarray,
	point,
	RawData;

Dim:=EvaluateProperty(T,"dimension");
lp:=Dim-1;
ArraySize:=EvaluateProperty(T,"arraySize");
ArraySizePO:=ArraySize+List([1..Length(ArraySize)],i->1);
Tarray:=Cartesian(List(ArraySize,x->[1..x]));

C:=ComplementTopologicalManifold(T);
C:=ThickenedTopologicalManifold(C);
RawData:=C!.BinaryList;

for point in Tarray do
if EvalT(T!.BinaryList,ArraySizePO,point)=0 then
SetT(RawData,ArraySizePO,point,lp,0);
fi;
od;

C!.BinaryList:=RawData;
return C;
end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(ContractTopologicalManifold,
function(T);

if not EvaluateProperty(T,"dimension")=2 then
Print("This function is not yet implemented for manifolds of dimension ",
EvaluateProperty(T,"dimension"),"\n");
return fail; fi;

ContractMatrix(T!.BinaryList);

end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(SingularChainComplex,
function(T);

if not EvaluateProperty(T,"dimension")=2 then 
Print("This function is not yet implemented for manifolds of dimension ",
EvaluateProperty(T,"dimension"),"\n");
return fail; fi;

return MatrixToChainComplex(T!.BinaryList);
end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(BoundarySingularities,
function(T)

local A;

if not EvaluateProperty(T,"dimension")=2 then
Print("This function is not yet implemented for manifolds of dimension ",
EvaluateProperty(T,"dimension"),"\n");
return fail; fi;

A:=SingularityMatrix(T!.BinaryList);
return 

Objectify(HapTopologicalManifold,
          rec(
               BinaryList:=A,
               properties:=[
             ["dimension",2],
             ["arraySize",[Length(A),Length(A[1])]],
             ["genericCWtype","cubical"],
             ["CW",false]]
		                 ));

end);
######################################################################
######################################################################

