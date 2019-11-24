#########################################
#########################################
InstallGlobalFunction(TruncatedRegularCWComplex,
function(Y,nn)
local B, n;

n:=Minimum(nn,Dimension(Y));
B:=Y!.boundaries;
B:=B{[1..n+1]};
B[n+2]:=[];

return RegularCWComplex(B);
end);
#########################################
#########################################

#########################################
#########################################
InstallGlobalFunction(PiZeroOfRegularCWComplex,
function(Y)
local T, cells, C, fn;

if Dimension(Y)=0 then 

###########################
fn:=function(i);
return i;
end;
###########################
return [[1..Y!.nrCells(0)], fn];
fi;

T:=TruncatedRegularCWComplex(Y,2);
cells:=CocriticalCellsOfRegularCWComplex(T,2);
cells:=Filtered(cells,x->x[1]=0);
cells:=List(cells,x->x[2]);

C:=ChainComplex(T);

###########################
fn:=function(i);
return C!.deform(0,i)[1];
end;
###########################

return [cells,fn];

end);
#########################################
#########################################

#########################################
#########################################
InstallMethod(PiZero,
"PiZero for RegularCWComplexes",
[IsHapRegularCWComplex],

function(Y)
return PiZeroOfRegularCWComplex(Y);
end);
#########################################
#########################################

#########################################
#########################################
InstallMethod(PiZero,
"PiZero for SimplicialComplexes",
[IsHapSimplicialComplex],

function(K)
local Y;
Y:=SimplicialComplexToRegularCWComplex(K);
return PiZeroOfRegularCWComplex(Y);
end);
#########################################
#########################################

#########################################
#########################################
InstallMethod(PiZero,
"PiZero for graphs",
[IsHapGraph],

function(G)
local Y;
Y:=SimplicialNerveOfGraph(G,2);
Y:=SimplicialComplexToRegularCWComplex(Y);
return PiZeroOfRegularCWComplex(Y);
end);
#########################################
#########################################

##################################################
##################################################
InstallGlobalFunction(DisplayDendrogramMat,
function(A,t,n);


DisplayDendrogram(DendrogramMat(A,t,n));

end);
##################################################
##################################################

##################################################
##################################################
InstallGlobalFunction(DendrogramMat,
function(A,tt,n)
local L, Inv, D, i, t;

t:=tt+1;

L:=[];
for i in [1..t+1] do
L[i]:=PiZero(SymmetricMatrixToGraph(A,(i-1)*n));;
od;

Inv:=List(L,l->function(i) return l[2](i); end);
L:=List(Inv,x->Classify([1..5],x));

D:=List([1..Length(L)-1],l->
        List(L[l],x->PositionProperty(L[l+1],y->x[1] in y)));

return D;

end);
##################################################
##################################################

