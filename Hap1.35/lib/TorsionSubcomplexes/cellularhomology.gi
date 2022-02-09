InstallGlobalFunction(IntegralCellularHomology,
function(C)
local lnth, Cells, nrCells,x,w,j,s,id, Mat,t,M,b,p,h,r,
      Mult, CLeftCosetElt, pos, hom, Elts, MatRank, y, BdMat,
      KillGAction, Reduction;




##################################################################
MatRank:=function(g)
if not (IsBound(g[1]) and IsBound(g[1][1])) then return 0;
else return RankMat(g);fi;
end;

if IsBound(C!.length) then lnth:=C!.length;
else
    lnth:=0;
    while not C!.dimension(lnth)=0 do
        lnth:=lnth+1;
    od;
    lnth:=lnth-1;
fi;
###################################################################
KillGAction:=function(L)
  return List(L,x->x[1]);
end;
###################################################################
Reduction:=function(w)
local v,x,k;
v := Filtered( w, x->x>0 );
 for x  in w  do
     if x < 0  then
         k := Position( v, -x );
         if k = fail  then
             Add( v, x );
         else
             Unbind( v[k] );
         fi;
     fi;
 od;
 return Filtered( v, function ( a )
         return IsBound( a);end);
end;
###################################################################
Mat:=[];
BdMat:=[];
for j in [2..lnth+1] do
  M:=[];
  for s in [1..C!.dimension(j-1)] do
    M[s]:=[];
    for t in [1..C!.dimension(j-2)] do
      M[s][t]:=0;
    od;
  od;
#Print(1/0);
  for s in [1..C!.dimension(j-1)] do

    b:=Collected(Reduction(KillGAction(C!.boundary(j-1,s))));

    for t in b do

        M[s][AbsInt(t[1])]:=SignInt(t[1])*t[2];
    od;
  od;
#  Print("Representation matrix at degree ",j-1,"  ",M,"\n");
  Mat[j]:=SmithNormalFormIntegerMat(TransposedMat(M));
  BdMat[j]:=M;


od;

# Create Mat[1] for d_0
M:=[];
for s in [1..C!.dimension(0)] do
    M[s]:=0;
od;
Mat[1]:=TransposedMat([M]);
BdMat[1]:=TransposedMat([M]);

# Create Mat[lnth+1] for d_(lnth+1)
M:=[];
for s in [1..C!.dimension(lnth)] do
    M[s]:=0;
od;
Mat[lnth+2]:=[M];
BdMat[lnth+2]:=[M];
hom:=[];
for j in [1..lnth+1] do
    h:=[];
    r:=C!.dimension(j-1)-MatRank(Mat[j])-MatRank(Mat[j+1]);
    for s in [1..r] do
  h[s]:=0;
    od;
    for s in [1..MatRank(Mat[j+1])] do
        if not Mat[j+1][s][s]=1 then
            Add(h,Mat[j+1][s][s]);
        fi;
    od;
    Add(hom,h);
od;
#Print(1/0);

 return hom;



end);

########################################
