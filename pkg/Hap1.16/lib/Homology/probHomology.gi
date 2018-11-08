#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(HomologyPb,
function(C,n)
local
        M1, M2,
        Dimension, Boundary,
        dr, drtmp,
	i,j,ans,dim1,dim2,v,x,AppendTo,PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

if n <=0 then return Homology(C,n); fi;
#if n=0 then return [0]; fi;  #THIS IS MATHEMATICALLY WRONG!!!

if not "snf" in NamesOfComponents(C) then
C!.snf:=[1..EvaluateProperty(C,"length")];
fi;

Dimension:=C!.dimension;
if Dimension(n)=0 then return []; fi;
Boundary:=C!.boundary;
M1:=[];
M2:=[];ans:=[];

if IsInt(C!.snf[n]) then
drtmp:=DirectoryTemporary();
dr:=Filename(drtmp,"HAPtmp");
AppendTo(dr,Dimension(n)," ",Length(Boundary(n,1))," ","M \n");
for i in [1..Dimension(n)] do
#M1[i]:=Boundary(n,i);
v:=Boundary(n,i); 
for x in [1..Length(v)] do
if not v[x]=0 then
AppendTo(dr,i," ",x," ",v[x], " \n");
fi;
od;
od;

AppendTo(dr,0," ",0," ",0, " \n");
M1:=SMInvariantFactors(dr);
dim1:=M1[1]-Sum(List([3..Length(M1)-1],a->M1[a][2]));
C!.snf[n]:=M1;
dr:=Filename(drtmp,"");
Exec(Concatenation("rm -r ",dr));

else
M1:=C!.snf[n];
dim1:=M1[1]-Sum(List([3..Length(M1)-1],a->M1[a][2]));

fi;

if IsInt(C!.snf[n+1]) then

if Dimension(n+1)=0 then M2:=[0,0]; C!.snf[n+1]:=M2;

else
drtmp:=DirectoryTemporary();
dr:=Filename(drtmp,"HAPtmp");
AppendTo(dr,Dimension(n+1)," ",Length(Boundary(n+1,1))," ","M \n");
for i in [1..Dimension(n+1)] do
#M2[i]:=Boundary(n+1,i);
v:=Boundary(n+1,i);
for x in [1..Length(v)] do
if not v[x]=0 then
AppendTo(dr,i," ",x," ",v[x], " \n");
fi;
od;
od;
AppendTo(dr,0," ",0," ",0, " \n");

M2:=SMInvariantFactors(dr);
dr:=Filename(drtmp,"");
Exec(Concatenation("rm -r ",dr));

fi;
C!.snf[n+1]:=M2;
dim2:=0;

for i in [3..Length(M2)] do
if M2[i][1]>0 then dim2:=dim2+M2[i][2]; fi;
if M2[i][1]>1 then Append(ans,List([1..M2[i][2]],x->M2[i][1])); fi;
od;

else
M2:=C!.snf[n+1];
dim2:=0;

for i in [3..Length(M2)] do
if M2[i][1]>0 then dim2:=dim2+M2[i][2]; fi;
if M2[i][1]>1 then Append(ans,List([1..M2[i][2]],x->M2[i][1])); fi;
od;
fi;

for i in [1..dim1-dim2] do
Append(ans,[0]);
od;
return ans;
end );
#####################################################################
