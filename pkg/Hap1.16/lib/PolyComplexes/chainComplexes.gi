#(C) Graham Ellis

#################################################################
#################################################################
InstallGlobalFunction(CoreducedChainComplex,
function(arg)
local
        C,D,
        S,T,
        Lngth,
        NewBoundary,
        Bnd,
        FinalBoundary,
        NewDimension,
        n,i,pp;

C:=arg[1];
Bnd:=[];

##########################
if Length(arg)=2 then
if arg[2]=2 then return 2CoreducedChainComplex(C);fi;
fi;
##########################

Lngth:=EvaluateProperty(C,"length");
S:=[];T:=[];
for n in [0..Lngth] do
S[n+1]:=[1..C!.dimension(n)];
T[n+1]:=[];
od;

############################
NewDimension:=function(n);
return Length(S[n+1]);
end;
############################

############################
NewBoundary:=function(n,i)
local v,j;
if n=0 then return C!.boundary(n,i);fi;
v:= StructuralCopy(C!.boundary(n,i));
for j in T[n] do
v[j]:=0;
od;
return v;
end;
############################

############################
FinalBoundary:=function(n,i) local j;

if not IsBound(Bnd[n]) then Bnd[n]:=[]; fi;
if not IsBound(Bnd[n][i]) then
Bnd[n][i]:= NewBoundary(n,S[n+1][i]); 
for j in T[n] do
Bnd[n][i][j]:=true;
od;
Bnd[n][i]:=Filtered(Bnd[n][i],a->not a=true);

fi;

return Bnd[n][i];

end;
############################

############################
for n in [1..Lngth] do
for i in [1..C!.dimension(n)] do
if AbsInt(Sum(List(NewBoundary(n,i),x->AbsInt(x))))=1 then 

RemoveSet(S[n+1],i);
pp:=PositionProperty(NewBoundary(n,i),a->not a=0);
RemoveSet(S[n],pp);
Add(T[n],pp);
Add(T[n+1],i);

fi;
od;
od;
############################

D:=rec( dimension:=NewDimension,
        boundary:=FinalBoundary,
        properties:=StructuralCopy(C!.properties)
        );

D:=Objectify(HapChainComplex, D);

return D;
end);
#################################################################
#################################################################



#################################################################
#################################################################
InstallGlobalFunction(2CoreducedChainComplex,
function(C)
local
        D,
        S,T,TCopy,
        Lngth,
        NewBoundary,
        FinalBoundary,
        NewDimension,
	Dels,DELS,Smt,
        n,i,j,x,pp,ppp;

Lngth:=EvaluateProperty(C,"length");
S:=[];T:=[];
for n in [0..Lngth] do
S[n+1]:=[1..C!.dimension(n)];
T[n+1]:=[];
od;

############################
NewDimension:=function(n);
return Length(S[n+1]);
end;
############################

############################
NewBoundary:=function(n,i)
local v,j;
if n=0 then return C!.boundary(n,i);fi;
v:= StructuralCopy(C!.boundary(n,i));
for j in T[n] do
v[j]:=0;
od;
return v;
end;
############################

############################
FinalBoundary:=function(n,i);

return NewBoundary(n,S[n+1][i]);

end;
############################

############################
for n in [1..Lngth] do
Dels:=[];
DELS:=[];

for i in [1..C!.dimension(n)] do
if Length(Filtered(NewBoundary(n,i),x->not x=0))=2 then Add(Dels,i);fi;
od;


for i in Dels do
for j in Dels{[Position(Dels,i)+1..Length(Dels)]} do
if 
Length(Filtered(
 List(NewBoundary(n,i),x->AbsInt(x))+List( NewBoundary(n,j),x->AbsInt(x)),
x->not x=0))=2 
then Add(DELS,[i,j]);
fi;
od;
od;


for x in DELS do
i:=x[1];j:=x[2];
Smt:=SmithNormalFormIntegerMat([NewBoundary(n,i),NewBoundary(n,j)]);
Smt[1]:=List(Smt[1],x->AbsInt(x));
Smt[2]:=List(Smt[2],x->AbsInt(x));
if Sum(Smt[1])=1 and Sum(Smt[2])=1 then

RemoveSet(S[n+1],i);
RemoveSet(S[n+1],j);
pp:=Position(Smt[1],1);
RemoveSet(S[n],pp);
ppp:=Position(Smt[2],1);
RemoveSet(S[n],ppp);

Add(T[n],pp);
Add(T[n],ppp);
Add(T[n+1],i);
Add(T[n+1],j);
fi;
od;

od;
############################

D:=rec( dimension:=NewDimension,
        boundary:=FinalBoundary,
        properties:=StructuralCopy(C!.properties)
        );

D:=Objectify(HapChainComplex, D);

return D;
end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(SuspendedChainComplex,
function(C)
local 
	D,
	L;

D:=Objectify(HapChainComplex,
rec( properties:=StructuralCopy(C!.properties) ));

######################
D!.dimension:=function(n);
if n=0 then return 0; fi;
return C!.dimension(n-1);
end;
######################

######################
D!.boundary:=function(n,i);
if n=0 then return C!.boundary(0,i); fi;
return C!.boundary(n-1,i);
end;
######################

L:=PositionProperty(C!.properties,x->"length" in x);
D!.properties[L][2]:=C!.properties[L][2]+1;
return D;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ReducedSuspendedChainComplex,
function(C)
local
        D,
        L;

D:=Objectify(HapChainComplex,
rec( properties:=StructuralCopy(C!.properties) ));

######################
D!.dimension:=function(n);
if n=0 then return 0; fi;
return C!.dimension(n-1);
end;
######################

######################
D!.boundary:=function(n,i);
if n=0 then return C!.boundary(0,i); fi;
if n=1 then return [1]; fi;
return C!.boundary(n-1,i);
end;
######################

L:=PositionProperty(C!.properties,x->"length" in x);
D!.properties[L][2]:=C!.properties[L][2]+1;
return D;
end);
#################################################################
#################################################################

