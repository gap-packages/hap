# preserved most of notations from homToZmodule.gi

#####################################################################
#####################################################################
InstallGlobalFunction(TensorWithModPModule,
function(arg)
local R, f, g, p, k, one, TensorObj, TensorArr;

R:=arg[1];
f:=arg[2];
k:=FieldOfMatrixGroup(Target(f));
p:=Characteristic(k);
if not k=GF(p) then
Print("The ground field is not of prime order.\n");
return fail; fi;
one:=One(k);


TensorObj:=function(R,f)
local
	DimensionC,
	BoundaryC,
	LengthC,
	BoundaryOfElt,
	Cache,
	LA,
	IntToPair;

LA:=Length(Identity(Image(f)));
LengthC:=EvaluateProperty(R,"length");

#####################################################################
DimensionC:=function(n);
return LA*R!.dimension(n);
end;
#####################################################################

####################################################################
IntToPair:=function(i)
local q, r;
r:=i mod LA;
q:=(i-r)/LA;

if r>0 then return [q+1,r];
else return [q, LA]; fi;
end;
####################################################################

#####################################################################
BoundaryOfElt:=function(n,k)	#Only use this for k>0
local row, kq, kr, i, x, r, v;

if n<=0 then
	return [];
fi;

x:=IntToPair(k);
kq:=x[1]; kr:=x[2];

row:= [1..LA * R!.dimension(n-1)]*0;

for x in R!.boundary(n,kq) do
	i := AbsoluteValue(x[1]);

	# range of generators of tensored complex which coresponds to generator i of resolution
	r := [ (i-1)*LA+1 .. i*LA ];
        v:=r*0; v[kr]:=1;

	# It is left action.
	# Would be faster:
	# for ?
	# list of lists and then Flat() ?
	#row{r} := row{r} + SignInt(x[1]) * Image(f,(R!.elts[x[2]])^-1){[1..LA]}[kr];
row{r} := row{r} + SignInt(x[1]) * v*TransposedMat((Image(f,R!.elts[x[2]]^-1)));

od;

return row*one;
end;
#####################################################################

#####################################################################
BoundaryC:=function(n,k)

if n <= 0 then
	return [];
fi;

if not IsBound(Cache) then
	Cache := [];
fi;
if not IsBound(Cache[n]) then
	Cache[n] := [];
fi;
if not IsBound(Cache[n][k]) then
	Cache[n][k] := BoundaryOfElt(n,k);
	MakeImmutable(Cache[n][k]); # is it working as I want ?
fi;

return Cache[n][k];
end;
#####################################################################


return Objectify(HapChainComplex,
		rec(
                dimension:=DimensionC,
                boundary:=BoundaryC,
                IntToPair:=IntToPair,
                properties:=
                [["length",LengthC],
                ["connected",true],
                ["type", "chainComplex"], 
                ["characteristic",p]
                 ]));
end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
TensorArr:=function(F,f)
local
                R,S,g,                    #R->S is an equivariant chain
                C,D,ChomD,              #map.
                DimensionS,
                DimensionC,
                x, gg, mapgen, LA;
R:=F!.source;
S:=F!.target;
g:=F!.conjugator;
DimensionS:=S!.dimension;
C:=TensorObj(R,f);
D:=TensorObj(S,f);
DimensionC:=C!.dimension;
LA:=Length(Identity(Image(f)));


#####################################################################
mapgen:=function(n,k)
local w, v, x, y, z, u, i, zz;
w:=[1..D!.dimension(n)]*0;
x:=C!.IntToPair(k);
v:=[1..LA]*0; v[x[2]]:=1;
y:=F!.mapping([[x[1],1]],n);

gg:=TransposedMat(Image(f,g))^-1;
for z in y do
u:=v*gg*TransposedMat(Image(f,S!.elts[z[2]]^-1));
u:=SignInt(z[1])*u;
zz:=AbsInt(z[1]);
for i in [(zz-1)*LA+1..zz*LA] do
w[i]:=w[i]+u[i-(zz-1)*LA];
od;
od;
return w;
end;
#####################################################################


#####################################################################
ChomD:=function(v,n)
local w,k;
w:=[1..D!.dimension(n)]*0;
for k in [1..Length(v)] do
if not v[k]=0 then
w:=w+v[k]*mapgen(n,k);
fi;
od;
return w;
end;
#####################################################################


return Objectify(HapChainMap,
        rec(
           source:=C,
           target:=D,
           mapping:=ChomD,
           properties:=[ ["type","chainMap"],
           ["characteristic", EvaluateProperty(C,"characteristic")] 
           ]));
end;
#####################################################################
#####################################################################


if IsHapResolution(R) then
return TensorObj(R,f);
fi;
if IsHapEquivariantChainMap(R) then
return TensorArr(R,f);
fi;

return fail;
						
end);
#####################################################################
#####################################################################

