# preserved most of notations from homToZmodule.gi

#####################################################################
#####################################################################
InstallGlobalFunction(TensorWithIntegralModule,
function(R,f)
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
local row, kq, kr, i, x, r;

if n<=0 then
	return [];
fi;

x:=IntToPair(k);
kq:=x[1]; kr:=x[2];

row:= ListWithIdenticalEntries(LA * R!.dimension(n-1), 0);

for x in R!.boundary(n,kq) do
	i := AbsoluteValue(x[1]);

	# range of generators of tensored complex which coresponds to generator i of resolution
	r := [ (i-1)*LA+1 .. i*LA ];

	# It is left action.
	# Would be faster:
	# for ?
	# list of lists and then Flat() ?
	row{r} := row{r} + SignInt(x[1]) * Image(f,(R!.elts[x[2]])^-1){[1..LA]}[kr];
od;

return row;
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
                properties:=
                [["length",LengthC],
                ["connected",true],
                ["type", "chainComplex"], 
                ["characteristic",
                EvaluateProperty(R,"characteristic")] ]));
						
end);
#####################################################################
#####################################################################

