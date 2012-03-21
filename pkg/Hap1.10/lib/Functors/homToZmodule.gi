#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(HomToIntegralModule,
function(R,f)
local
	 DimensionC,
	 BoundaryC,
         LengthC,
	 BoundaryOfElt,
         M, M0,
	 LA,
	 IntToPair;

LA:=Length(Identity(Image(f)));
LengthC:=EvaluateProperty(R,"length");
M:=[1..LengthC]; M0:=0;

#####################################################################
DimensionC:=function(n);
return LA*R!.dimension(n);
end;
#####################################################################

#####################################################################
IntToPair:=function(i)
local q, r;
r:=i mod LA;
q:=(i-r)/LA;

if r>0 then return [q+1,r];
else return [q, LA]; fi;
end;
#####################################################################

#####################################################################
BoundaryOfElt:=function(n,k)	#Only use this for k>0
local
        row, kq, kr, i, x, sum;

if n<0 then return List([1..DimensionC(0)],a->0); fi;

#if n=0 then return List([1..DimensionC(1)],x->0); fi;

x:=IntToPair(k);
kq:=x[1]; kr:=x[2];

row:=[];
for i in [1..R!.dimension(n+1)] do
	sum := ListWithIdenticalEntries(LA, 0);

	for x in R!.boundary(n+1,i) do
		if AbsoluteValue(x[1])=kq then 
			# It is left action
			sum := sum + SignInt(x[1])*Image(f,R!.elts[x[2]]){[1..LA]}[kr];
		fi;
	od;
	Append(row,sum);

od;


if Length(row)>0 then
return row;
else
return [0];
fi;
end;
#####################################################################

#####################################################################
BoundaryC:=function(n,k)
local Mt,i,row,j;

############ CASE n=0 ####################
if n=0 then 
if M0=0 then 
Mt:=[];
for i in [1..DimensionC(n)] do
Append(Mt, [BoundaryOfElt(n,i)]);
od;
M0:=Mt;
fi;
return M0[k];
fi;
########### CASE n=0 DONE ###############



if  M[n]=n then
Mt:=[];


	for i in [1..DimensionC(n)] do
	Append(Mt, [BoundaryOfElt(n,i)]);
	od;
	M[n]:=Mt;

fi;



return M[n][k];
end;
#####################################################################


return Objectify(HapCochainComplex,
		rec(
                dimension:=DimensionC,
                boundary:=BoundaryC,
                properties:=
                [["length",LengthC],
                ["connected",true],
                ["type", "cochainComplex"],
                ["characteristic",
                EvaluateProperty(R,"characteristic")] ]));
						
end);
#####################################################################
#####################################################################

