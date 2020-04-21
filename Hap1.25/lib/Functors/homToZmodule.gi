#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(HomToIntegralModule,
function(R,f)
local HomObj, HomArr,Image;

####################################
Image:=function(f,x);
return f!.fun(x);
end;
####################################

#####################################################################
#####################################################################
HomObj:=function(R,f)
local
	 DimensionC,
	 BoundaryC,
         LengthC,
	 BoundaryOfElt,
         M, M0,
	 LA,
	 IntToPair;

LA:=Length(Identity(Range(f)));
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
                intToPair:=IntToPair,
                properties:=
                [["length",LengthC],
                ["connected",true],
                ["type", "cochainComplex"],
                ["characteristic",
                EvaluateProperty(R,"characteristic")] ]));
end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
HomArr:=function(map,f)
local R, S, C, D, DhomC, mapping,IntToPair,LA;

R:=Source(map);
S:=Target(map);
mapping:=map!.mapping;

#D:=f(S) --> f(R)=:C

C:=HomToIntegralModule(R,f);
IntToPair:=C!.intToPair;
LA:=Length(Identity(Range(f)));
D:=HomToIntegralModule(S,f);

##################################
##################################
DhomC:=function(v,n)
local u, m,i,j,k, bnd, x,z,zz,posD, posC;
u:=[1..C!.dimension(n)]*0;


for j in [1..C!.dimension(n)/LA] do
bnd:=mapping([[j,1]],n);
   for x in bnd do
   z:=SignInt(x[1])*v{[(AbsInt(x[1])-1)*LA+1..AbsInt(x[1])*LA  ]};
   z:=TransposedMat([z]);
   zz:=Image(f,S!.elts[x[2]])*z;
   zz:=TransposedMat(zz); zz:=zz[1];
      for k in [1..LA] do
      posC:=(j-1)*LA + k;
      u[posC] := u[posC] + zz[k];
      od;
   od;
od;



return u;
end;
##################################
##################################

return Objectify(HapCochainMap,
        rec(
           source:=D,
           target:=C,
           mapping:=DhomC,
           properties:=[ ["type","cochainMap"],
           ["characteristic", Maximum(
           EvaluateProperty(D,"characteristic"),
           EvaluateProperty(C,"characteristic"))]
           ]));

end;
#####################################################################
#####################################################################


if IsHapResolution(R) then
return HomObj(R,f);
fi;

if IsHapEquivariantChainMap(R) then
return HomArr(R,f);
fi;


end);
#####################################################################
#####################################################################

						
