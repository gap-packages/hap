#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(HomToModPModule,
function(arg)
local R,f,p,k,ag, HomObj, HomArr,Image;

R:=arg[1];
f:=arg[2];
if Length(arg)=3 then ag:=arg[3]; else ag:=fail; fi;
k:=FieldOfMatrixGroup(Target(f));
p:=Characteristic(k);
if not k=GF(p) then
Print("The ground field is not of prime order.\n");
return fail; fi;

if IsBound(f!.fun) then
####################################
Image:=function(f,x);
return f!.fun(x);
end;
####################################
#else Image:=ImageElm; 
else Image:=ImagesRepresentative;
fi;

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
	#sum := ListWithIdenticalEntries(LA, 0);
         sum:=[1..LA]*0;

	for x in R!.boundary(n+1,i) do   ##Stupidly inefficient!!!
		if AbsoluteValue(x[1])=kq then 
			# It is left action
			sum := sum + SignInt(x[1])*Image(f,R!.elts[x[2]]){[1..LA]}[kr];
		fi;
	od;
	Append(row,sum);

od;


if Length(row)>0 then
return row*One(GF(p));
else
return [0]*One(GF(p));
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
                ["characteristic",p]
                ]));
end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
HomArr:=function(map,f)
local R, S, C, D, DhomC, mapping,IntToPair,LA, hom, ff, rnk,
      Elts, EltsSS, recf;

hom:=map!.originalHom;
R:=Source(map);
S:=Target(map);
mapping:=map!.mapping;

#D:=f(S) --> f(R)=:C

ff:=Compose(f,hom);
if ag=fail then
C:=HomToIntegralModule(R,ff);
else
C:=HomToIntegralModule(R,f);   ##This should usually be ff
fi;
IntToPair:=C!.intToPair;
LA:=Length(Identity(Range(f)));
D:=HomToIntegralModule(S,f);
if not ag=fail then ag:=ImagesRepresentative(f,ag); fi;

Elts:=[];
EltsSS:=[];
recf:=[];
##################################
##################################
DhomC:=function(v,n)
local u,uu, m,i,j,k, bnd, bnd2, x,z,zz,posD, posC, pos;

u:=[1..C!.dimension(n)]*0;

for j in [1..C!.dimension(n)/LA] do

bnd:=mapping([[j,1]],n);

   for x in bnd do
   z:=SignInt(x[1])*v{[(AbsInt(x[1])-1)*LA+1..AbsInt(x[1])*LA  ]};

if not x[2] in EltsSS then AddSet(EltsSS,x[2]); Add(Elts,x[2]);
   pos:=Length(Elts);recf[pos]:=Image(f,S!.elts[x[2]]);
   recf[pos]:=TransposedMat(recf[pos]);
else pos:=Position(Elts,x[2]); 
fi;

   #zz:=Image(f,S!.elts[x[2]])*z;
   zz:=z*recf[pos];

      for k in [1..LA] do
      posC:=(j-1)*LA + k;
      u[posC] := u[posC] + zz[k];
      od;
   od;
od;

if not ag=fail then uu:=[];
   for j in [1..C!.dimension(n)/LA] do
   Append(uu,  u{[(j-1)*LA+1..j*LA]}*ag);
   od;
u:=uu;
fi;

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
           ["characteristic", p]
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

						
