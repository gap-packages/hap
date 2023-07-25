

######################################################
######################################################
InstallGlobalFunction(HAP_ChainComplexToEquivariantChainComplex,
function(C)
local R, G, Elts, Boundary, i, x, m;

G:=Group( () );
Elts:=Elements(G);

Boundary:=function(n,k)
local v, bnd;
v:=C!.boundary(n,k);
bnd:=[];

for i in [1..Length(v)] do
if v[i]<>0 then x:=[SignInt(v[i])*i,1];
for m in [1..AbsInt(v[i])] do
Add(bnd,x);
od;
fi;
od;
return bnd;
end;

R:=Objectify(HapResolution,
                rec(
                dimension:=C!.dimension,
                boundary:=Boundary,
                homotopy:=fail,
                elts:=Elts,
                group:=G,
                vectorField:=fail,
                properties:=
                   [["length",EvaluateProperty(C,"length")],
                    ["reduced",true],
                    ["type","resolution"],
                    ["characteristic",0]  ]));


return R;
end);
######################################################
######################################################


###################################################################
###################################################################
InstallGlobalFunction(HAP_bockstein,
function(A)
local G,i,prime,x,a,b,c,h,s,bhomc,B,C,psi,R,dim,cons,homs,bas,bok,gens;

G:=Group(());
prime:=Characteristic(A);
x:=[2..prime^2];
Add(x,1);
x:=PermList(x);

a:=Group(x^prime);;
b:=Group(x);;
c:=Group(x^prime);;

#bhomc:=NaturalHomomorphismByNormalSubgroup(b,a);
B:=TrivialGModuleAsGOuterGroup(G,b);
#C:=TrivialGModuleAsGOuterGroup(G,Image(bhomc));
C:=TrivialGModuleAsGOuterGroup(G,c);
bhomc:=GroupHomomorphismByImages(b,c,[x],[x^prime]);
psi:=GOuterGroupHomomorphism();
psi!.Source:=B;
psi!.Target:=C;
psi!.Mapping:=bhomc;

C:=A!.chainComplex;
dim:=Length(C);
R:=HAP_ChainComplexToEquivariantChainComplex(C);
R!.properties[1]:=["length",dim+1];

cons:=[];
for i in [1..dim-1] do
Add(cons, ConnectingCohomologyHomomorphism(psi,i,R));
od;
homs:=[];
for h in cons do
s:=Source(h);
Add(homs,s!.nat);
od;
s:=Target(h);
Add(homs,s!.nat);

bas:=Basis(A);
gens:=List(homs,h->Target(h!.Mapping));
gens:=List(gens,g->GeneratorsOfGroup(g));

########################################
bok:=function(x)
local w,ww,h,z,cnt;

w:=Coefficients(bas,x);
ww:=[];
for i in [2..Length(w)] do
Add(ww,[IntFFE(w[i]),A!.int2pair[i]]- [0,[1,0]]);
od;


ww:=Filtered(ww,x->x[2][1]<dim);

h:=[];
for i in [1..Length(ww)] do
Add(h,[ww[i][2][1]+1 ,
 
 Image(cons[ww[i][2][1]]!.Mapping, gens[ww[i][2][1]][ww[i][2][2]])^ww[i][1] 

] );
od;

h:=List(h,x->[x[1], ExtRepOfObj(
Factorization(Source(homs[x[1]]!.Mapping),
PreImagesRepresentative(homs[x[1]]!.Mapping, x[2]))
)
]);

h:=Filtered(h,z->Length(z[2])>0);
w:=0*A.1;
for z in h do
cnt:=A!.pair2int[z[1]+1][1]-1;

w:=w+
Sum(List([1..Length(z[2])/2],k-> bas[ cnt+ z[2][2*k-1]  ]*z[2][2*k]));
;
od;

return w;

end;
########################################

return bok;
end);
###################################################################
###################################################################





