
##########################################################
##########################################################
InstallMethod(HAP_MultiplicativeGenerators,
" internal method for expresing a basis in terms of ring generators",
[IsAlgebra],
function(A)
local gensA, BasA, BasAexp, grgensA, x, i, W, n, P, w, v,
      AsGens, mx; 

#A:=ModPCohomologyRing(G,n) for some p-group G

gensA:=ModPRingGenerators(A);
BasA:=Filtered(gensA,x->A!.degree(x)=0);
BasAexp:=[[BasA[1]]];
mx:=Maximum(List(gensA,A!.degree));

grgensA:=List([1..Maximum(List(gensA,A!.degree))],i->[]);

for x in gensA do
    i:=A!.degree(x);
    if i>0 then Add(grgensA[i],x); fi;
od;

W:=Submodule(A,BasA);


n:=0;
while Dimension(W)<Dimension(A) do
    n:=n+1;
    for P in Partitions(n) do
        if Maximum(Flat(P)) <= mx then
        for w in Cartesian(grgensA{P}) do
            v:=Product(w);
            if not v in W then Add(BasA,v); W:=Submodule(A,BasA); 
            Add(BasAexp,w); fi;
        od;
        fi;
    od;
od;

BasA:=Basis(A,BasA);

###################################
AsGens:=function(v)
local w;
w:=Coefficients(BasA,v);
w:=Filtered([1..Length(BasA)],i->not IsZero(w[i]));
return BasAexp{w};
end;
###################################

return [BasAexp,AsGens];
end);
##################################################################
##################################################################

