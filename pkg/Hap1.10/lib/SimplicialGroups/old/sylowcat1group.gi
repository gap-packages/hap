InstallGlobalFunction(SylowSubgroupOfCatOneGroup_alt,
function(C,p)
local
    G,P,Gens,s,t,sp,tp,
	Num,i,k;
	
    s:=C!.sourceMap;
    t:=C!.targetMap;
    G:=Source(s);
    P:=SylowSubgroup(Image(s,G),p);
	k:=1;
    Num:=Factors(Order(G)); 
    for i in Num do
        if i=p then
            k:=k*p;
        fi;
    od;
    while Order(P)<k do
        P:=SylowSubgroup(Normalizer(G,P),p);
    od;
	Gens:=GeneratorsOfGroup(P);
    sp:=GroupHomomorphismByImages(P,P,Gens,List(Gens,x->Image(s,x)));
    tp:=GroupHomomorphismByImages(P,P,
    Gens,List(Gens,x->Image(t,x)));
    return Objectify(HapCatOneGroup, rec(
                         sourceMap:=sp,
                         targetMap:=tp,
                         ));
end);




