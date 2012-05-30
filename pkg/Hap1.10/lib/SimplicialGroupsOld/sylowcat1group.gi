InstallGlobalFunction(SylowSubgroupOfCatOneGroup2,
function(CC,p)
local
    C,G,P,gens,s,t,sp,tp,
	Num,i,k;
    C:=XmodToHAP(CC);
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
	gens:=GeneratorsOfGroup(P);
    sp:=GroupHomomorphismByImages(P,P,gens,List(gens,x->Image(s,x)));
    tp:=GroupHomomorphismByImages(P,P,
    gens,List(gens,x->Image(t,x)));
    return Objectify(HapCatOneGroup, rec(
                         sourceMap:=sp,
                         targetMap:=tp,
                         ));
end);