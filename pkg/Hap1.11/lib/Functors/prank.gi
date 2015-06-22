Prank:=function(G)
local S, p;

p:=SSortedList(Factors(Order(G)));

if not Length(p)=1 then 
Print("G must be a finite p-group. \n");
return fail;
fi;

p:=p[1];
S:=LatticeSubgroups(G);
S:=ConjugacyClassesSubgroups(S);
S:=List(S,x->ClassElementLattice(x,1));
S:=Filtered(S,x->IsElementaryAbelian(x));
S:=List(S,x->Order(x));

return Log(Maximum(S),p);

end;
