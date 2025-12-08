##################################################
##################################################
4x4to2x2:=function(R,d)
local k, i, Q, OQ, I, G, n;
n:=Length(R);

  for k in [1..n] do
    for i in [1..R!.dimension(k)] do
    R!.boundary(k,i);
    od;
  od;
  Q:=QuadraticNumberField(d);;OQ:=RingOfIntegers(Q);;I:=QuadraticIdeal(OQ,1);;
  Apply(R!.elts,x->HAP_4x4MatTo2x2Mat(x,d));
  G:=HAP_CongruenceSubgroupGamma0(I);;
  G!.tree:=true;
  R!.group:=G;

end;
###################################################
###################################################
