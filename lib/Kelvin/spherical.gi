
##############################################################
##############################################################
## inputs an arc presentation and returns the complement of the link
## in the 3-sphere S^3 (rather than the complement in R^3).
InstallGlobalFunction(SphericalKnotComplement,
function(arg)
local arc, K, bm, lst, paths, Smap, S, B, bnds;

arc:=arg[1];
if Length(arg)=1 then
   K:=KnotComplement(arc);
else  #This seems to produce a slightly smaller CW structure
   Print("Using alternative method. \n");
   K:=ArcPresentationToKnottedOneComplex(arc);
   K:=RegularCWComplexComplement(K);
fi;
K:=SimplifiedComplex(K);
bm:=BoundaryMap(K);
lst:=RegularCWMapToCWSubcomplex(bm);
paths:=PathComponentsCWSubcomplex(lst);
paths:=List(paths,CWSubcomplexToRegularCWMap);

for Smap in paths do
if Homology(Source(Smap),1)=[] then break; fi;
od;

#Smap is a map from the sphere to the boundary of the knot complement K

S:=Source(Smap);
B:=List([1..S!.nrCells(2)], i->Smap!.mapping(2,i));
B:=Concatenation([Length(B)],B);

bnds:=1*K!.boundaries;
Add(bnds[4],B);

return SimplifiedComplex(RegularCWComplex(bnds));
end);
##############################################################
##############################################################

