#C 2008 Graham Ellis

#####################################################################
InstallMethod(HAPPrintTo,
      "Method for writing HAP resolutions to a file",
	[IsString,IsHapResolution],

function(file,R)
local n,k,g;

#Create a gap function
PrintTo(file,"HAPTEMPORARYFUNCTION:=function() local TYPE,RANKS, HIGHEST_DEGREE,BOUNDARIES,CONTRACTING_HOMOTOPY,ELEMENTS,SPECIAL_PROPERTIES;\n\n");
#Export type of HAP object 
AppendTo(file,"TYPE:=\nHapResolution;\n\n");

#Export length of resolution
AppendTo(file,"HIGHEST_DEGREE:= \n");
AppendTo(file,Length(R),";\n\n");

#Export ranks of the resolution in each degree.
AppendTo(file,"RANKS:=\n");
AppendTo(file,List([0..Length(R)],i->Dimension(R)(i)),";\n\n");

#Export boundary words
AppendTo(file,"#BOUNDARIES[n][k] is the boundary of the kth free generator in degree n\n");

AppendTo(file,"BOUNDARIES:=","\n[\n");
for n in [1..Length(R)] do
AppendTo(file,"#Degree ", n,"\n","[\n\n"  );
for k in [1..R!.dimension(n)] do
AppendTo(file,R!.boundary(n,k));
if k<R!.dimension(n) then AppendTo(file,",\n\n");
else AppendTo(file,"\n\n");
fi;
od;
if n<Length(R) then
AppendTo(file,"],\n\n");
else
AppendTo(file,"]\n];\n\n");
fi;
od;

#Export the contracting homotopy
AppendTo(file,"#CONTRATCTING_HOMOTOPY[n][k][g] is the image of the free\n");
AppendTo(file,"#abelian group generator g(e_k)^n in degree n.\n");
if R!.homotopy=fail then
AppendTo(file, "CONTRACTING_HOMOTOPY:=\nfail;\n\n");
else
AppendTo(file, "CONTRACTING_HOMOTOPY:=\n[\n\n");
for n in [0..Length(R)-2] do
AppendTo(file,"#Degree ", n,"\n","[\n\n"  );
for k in [1..R!.dimension(n)] do
AppendTo(file,"#Free Generator ",k," in degree ",n," \n[\n\n");
  for g in [1..Length(R!.elts)] do
AppendTo(file,R!.homotopy(n,[k,g]));
AppendTo(file,",\n\n");
od;
AppendTo(file,"\n],\n\n");
od;
if n<Length(R)-2 then
AppendTo(file,"],\n\n");
else
AppendTo(file,"]\n];\n\n");
fi;
od;
fi;

#Export group elements
AppendTo(file, "ELEMENTS:=\n");
AppendTo(file, R!.elts,";\n\n");

#Export special properties of the resolution
AppendTo(file,"SPECIAL_PROPERTIES:=\n", R!.properties,";\n\n");

AppendTo(file, "return rec( type:=TYPE, highest_degree:=HIGHEST_DEGREE, ranks:=RANKS, boundaries:=BOUNDARIES,contracting_homotopy:=CONTRACTING_HOMOTOPY, elements:=ELEMENTS,special_properties:=SPECIAL_PROPERTIES);\n\n");

AppendTo(file,"end;\n");
end);
