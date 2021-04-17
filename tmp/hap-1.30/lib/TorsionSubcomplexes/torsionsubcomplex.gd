
DeclareGlobalFunction("TorsionSubcomplex");
DeclareGlobalFunction("DisplayAvailableCellComplexes");
DeclareGlobalFunction("VisualizeTorsionSkeleton");
DeclareGlobalFunction("ReduceTorsionSubcomplex");
DeclareGlobalFunction("GetTorsionSubcomplex");
DeclareGlobalFunction("IsPNormal");
DeclareGlobalFunction("RigidFacetsSubdivision");
DeclareGlobalFunction("EquivariantEulerCharacteristic");

DeclareGlobalFunction("CountingCellsOfACellComplex");
DeclareGlobalFunction("CountingControlledSubdividedCells");
DeclareGlobalFunction("CountingBaryCentricSubdividedCells");
DeclareGlobalFunction("GroupHomomorphismToMatrix");

DeclareGlobalFunction("EquivariantSpectralSequencePage");
DeclareGlobalFunction("ExportHapCellcomplexToDisk");
DeclareGlobalFunction("CountingCellsOfBaryCentricSubdivision");
DeclareGlobalFunction("QuotientByTorsionSubcomplex");
DeclareGlobalFunction("ExtractTorsionSubcomplex");
DeclareGlobalFunction("IntegralCellularHomology");
DeclareGlobalFunction("GammaSubgroupInSL3Z");


#####################################################################
DeclareCategory("IsHapTorsionSubcomplex",IsObject);

DeclareRepresentation(	"IsHapTorsionSubcomplexRep",
			IsComponentObjectRep,
			["torsion",
			 "groupname",
			  "reducedTorsionCells",
			  "celldata"]);

HapTorsionSubcomplexFamily:=NewFamily(	"HapTorsionSubcomplexFamily",
				IsHapTorsionSubcomplex,
				IsHapTorsionSubcomplex);

HapTorsionSubcomplex:=NewType(HapTorsionSubcomplexFamily, IsHapTorsionSubcomplexRep);

InstallMethod( ViewObj,
"for HapTorsionSubcomplex",
[IsHapTorsionSubcomplex],
 function(R)
if not R!.groupname = fail then
Print("Reduced ",R!.torsion,"-torsion subcomplex of the given cell complex for the group ",R!.groupname);
else
Print("Reduced ",R!.torsion,"-torsion subcomplex of the given cell complex for the given group");
fi;
 end);

InstallMethod( PrintObj,
"for HapTorsionSubcomplex",
[IsHapTorsionSubcomplex],
function(R)
local N, reducedTorsionCells,J,n,j,G,B,b, ReduceModP, p, celldata;

Print("Reduced ",R!.torsion,"-torsion subcomplex of the given cell complex for the group ",R!.groupname);

reducedTorsionCells:=R!.torsionCells;
#originalData:=R!.originalData;
celldata:=R!.celldata;
p:=R!.torsion;
#########
ReduceModP := function( G, p)
local K, h, H, Q, P;

	Q := G;
	K := ConjugacyClassesSubgroups(G);
	for h in K do
          if Size(h)=1 then
	       H := Representative(h);
	           if Size(H) > 1 then
	     	       if Gcd(Size(H),p) = 1 then
            	       # We pass to the quotient Q in the extension
		       # 1 -> H -> G -> Q -> 1,
		       # because H has trivial mod p cohomology.
		       # Our loop ends up using the normal subgroup H of maximal size,
		       # because the list K is ordered from small to large.
		           Q := Image( NaturalHomomorphismByNormalSubgroup( G,H ));
	              fi;
	           fi;
	    fi;
	od;
       P:=Q;
       if IsPNormal(Q,p) then P:=Normalizer(Q,Center(SylowSubgroup(Q,p)));fi;
  return P;
end;
#########
for N in [1..Size(reducedTorsionCells)] do
  for J in reducedTorsionCells[N] do
	n := J[1];
	j := J[2];
	G := celldata[n+1][j]!.TheMatrixStab;;

	Print("\n",n,"-cell number ",j,": ",IdSmallGroup(G)," is controlled by ",IdSmallGroup(ReduceModP(G,p))," and has boundary ");
 	B := celldata[n+1][j]!.BoundaryImage!.ListIFace;
	for b in B do
		G := celldata[n][b]!.TheMatrixStab;;
		Print(n-1,"-cell number ",b,": ",IdSmallGroup(G)," is controlled by ",
			IdSmallGroup(ReduceModP(G,p)),"\n");
	od;
  od;
Print("\n");
od;
end);
#####################################################################

#####################################################################
DeclareCategory("IsHapEquivariantSpectralSequencePage",IsObject);

DeclareRepresentation(	"IsHapEquivariantSpectralSequencePageRep",
			IsComponentObjectRep,
			["page",
			 "differential",
			  "groupname",
			  "torsion"]);

HapEquivariantSpectralSequencePageFamily:=NewFamily(	"HapEquivariantSpectralSequencePageFamily",
				IsHapEquivariantSpectralSequencePage,
				IsHapEquivariantSpectralSequencePage);

HapEquivariantSpectralSequencePage:=NewType(HapEquivariantSpectralSequencePageFamily, IsHapEquivariantSpectralSequencePageRep);

InstallMethod( ViewObj,
"for HapEquivariantSpectralSequencePage",
[IsHapEquivariantSpectralSequencePage],
 function(R)
Print("Equivariant Spectral Sequence for the group ",R!.groupname);
 end);

InstallMethod( PrintObj,
"for HapEquivariantSpectralSequencePage",
[IsHapEquivariantSpectralSequencePage],
function(R)
Print("Equivariant Spectral Sequence for the group ",R!.groupname);
end);
#####################################################################

DeclareGlobalFunction("KernelOfMap");
DeclareGlobalFunction("ImageOfMap");
#####################################################################
InstallGlobalFunction(KernelOfMap,
function(M)

if IsTrivial(Source(M)) then return Source(M);
else
return Kernel(M);
fi;
end);
#####################################################################
InstallGlobalFunction(ImageOfMap,
function(M)

if IsTrivial(Source(M)) then return One(Target(M));
else
return Image(M);
fi;
end);
#####################################################################
