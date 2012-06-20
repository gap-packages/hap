#(C) Graham Ellis, 2005-2006


#####################################################################
#####################################################################
DeclareCategory("IsHapSparseMat",IsObject);

DeclareRepresentation(	"IsHapSparseMatRep",
			IsComponentObjectRep,
			["rows",
			 "cols",
			  "characteristic",
			  "mat",
			]);

HapSparseMatFamily:=NewFamily(	"HapSparseMatFamily",
				IsHapSparseMat,
				IsHapSparseMat);

HapSparseMat:=NewType(HapSparseMatFamily,IsHapSparseMatRep);

InstallMethod( ViewObj,
"for HapSparseMat",
 [IsHapSparseMat],
 function(M)
 Print("Sparse matrix with ", M!.rows, " rows and ", M!.cols, " columns in characteristic ", M!.characteristic,"\n");
 end);

InstallMethod( PrintObj,
"for HapSparseMat",
 [IsHapSparseMat],
 function(M)
 Print("Sparse matrix with ", M!.rows, " rows and ", M!.cols, " columns in characteristic ", M!.characteristic,"\n");
 end);

#####################################################################
#####################################################################


