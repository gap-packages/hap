
InstallGlobalFunction("ExportHapCellcomplexToDisk", function( C,groupName) 
local GcomplexesPath, exportFile, ExportToHAP, i, j, celldata, torsion, m,n,k, cell, boundary; 


if IsHapTorsionSubcomplex(C) then
celldata:=C!.celldata;
torsion:=C!.torsion;
else
       celldata:=[];
       i:=0;
       while C!.dimension(i) > 0 do
           cell:=[];
           for j in [1..C!.dimension(i)] do
               if not i=0 then 
               boundary:=C!.boundary(i,j);
               Add(cell,rec(TheMatrixStab :=C!.stabilizer(i,j),
                           TheRotSubgroup:=C!.stabilizer(i,j),
                           BoundaryImage :=rec(
                                 ListIFace:=List(boundary,w->AbsInt(w[1])),
                                 ListSign:=List(boundary,w->SignInt(w[1])),
                                 ListElt:=List(boundary,w->C!.elts[w[2]])
                           )
                       )
               );
               else
               Add(cell,rec(TheMatrixStab :=C!.stabilizer(i,j),
                           TheRotSubgroup:=C!.stabilizer(i,j),
                           BoundaryImage :=rec(
                                 ListIFace:=[],
                                 ListSign:=[],
                                 ListElt:=[]
                           )
                       )
               );                   
               fi;
           od;
           Add(celldata,cell);
           i:=i+1;
       od; 
fi;

 GcomplexesPath := Concatenation(DirectoriesPackageLibrary("HAP")[1]![1],"Perturbations/Gcomplexes/");
 exportFile := Concatenation("./", groupName);


ExportToHAP := function(i,j)
local currentCell, currentStab;


  currentCell := celldata[i][j];
  currentStab:=currentCell!.TheMatrixStab;
  if IdSmallGroup(currentStab) = [1,1] then
  	AppendTo(  exportFile, "\n  rec( TheMatrixStab := Group(");
  	AppendTo( exportFile,  Elements(currentStab));
  	AppendTo(  exportFile, "), \n       TheRotSubgroup := Group(");
  	AppendTo( exportFile, Elements(currentStab) );
  else	
  	AppendTo(  exportFile, "\n  rec( TheMatrixStab := Group(");
  	AppendTo( exportFile, GeneratorsOfGroup(currentStab) );
  	AppendTo(  exportFile, "), \n       TheRotSubgroup := Group(");
  	AppendTo( exportFile, GeneratorsOfGroup(currentStab) );
  fi;
       AppendTo( exportFile,  "), \n       BoundaryImage := rec("); 
       AppendTo(  exportFile,  "\n          ListIFace := ");
       AppendTo( exportFile, currentCell!.BoundaryImage!.ListIFace );
       AppendTo( exportFile, ",");
       AppendTo(  exportFile,  "\n          ListSign := ");
       AppendTo( exportFile, currentCell!.BoundaryImage!.ListSign );
       AppendTo( exportFile, ",");
       AppendTo(  exportFile,  "\n          ListElt := ");
       AppendTo( exportFile, currentCell!.BoundaryImage!.ListElt );
       AppendTo( exportFile, "  ) \n  )");
end;
# end of function ExportToHAP
 
  AppendTo( exportFile, "HAP_GCOMPLEX_SETUP:=[false]; \n");   
#  AppendTo( exportFile, "torsion:=",torsion,"; \n");
  AppendTo( exportFile, "HAP_GCOMPLEX_LIST := [ \n");
  n:=Length(celldata);
  for i in [1..(n-1)] do
        m:=Length(celldata[i]);
  	AppendTo( exportFile, " ["); 

	for j in [1..(m-1)] do
		ExportToHAP(i,j);
		AppendTo( exportFile, ",");
	od;

        ExportToHAP(i,m);
  	AppendTo( exportFile, "\n ], \n"); 
  od;
        m:=Length(celldata[n]);
  	AppendTo( exportFile, " ["); 

	for j in [1..(m-1)] do
		ExportToHAP(n,j);
		AppendTo( exportFile, ",");
	od;

        ExportToHAP(n,m);
  	AppendTo( exportFile, "\n ] \n"); 
  AppendTo( exportFile,  "]; \n");
  
Print("Cell complex written in file ./bin/",exportFile,", you have to move it to \n",GcomplexesPath);
end);
  
# end of program HAPtoDisk


