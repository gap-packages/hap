


InstallGlobalFunction("IsPNormal", function( G, p)
#########################################################################
##  Check if the group G is p-normal for the prime p.                  ##
## Zassenhaus defines a finite group to be p-normal if the center of   ##
## one of its Sylow p-groups is the center of every Sylow p-group in   ##
## which it is contained.                                              ##
#########################################################################
local SylowCenter, SylowGroups, N, isPnormal, k, j;
isPnormal := true;

  SylowGroups := ConjugateSubgroups( G, SylowSubgroup( G, p));
  N := Size( SylowGroups);

  for j in [1..N] do

	SylowCenter := Center( SylowGroups[j]);

	for k in [1..N] do
 	  if k <> j then

		if IsSubset(SylowGroups[k], SylowCenter) then

		  if SylowCenter = Center( SylowGroups[k]) then
			; ## the condition is fulfilled for this pair ##
		  else
			isPnormal := false;
#			Print("The group ",G," is not p-normal.");
		  fi;
		fi;
	  fi;
	od;
  od;
  return isPnormal;
end);

InstallGlobalFunction("DisplayAvailableCellComplexes", function()
Exec( Concatenation("ls ",DirectoriesPackageLibrary("HAP")[1]![1], "Perturbations/Gcomplexes/"));
end);

InstallGlobalFunction( "GetTorsionSubcomplex", function(C, p)
#####################################################################
## Here, p is the prime for which to take the torsion subcomplex.  ##
## We extract the cells the stabilizer of which contains p-torsion.##
#####################################################################
local vcd, stabilizerCardinalities, celldata, data, torsionPosition,
torsionCells, numberOfTorsionCells, n, j, returnedData, warned, groupname, admissibilityCheck, x, i, b, tmpCell, cell, boundary, groupName;

    admissibilityCheck := function(celldata)
    #########################################################
    ## A cell complex is admissible in the sense of Brown, ##
    ## if each cell stabilizer fixes its cell pointwise.   ##
    ## Additionally,				       ##
    ## we gather the cardinalities of the stabilizers.     ##
    #########################################################
    local stabilizerCardinalities, G, card, n, j, R, vcd, warned;
       warned := false;
       stabilizerCardinalities := [];
       vcd := Length(celldata)-1;

       for n in [0..vcd] do
	    stabilizerCardinalities[n+1] := [];
	    for j in [1..Length(celldata[n+1])] do
	       G :=   celldata[n+1][j]!.TheMatrixStab;
	       if IsFinite(G) then
	          card := Order(G);
	          stabilizerCardinalities[n+1][j] := card;
	          ## *** Now we have to compare              *** ##
	          ## *** with the order of "TheRotSubgroup"  *** ##
	          R := celldata[n+1][j]!.TheRotSubgroup;
	          if card > Order(R) and warned = false then
		    Print("****Warning: cell complex not admissible ",
			    "in the sense of Brown!****\n",
		    " Torsion subcomplex reduction requires cell subdivision.\n");
		    warned := true;
	          fi;
	       fi;
	    od;
       od;
       return [stabilizerCardinalities, warned];
   end;

   # Case 1: the input is a group name
   if IsString(C) then
   groupName:=C;
   groupname := Filtered( C, function ( x )
            return not (x = '(' or x = ')' or x = ',' or x = '[' or x = ']');
   end );
   Read(Concatenation( 	DirectoriesPackageLibrary("HAP")[1]![1],
			"Perturbations/Gcomplexes/",groupname));
   celldata := StructuralCopy(HAP_GCOMPLEX_LIST);

   # Case 2: the input is a variable
   else
       groupName:=fail;
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
   vcd := Length(celldata) -1;
#   Print("Extracting the ",p,"-torsion subcomplex of the ",
#		vcd,"-dimensional ",groupName,"-cell complex ... \n");
   returnedData := admissibilityCheck(celldata);
   stabilizerCardinalities := returnedData[1];
   warned := returnedData[2];
   torsionCells := [];
   numberOfTorsionCells := [];
   for n in [0..vcd] do
	torsionCells[n+1] := [];
	numberOfTorsionCells[n+1] := 0;
	for j in [1..Length(celldata[n+1])] do
	   ## Check if the stabilizer contains p-torsion ##
	   if stabilizerCardinalities[n+1][j] mod p = 0 then
#		Print("Extracted ",n,"-cell numero ",j,
#			" of stabilizer cardinality ",
#			stabilizerCardinalities[n+1][j],".\n");
		numberOfTorsionCells[n+1]
			:= numberOfTorsionCells[n+1]+1;
	        torsionCells[n+1][numberOfTorsionCells[n+1]]
			:=[n, j];
	   fi;
	od;
   od;
#   return
#     [torsionCells, numberOfTorsionCells, celldata, stabilizerCardinalities, warned];
  data:=[];
#  Print("torsionCells",torsionCells,"\n");
  for i in [1..Length(torsionCells)] do
      data[i]:=[];
      for x in torsionCells[i] do
          Add(data[i],celldata[i][x[2]]);
      od;
  od;
for j in [2..Size(data)] do
  for i in [1..Size(data[j])] do
      tmpCell:=StructuralCopy(data[j][i]!.BoundaryImage);
      b:=List(tmpCell!.ListIFace,w->Position(torsionCells[j-1], [j-2,w]));
      tmpCell!.ListIFace:=b;
      data[j][i]!.BoundaryImage:=tmpCell;
  od;
od;
  torsionPosition:=torsionCells;
  torsionCells:=[];
  for i in [1..Size(data)] do
     torsionCells[i]:=[];
     for j in [1..Size(data[i])] do
         torsionCells[i][j]:=[i-1,j];
     od;
  od;

  return Objectify(HapTorsionSubcomplex,
            rec(

            torsion:=p,
            groupname:=groupName,
            torsionCells:=torsionCells,
            torsionPosition:=torsionPosition,
            originalData:=celldata,
            celldata:= data,
            numberOfTorsionCells:= numberOfTorsionCells,
            stabilizerCardinalities:= stabilizerCardinalities,
            warned:= warned ));
end);

DeclareGlobalFunction("GetTorsionPowerSubcomplex");
InstallGlobalFunction( "GetTorsionPowerSubcomplex", function(C, p, rk)
#####################################################################
## Here, p is the prime for which to take the torsion subcomplex.  ##
## We extract the cells the stabilizer of which contains p-torsion.##
#####################################################################
local vcd, stabilizerCardinalities, celldata, data, torsionPosition,
torsionCells, numberOfTorsionCells, n, j, returnedData, warned, groupname, admissibilityCheck, x, i, b, tmpCell, cell, boundary, groupName, elementaryAbelianID, G, H, subgroupTypesOfG, numberOfSubgroups;

    admissibilityCheck := function(celldata)
    #########################################################
    ## A cell complex is admissible in the sense of Brown, ##
    ## if each cell stabilizer fixes its cell pointwise.   ##
    ## Additionally,				       ##
    ## we gather the cardinalities of the stabilizers.     ##
    #########################################################
    local stabilizerCardinalities, G, card, n, j, R, vcd, warned;
       warned := false;
       stabilizerCardinalities := [];
       vcd := Length(celldata)-1;

       for n in [0..vcd] do
	    stabilizerCardinalities[n+1] := [];
	    for j in [1..Length(celldata[n+1])] do
	       G :=   celldata[n+1][j]!.TheMatrixStab;
	       if IsFinite(G) then
	          card := Order(G);
	          stabilizerCardinalities[n+1][j] := card;
	          ## *** Now we have to compare              *** ##
	          ## *** with the order of "TheRotSubgroup"  *** ##
	          R := celldata[n+1][j]!.TheRotSubgroup;
	          if card > Order(R) and warned = false then
		    Print("****Warning: cell complex not admissible ",
			    "in the sense of Brown!****\n",
		    " Torsion subcomplex reduction requires cell subdivision.\n");
		    warned := true;
	          fi;
	       fi;
	    od;
       od;
       return [stabilizerCardinalities, warned];
   end;

   # Case 1: the input is a group name
   if IsString(C) then
   groupName:=C;
   groupname := Filtered( C, function ( x )
            return not (x = '(' or x = ')' or x = ',' or x = '[' or x = ']');
   end );
   Read(Concatenation( 	DirectoriesPackageLibrary("HAP")[1]![1],
			"Perturbations/Gcomplexes/",groupname));
   celldata := StructuralCopy(HAP_GCOMPLEX_LIST);

   # Case 2: the input is a variable
   else
       groupName:=fail;
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
   vcd := Length(celldata) -1;
#   Print("Extracting the ",p,"-torsion subcomplex of the ",
#		vcd,"-dimensional ",groupName,"-cell complex ... \n");
   elementaryAbelianID := IdSmallGroup(DirectProduct(List([1..rk], i -> CyclicGroup(p))));
   returnedData := admissibilityCheck(celldata);
   stabilizerCardinalities := returnedData[1];
   warned := returnedData[2];
   torsionCells := [];
   numberOfTorsionCells := [];
   for n in [0..vcd] do
	torsionCells[n+1] := [];
	numberOfTorsionCells[n+1] := 0;
	for j in [1..Length(celldata[n+1])] do
	   ## Check if the stabilizer contains p-torsion ##
	   if stabilizerCardinalities[n+1][j] mod p^rk = 0 then

		G := C!.stabilizer(n,j);
		subgroupTypesOfG := []; numberOfSubgroups := 0;
		for H in ConjugacyClassesSubgroups(G) do
			numberOfSubgroups := numberOfSubgroups+1;
			subgroupTypesOfG[numberOfSubgroups] := IdSmallGroup(Representative(H));
		od;
		if IsInt(Position(subgroupTypesOfG,elementaryAbelianID)) then

#			Print("Extracted ",n,"-cell numero ",j,
#			" of stabilizer cardinality ",
#			stabilizerCardinalities[n+1][j],".\n");
			numberOfTorsionCells[n+1]
				:= numberOfTorsionCells[n+1]+1;
	        	torsionCells[n+1][numberOfTorsionCells[n+1]]
				:=[n, j];
		fi;
	   fi;
	od;
   od;
#   return
#     [torsionCells, numberOfTorsionCells, celldata, stabilizerCardinalities, warned];
  data:=[];
  for i in [1..Length(torsionCells)] do
      data[i]:=[];
      for x in torsionCells[i] do
          Add(data[i],celldata[i][x[2]]);
      od;
  od;

for j in [2..Size(data)] do
  for i in [1..Size(data[j])] do
      tmpCell:=StructuralCopy(data[j][i]!.BoundaryImage);
      b:=List(tmpCell!.ListIFace,w->Position(torsionCells[j-1], [j-2,w]));
      tmpCell!.ListIFace:=b;
      data[j][i]!.BoundaryImage:=tmpCell;
  od;
od;
  torsionPosition:=torsionCells;
  torsionCells:=[];
  for i in [1..Size(data)] do
     torsionCells[i]:=[];
     for j in [1..Size(data[i])] do
         torsionCells[i][j]:=[i-1,j];
     od;
  od;

  return Objectify(HapTorsionSubcomplex,
            rec(

            torsion:=p,
            groupname:=groupName,
            torsionCells:=torsionCells,
            torsionPosition:=torsionPosition,
            originalData:=celldata,
            celldata:= data,
            numberOfTorsionCells:= numberOfTorsionCells,
            stabilizerCardinalities:= stabilizerCardinalities,
            warned:= warned ));
end);


InstallGlobalFunction("TorsionSubcomplex", function(groupName, p)
############################################
local torsionCells, numberOfTorsionCells, celldata, sortedData, warned, computeIncidenceMatrix;


    computeIncidenceMatrix := function(torsionCells, numberOfTorsionCells, celldata)
    ########################################################
    ## The incidence matrix of the 1-skeleton of the      ##
    ## p-torsion subcomplex is returned by this function. ##
    ########################################################
    local incidenceMatrix, j, k, q, endpoints, inverseIndex, ORIGIN, END;
       incidenceMatrix := [];
       for j in [1..numberOfTorsionCells[1]] do
        	incidenceMatrix[j] := [];
         	for k in [1..numberOfTorsionCells[1]] do
        		incidenceMatrix[j][k] := 0;
        	od;
       od;
#       Print("The edges in the quotient by the action on the torsion subcomplex are: ");
       ## Record the indices in the torsion subcomplex of vertices in  ##
       ## the cell complex, in order to assign the endpoints of an edge##
       inverseIndex := [];
       for q in [1..Length(celldata[1])] do
        	inverseIndex[q] := "error: not a p-torsion vertex";
       od;
       for j in [1..numberOfTorsionCells[1]] do
        	inverseIndex[torsionCells[1][j][2]] := j;
       od;
       for q in [1..numberOfTorsionCells[2]] do
        	endpoints := celldata[2]
			[torsionCells[2][q][2]]!.BoundaryImage.ListIFace;
	    ORIGIN := inverseIndex[endpoints[1]];
	    END    := inverseIndex[endpoints[2]];
#	    Print([ORIGIN,END]);
	    incidenceMatrix[ORIGIN][END] := incidenceMatrix[ORIGIN][END] +1;
	    ## If we do not want to have the incidence matrix symmetric, ##
	    ## then deactivate the following line. ##
	    incidenceMatrix[END][ORIGIN] := incidenceMatrix[END][ORIGIN] +1;
       od;
#       Print(".\n");
       return incidenceMatrix;
    end;


	sortedData := GetTorsionSubcomplex(groupName, p);
	torsionCells := sortedData!.torsionCells;
	numberOfTorsionCells := sortedData!. numberOfTorsionCells;
	celldata := sortedData!.celldata;
	warned := sortedData!.warned;
	if IsPrime(p) then
	  sortedData := [computeIncidenceMatrix(torsionCells, numberOfTorsionCells, celldata), warned];
	else
	  sortedData := "The number p must be prime in order to compute the incidence matrix.";
	fi;
   return
     sortedData;
end);



InstallGlobalFunction("VisualizeTorsionSkeleton", function(groupName, p)
##################################################
local incidenceMatrix, graphData, returnedData, warned;
    returnedData := TorsionSubcomplex(groupName, p);
    if IsPrime(p) then
	incidenceMatrix := returnedData[1];
	warned := returnedData[2];
	if warned = false then
#	    Print("The quotient by the action on the 1-skeleton of the p-torsion #subcomplex is ");
	else
#	    Print("As the cell stabilizers do not fix the cells pointwise, \n",
#		"we do not obtain the quotient by the action on the 1-skeleton of the p-#torsion subcomplex.\n",
#		" We obtain just some graph, and the latter is ");
	fi;
	if Size(incidenceMatrix) = 0 then
		Print("empty.\n");
	fi;
	if Size(incidenceMatrix) > 0 then
#		Print("displayed in a separate window on screen now.\n");
		graphData := StructuralCopy(IncidenceMatrixToGraph(incidenceMatrix));
		GraphDisplay(graphData);
	fi;
   else
#Print("The number p must be prime in order to compute the incidence matrix.\n");
   fi;
end);




InstallGlobalFunction("ReduceTorsionSubcomplex", function(arg) ## The input arg consists of groupname, p, and optionally a third argument rk,
#############################I######################
local torsionCells, numberOfTorsionCells, celldata, sortedData, stabilizerCardinalities, fusionCandidates, reducedTorsionCells, groupname,
terminaljMinus1cells, warned, extractPmultipleTorsionCells, cellcomplex,
bifurcationFreeCells, pairsIntersection, groupName, p, rk,
getIdentifier, checkStabilizerConjugacy, mergeCells, printIsotropyGroups, ReduceModP, getTerminalFacets, cutOffCells,
terminalFacets, data, i, x, tmpCell,
N, J, n, j, G, B, b;

groupName:=arg[1];
p:=arg[2];

extractPmultipleTorsionCells := function(torsionCells, numberOfTorsionCells, celldata,stabilizerCardinalities, p)
##################################################################
## Make lists of n-cells of each stabilizer cardinality         ##
## in the p-torsion subcomplex.                                 ##
##################################################################
local vcd, pMultipleTorsionCells, numberOfPmultipleTorsionCells, Pmultiples, PmultipleNumero, m, n, cell, j, getCardinals;


getCardinals := function(stabilizerCardinalities, p, n)
############################################################
## Extract the cardinalities of n-cell stabilisers        ##
## containing p-torsion.                                  ##
############################################################
local S, m, Cardinals, counter;
    Cardinals := []; counter := 0;
    S := stabilizerCardinalities[n+1];
    for m in S do
	if m mod p = 0 then
	  if counter > 0 then
	    if m in Cardinals then ;
	    else
		  counter := counter +1;
		  Cardinals[counter] := m;
	    fi;
	  else counter := 1;
	       Cardinals[counter] := m;
	  fi;
	fi;
    od;
    return Cardinals;
end;
#######################################


  vcd := Length(celldata) -1;
  pMultipleTorsionCells := [];
  numberOfPmultipleTorsionCells := [];
  Pmultiples := [];

  for n in [0..vcd] do

    pMultipleTorsionCells[n+1] := [];
    numberOfPmultipleTorsionCells[n+1] := [];
    Pmultiples[n+1] := getCardinals(stabilizerCardinalities, p, n);
    PmultipleNumero := 0;

    for m in Pmultiples[n+1] do

      PmultipleNumero := PmultipleNumero +1;
      numberOfPmultipleTorsionCells[n+1][PmultipleNumero] := 0;
      pMultipleTorsionCells[n+1][PmultipleNumero] := [];

      for cell in torsionCells[n+1] do

 	## Check if the stabilizer contains p-torsion of multiple m ##


	if stabilizerCardinalities[n+1][cell[2]] = m then

	  numberOfPmultipleTorsionCells[n+1] :=
		numberOfPmultipleTorsionCells[n+1] +1;
 	  pMultipleTorsionCells[n+1][PmultipleNumero]
		[numberOfPmultipleTorsionCells[n+1][PmultipleNumero]] :=
		cell;
	fi;
      od;
#      Print("There are ",numberOfPmultipleTorsionCells[n+1][PmultipleNumero],
#		" orbits of ",n,"-cells in the ",p,"-torsion subcomplex",
#		", the stabilizers of which are of cardinality ",m,
#		", namely the ones numero "
#	);
      for j in [1..Length(pMultipleTorsionCells[n+1][PmultipleNumero])-1] do
#	    Print(pMultipleTorsionCells[n+1][PmultipleNumero][j][2],", ");
      od;
#      Print(pMultipleTorsionCells[n+1][PmultipleNumero][
#	Length(pMultipleTorsionCells[n+1][PmultipleNumero])][2],".\n");
    od;
  od;
  return [pMultipleTorsionCells, Pmultiples];
end;

pairsIntersection := function(sortedData, celldata)
#########################################################################
## We establish a preliminary list of cell triples as candidates for   ##
## fusions, consisting of a pair of n-cells that admits precisely      ##
## one adjacent (n-1)-cell in common.                                  ##
## This pair and their common boundary cell constitute our cell triple.##
#########################################################################
local pMultipleTorsionCells, Pmultiples, n, j, f, s, vcd, fusionCandidates,
numberOfFusionCandidates, firstCell, secondCell, commonBoundary;
  pMultipleTorsionCells := sortedData[1];
  Pmultiples := sortedData[2];
  vcd := Length(celldata) -1;
  fusionCandidates := [];
  numberOfFusionCandidates := 0;

  for n in [0..vcd] do
    for j in [1..Length(Pmultiples[n+1])] do
      for f in [1..Length(pMultipleTorsionCells[n+1][j])] do
	firstCell := pMultipleTorsionCells[n+1][j][f];
        for s in [f+1..Length(pMultipleTorsionCells[n+1][j])] do
	  secondCell := pMultipleTorsionCells[n+1][j][s];
	  commonBoundary :=
		Set(celldata[n+1][firstCell[2]]!.BoundaryImage!.ListIFace);
	  IntersectSet( commonBoundary,
	    Set(celldata[n+1][secondCell[2]]!.BoundaryImage!.ListIFace)
	  );

	  if Size(commonBoundary)=1  then
#	    Print("The boundary of ",n,"-cell numero ",firstCell[2]," is ",
#	      celldata[n+1][firstCell[2]]!.BoundaryImage!.ListIFace,
#	      " and that of ",n,"-cell numero ",secondCell[2]," is ",
#	      celldata[n+1][secondCell[2]]!.BoundaryImage!.ListIFace,
#	      " so they have ",n-1,"-cell numero ",commonBoundary[1],
#		" in common.\n");
	    numberOfFusionCandidates := numberOfFusionCandidates +1;
	    fusionCandidates[numberOfFusionCandidates] :=
	      [firstCell, secondCell, commonBoundary[1]];
              return fusionCandidates;
	  fi;
        od;
      od;
    od;
  od;
  return fusionCandidates;
end;


bifurcationFreeCells := function(fusionCandidates, torsionCells, celldata)
##########################################################################
## We keep only those cell triples as candidates for a fusion,          ##
## that are not "bifurcated" in the p-torsion subcomplex :              ##
## Triples such that only n-cells adjacent to the (n-1)-cell            ##
## of the triple are the two n-cells of the triple.                     ##
##########################################################################
local firstCell, secondCell, boundaryCell, cellTriple, n, cell, torsionCell,
remainingFusionCandidates;
remainingFusionCandidates := Set(StructuralCopy(fusionCandidates));

  for cellTriple in fusionCandidates do

    firstCell := cellTriple[1];
    secondCell := cellTriple[2];
    n := firstCell[1];
    boundaryCell := cellTriple[3];

      for torsionCell in torsionCells[n+1] do

	cell := celldata[n+1][torsionCell[2]];

	if boundaryCell in cell!.BoundaryImage!.ListIFace then

	  if torsionCell in [firstCell, secondCell] then ;
	  	# the cell triple remains a candidate for fusion.
	  else
#		Print("The p-torsion subcomplex is bifurcated at the ",
#			n-1,"-cell numero ",boundaryCell,".\n");
	  	RemoveSet(remainingFusionCandidates, cellTriple);
	  fi;
	fi;
      od;
  od;
  return remainingFusionCandidates;
end;


getIdentifier := function( cell, boundaryCell, celldata)
##############################################################
## Retrieve the group element sending the cell to           ##
## a representative that is adjacent to the boundary cell.  ##
##############################################################
local n, boundaryIndex, j, identifiers, g;
    n := cell[1];
    identifiers := celldata[n+1][cell[2]]!.BoundaryImage!.ListElt;

    for j in [1..Length(identifiers)] do

      if celldata[n+1][cell[2]]!.BoundaryImage!.ListIFace[j]
	 = boundaryCell then

		boundaryIndex := j;
      fi;
    od;
    g := celldata[n+1][cell[2]]!.BoundaryImage!.ListElt[boundaryIndex];
  return g;
end;


checkStabilizerConjugacy := function(fusionCandidates, celldata, p)
##################################################################
## Eliminate cell triples from the fusion candidates list,      ##
## for which the stabilizers of the n-cells are not conjugate.  ##
## Therefore, we first conjugate the stabilizers of the n-cells ##
## into stabilizers of representatives                          ##
## adjacent to the (n-1)-cell.                                  ##
## Then we check if we can conjugate them within the stabilizer ##
## of the (n-1)-cell.     					##
##################################################################
local firstCell, secondCell, boundaryCell, firstStabilizer, secondStabilizer, boundaryStabilizer, cellTriple, n, g, validated, remainingFusionCandidates, controllingGroupID;
remainingFusionCandidates := StructuralCopy(fusionCandidates);

  for cellTriple in fusionCandidates do

    firstCell := cellTriple[1];
    secondCell := cellTriple[2];
    n := firstCell[1];
    boundaryCell := cellTriple[3];
    firstStabilizer := celldata[n+1][firstCell[2]]!.TheMatrixStab;
    secondStabilizer := celldata[n+1][secondCell[2]]!.TheMatrixStab;
    boundaryStabilizer := celldata[n][boundaryCell]!.TheMatrixStab;
    g := getIdentifier( firstCell, boundaryCell, celldata);
    firstStabilizer := ConjugateGroup(firstStabilizer, g);
    g := getIdentifier( secondCell, boundaryCell, celldata);
    secondStabilizer := ConjugateGroup(secondStabilizer, g);
    controllingGroupID := IdSmallGroup(ReduceModP(firstStabilizer,p));

    if controllingGroupID =IdSmallGroup(ReduceModP(secondStabilizer,p))
   and controllingGroupID =IdSmallGroup(ReduceModP(boundaryStabilizer,p))
       then ;
        else RemoveSet(remainingFusionCandidates, cellTriple);
    fi;
  od;
  return remainingFusionCandidates;
end;


mergeCells := function( celldata, fusionCandidates, torsionCells)
####################################################################
## Merge the cells in that have been definitively retained        ##
## as fusion candidates. As proven, this does not change 	  ##
## the p-primary equivariant Farrell cohomology.		  ##
####################################################################
local firstCell, secondCell, boundaryCell, cellTriple, n, b, t1, t2, BI1, BI2, x,y,t01,t02,
      i, j, reducedTorsionCells, vcd, mergedBoundary, data, p0, p1, p2, tmpCell, U, V, X;
  mergedBoundary := [];
  reducedTorsionCells := StructuralCopy( torsionCells);
  vcd := Length(reducedTorsionCells) -1;
  for n in [0..vcd] do
	reducedTorsionCells[n+1] := Set(reducedTorsionCells[n+1]);
	mergedBoundary[n+1] := [];
  od;

  # Extract the data of torsionCells from original cell data
  data:=[];
  for i in [1..Length(reducedTorsionCells)] do
      data[i]:=[];
      for x in reducedTorsionCells[i] do
          Add(data[i],celldata[i][x[2]]);
      od;
  od;
#Print("reducedTorsionCells", reducedTorsionCells,"\n");
#Print("data= ",data,"\n");
  for cellTriple in fusionCandidates do

    firstCell := cellTriple[1];
    secondCell := cellTriple[2];
    n := firstCell[1];
    boundaryCell := cellTriple[3];
    p0:=Position(reducedTorsionCells[n], [n-1,boundaryCell]);
    p1:=Position(reducedTorsionCells[n+1], [n,firstCell[2]]);
    p2:=Position(reducedTorsionCells[n+1], [n,secondCell[2]]);

#    Print("The ",n-1,"-cell numero ",boundaryCell,
#	" is amalgamated into the merging ",n,"-cells numero ",
#	firstCell[2]," and ",secondCell[2],".\n");
    RemoveSet(reducedTorsionCells[n], [n-1,boundaryCell]);
    Remove(data[n],p0);
    mergedBoundary[n+1][firstCell[2]] := Set(Concatenation(
	celldata[n+1][firstCell[2]]!.BoundaryImage!.ListIFace,
	celldata[n+1][secondCell[2]]!.BoundaryImage!.ListIFace
    ));

    RemoveSet(mergedBoundary[n+1][firstCell[2]], boundaryCell);
    t1:=Position(data[n+1][p1]!. BoundaryImage!.ListIFace,mergedBoundary[n+1][firstCell[2]][1]);
    if t1=fail then
        t1:=Position(data[n+1][p1]!. BoundaryImage!.ListIFace,mergedBoundary[n+1][firstCell[2]][2]);
        t01:=3-t1;
        t2:=Position(data[n+1][p2]!. BoundaryImage!.ListIFace,mergedBoundary[n+1][firstCell[2]][1]);
        t02:=3-t2;
    else t2:=Position(data[n+1][p2]!. BoundaryImage!.ListIFace,mergedBoundary[n+1][firstCell[2]][2]);
    fi;
    t01:=3-t1;
    t02:=3-t2;
    BI1:=StructuralCopy(celldata[n+1][firstCell[2]]!.BoundaryImage);
    BI2:=StructuralCopy(celldata[n+1][secondCell[2]]!.BoundaryImage);


    x:=BI2!.ListElt[t02]*BI2!.ListElt[t01]^-1;
    U:=ConjugateGroup(data[n+1][p2]!.TheMatrixStab,(x)^-1);
    V:=data[n+1][p1]!.TheMatrixStab;
    X:=ConjugateGroup(celldata[n][BI1!.ListIFace[t01]]!.TheMatrixStab,BI1!.ListElt[t01]);
    y:=RepresentativeAction(X,U,V);
    tmpCell:=rec(ListIFace:=[BI1!.ListIFace[t1],BI2!.ListIFace[t2]],
                        ListSign:=BI1!.ListSign,
                        ListElt:=[BI1!.ListElt[t1],y*x*BI2!.ListElt[t2]]
    );
    data[n+1][p1]!.BoundaryImage:=StructuralCopy(tmpCell);
    Remove(data[n+1],p2);
  od;
  for i in [1..Size(data[2])] do
      tmpCell:=StructuralCopy(data[2][i]!.BoundaryImage);
      b:=List(tmpCell!.ListIFace,w->Position(reducedTorsionCells[1], [0,w]));
      tmpCell!.ListIFace:=b;
      data[2][i]!.BoundaryImage:=tmpCell;
  od;
  reducedTorsionCells:=[];
  for i in [1..Size(data)] do
     reducedTorsionCells[i]:=[];
     for j in [1..Size(data[i])] do
         reducedTorsionCells[i][j]:=[i-1,j];
     od;
  od;
  return [reducedTorsionCells, data];
end;


printIsotropyGroups := function(reducedTorsionCells, celldata)
#########################################################
## Return the stabilizers of the cells in the          ##
## reduced torsion subcomplex.                         ##
#########################################################
local G, n, k, j, tcd;
   tcd := Length(reducedTorsionCells)-1;
   Print("The following cells and stabilizers represent ",
		"the reduced torsion subcomplex.\n");

   for n in [0..tcd] do
     for k in [1..Length(reducedTorsionCells[n+1])] do
	j := reducedTorsionCells[n+1][k][2];
	G :=   celldata[n+1][j]!.TheMatrixStab;
	if IsFinite(G) then
		Print(n,"-cell number ",j," has stabilizer ",G,"\n");
		Print("of Abelianization ",GroupHomology(G,1),".\n");
	fi;
     od;
   od;
end;


###################################################################
#                       CUT OFF CELLS PROCESS                     #
#                                                                 #
###################################################################
cutOffCells := function(terminalFacets, torsionCells, celldata, n)
####################################################################
## Cut off the n-cells that have a terminal facet which has been  ##
## spotted to have the same p-torsion controller. As proven, this ##
## does not change the p-primary equivariant Farrell cohomology.  ##
####################################################################
local facet, nCell, reducedTorsionCells, cutnCells, cellPair;

  cutnCells := [];
  reducedTorsionCells := StructuralCopy( torsionCells);
  reducedTorsionCells[n+1] := Set(reducedTorsionCells[n+1]);
  reducedTorsionCells[n] := Set(reducedTorsionCells[n]);

  for cellPair in terminalFacets do

    facet := cellPair[1];
    nCell := cellPair[2];
    if not nCell in cutnCells then
#    	Print("The terminal ",n-1,"-cell numero ",facet,
#	" is cut against ",n,"-cell numero ",nCell,".\n");
    	RemoveSet(reducedTorsionCells[n], [n-1,facet]);
    	RemoveSet(reducedTorsionCells[n+1], [n,nCell]);
	cutnCells := Concatenation(cutnCells,[nCell]);
    fi;
  od;
#  Print("Length of reducedTorsionCells after cut off ",n,"-cells",List(reducedTorsionCells,i->Length(i)),"\n");
  return reducedTorsionCells;
end;


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

getTerminalFacets := function(reducedTorsionCells, celldata, p, n)
##############################################################################
## We single out the "terminal" (n-1)-facets of the p-torsion subcomplex :  ##
## Those with exactly one adjacent n-cell.                                  ##
##############################################################################
local examinedncell, examinedFacet, examinedncellData, numberOfAdjacencies,
	adjacentncell, reductionCandidates, numberOfTerminalFacets, G, IdGred, H, IdHred,
  nplus1cell;
  reductionCandidates:= []; numberOfTerminalFacets := 0;

  for examinedFacet in reducedTorsionCells[n] do

    examinedFacet := examinedFacet[2];
    numberOfAdjacencies := 0;

    for examinedncell in reducedTorsionCells[n+1] do
	examinedncellData := celldata[n+1][examinedncell[2]];

	if examinedFacet in examinedncellData!.BoundaryImage!.ListIFace then

	  numberOfAdjacencies := numberOfAdjacencies +1;
	  adjacentncell := examinedncell[2];
	fi;
    od;
    if numberOfAdjacencies = 1 then
        nplus1cell:=0;
        if IsBound(reducedTorsionCells[n+2]) then
         for examinedncell in reducedTorsionCells[n+2] do
             examinedncellData := celldata[n+2][examinedncell[2]];
             if adjacentncell in examinedncellData!.BoundaryImage!.ListIFace then
                 nplus1cell := nplus1cell +1;
             fi;
          od;
        fi;
        if nplus1cell=0 then
	        G := celldata[n+1][adjacentncell]!.TheMatrixStab;;
	        IdGred := IdSmallGroup(ReduceModP(G,p));
	        H := celldata[n][examinedFacet]!.TheMatrixStab;;
	        IdHred := IdSmallGroup(ReduceModP(H,p));
#	     Print("The ",n-1,"-facet number ",examinedFacet," is terminal.",
#		" Its stabilizer ",IdSmallGroup(H)," is controlled by ", IdHred,
#		". It is adjacent to the ",n,"-cell number ",adjacentncell,
#		" of stabilizer ",IdSmallGroup(G)," controlled by ",IdGred,".\n");
	         if IdGred = IdHred then
#	     Print("Reduce it against the adjacent ",n,"-cell number ",adjacentncell,
#		" of stabilizer ",IdSmallGroup(G)," controlled by ",IdGred,".\n");
		           numberOfTerminalFacets := numberOfTerminalFacets +1;
		           reductionCandidates[numberOfTerminalFacets] :=
				       [examinedFacet, adjacentncell];
	          fi;
         fi;
      fi;
  od;
#Print("Reduction Candidates of dimension",n,"    ",reductionCandidates,"\n");
  return reductionCandidates;
end;
#############################
   if Length(arg)=2 then
     sortedData := GetTorsionSubcomplex(groupName, p);
   else
	rk := arg[3];
        sortedData := GetTorsionPowerSubcomplex(groupName, p, rk);
   fi;
   torsionCells := sortedData!.torsionCells;
   numberOfTorsionCells := sortedData!. numberOfTorsionCells;
   celldata := sortedData!.celldata;
   stabilizerCardinalities := sortedData!. stabilizerCardinalities;
   warned := sortedData!.warned;
   groupname:=sortedData!.groupname;
#   originalData:=sortedData!.originalData;
#Print("Before cutting off: ", torsionCells,"\n");
#   torsionCells := reducedTorsionCells;
#   numberOfTorsionCells := sortedData[2];
#   celldata := sortedData[3];
#   stabilizerCardinalities := sortedData[4];
#   warned := sortedData[5];
N:=1;
while IsBound(torsionCells[N]) and (not IsEmpty(torsionCells[N])) do
	N:=N+1;
od;
N:=N-2;
j:=N;

reducedTorsionCells :=StructuralCopy(torsionCells);
while j>0 do
#Print("Dimension ",j,"\n");
terminalFacets := getTerminalFacets(reducedTorsionCells, celldata, p,j);

while not terminalFacets = [] do
	reducedTorsionCells := cutOffCells(terminalFacets, reducedTorsionCells, celldata, j);
	terminalFacets := getTerminalFacets(reducedTorsionCells, celldata, p,j);
od;

j:=j-1;
od;

  # Extract the data of torsionCells from original cell data
  data:=[];
  for i in [1..Length(reducedTorsionCells)] do
      data[i]:=[];
      for x in reducedTorsionCells[i] do
          Add(data[i],celldata[i][x[2]]);
      od;
  od;
#Print("reducedTorsionCells",reducedTorsionCells,"\n");
  for j in [2..Length(data)] do
    for i in [1..Size(data[j])] do
      tmpCell:=StructuralCopy(data[j][i]!.BoundaryImage);
      b:=List(tmpCell!.ListIFace,w->Position(reducedTorsionCells[j-1], [j-2,w]));
      tmpCell!.ListIFace:=b;
      data[j][i]!.BoundaryImage:=tmpCell;
    od;
  od;
  reducedTorsionCells:=[];
  for i in [1..Size(data)] do
     reducedTorsionCells[i]:=[];
     for j in [1..Size(data[i])] do
         reducedTorsionCells[i][j]:=[i-1,j];
     od;
  od;
#  originalData:=celldata;
  celldata:=StructuralCopy(data);
#########################End of cutting off cells#################################

 torsionCells := reducedTorsionCells;

 if warned = false then
   sortedData := extractPmultipleTorsionCells(torsionCells,
		numberOfTorsionCells, celldata,stabilizerCardinalities, p);
#Print("torsionCells  ", torsionCells,"\n");
#Print("stabilizerCardinalities  ", stabilizerCardinalities,"\n");
#Print("stabilizerCardinalities  ", stabilizerCardinalities,"\n");

   fusionCandidates := pairsIntersection(sortedData, celldata);
   fusionCandidates := bifurcationFreeCells(fusionCandidates, torsionCells,
						 celldata);
   fusionCandidates := checkStabilizerConjugacy(fusionCandidates, celldata,
						p);
   while not fusionCandidates=[] do

#Print("fusionCandidates  ", fusionCandidates,"\n");
   cellcomplex:= mergeCells( celldata, fusionCandidates,
					torsionCells);
    reducedTorsionCells:=cellcomplex[1];
#Print("reducedTorsionCells  ", reducedTorsionCells,"\n");
    celldata:=cellcomplex[2];
    numberOfTorsionCells:=List([1..Size(reducedTorsionCells)],
             i->Size(reducedTorsionCells[i]));
    stabilizerCardinalities:=List([1..Size(reducedTorsionCells)],
             i->List([1..Size(reducedTorsionCells[i])],j->Order(celldata[i][j]!.TheMatrixStab)));
    sortedData := extractPmultipleTorsionCells(reducedTorsionCells,
		numberOfTorsionCells, celldata,stabilizerCardinalities, p);
    fusionCandidates := pairsIntersection(sortedData, celldata);
   fusionCandidates := bifurcationFreeCells(fusionCandidates, torsionCells,
						 celldata);
   fusionCandidates := checkStabilizerConjugacy(fusionCandidates, celldata,
						p);
#Print("fusionCandidates:  ", fusionCandidates,"\n");
    od;
#   terminaljMinus1cells := getterminaljMinus1cells(reducedTorsionCells, celldata);
#   terminaljMinus1cells := checkGruenSwan( terminaljMinus1cells, celldata, p);
#Print("\n At the following terminal vertices, the adjacent edge can be cut off without #changing the equivariant Farrell cohomology : ",terminaljMinus1cells,"\n");
#   printIsotropyGroups(reducedTorsionCells, celldata);
 fi;

N:=Size(celldata);
while N>0 do
    if  celldata[N]=[] then
        Remove(reducedTorsionCells,N);
        Remove(celldata,N);
    fi;
    N:=N-1;
od;

#####################################################################################


  return Objectify(HapTorsionSubcomplex,
            rec(

            torsion:=p,
            groupname:=groupname,
            torsionCells:=torsionCells,
#            originalData:=data,
            celldata:= celldata ));

end);

#####################################################################################
