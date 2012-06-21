############################################################
##                                                        ##
##               TorsionSubcomplex.g                      ##
## HAP subpackage for GAP (Groups Algorithms Programming) ##
##         under the GNU GPL license (v. 3),  2012        ##
##               by Alexander D. Rahm                     ##
##          IRCSET post-doctoral Research fellow          ##
##	   National University of Ireland at Galway       ##
##                                                        ##
############################################################   	

###############################################################################
# HAP subpackage by Alexander D. Rahm for extracting torsion subcomplexes     #
# from cellular complexes with group actions. Version 1.0 of June 19th, 2012. # # The syntax is TorsionSubcomplex(groupName,p);                               #
#  where p is the prime for which to extract the p-torsion subcomplex,        #
# and groupName must be a string holding between \"quotes\" an entry of the   # # HAP library of cellular complexes with group actions,                       # 
# which admits the entries displayed with the following command.              #
# DisplayAvailableCellComplexes();					      #
# The incidence matrix of the 1-skeleton of the p-torsion subcomplex          # # is returned by this function. To visualize this 1-skeleton as a graph,      # # type instead VisualizeTorsionSkeleton(groupName,p).                         #
# If the cell complex is admissible in the sense of Brown,                    # # the Farrell cohomology of the group equals the equivariant Farrell cohomology
# of the reduced torsion subcomplex.  The syntax is 
# ReduceTorsionSubcomplex(groupName,p); 
# and it returns us the cells which are to merge in the torsion subcomplex.   #
###############################################################################

InstallGlobalFunction("IsPnormal", function( G, p)
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
			Print("The group ",G," is not p-normal.");
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

admissibilityCheck := function(celldata)
#########################################################
## A cell complex is admissible in the sense of Brown, ##
## if each cell stabilizer fixes its cell pointwise.   ##
## Additionally,				       ##
## we gather the cardinalities of the stabilizers.     ##
#########################################################
local stabilizerCardinalities, G, card, n, j, R, vcd;
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
	      if card > Order(R) then
		Print("****Warning: cell complex not admissible ",
			"in the sense of Brown!****",
		" Torsion subcomplex reduction requires cell subdivision.");
	      fi;	
	   fi;
	od;
   od;
   return stabilizerCardinalities;
end;	   


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
   Print("The 1-skeleton of the torsion subcomplex is: ");
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
	Print([ORIGIN,END]);
	incidenceMatrix[ORIGIN][END] := incidenceMatrix[ORIGIN][END] +1;
	## If we do not want to have the incidence matrix symmetric, ##
	## then deactivate the following line. ##
	incidenceMatrix[END][ORIGIN] := incidenceMatrix[END][ORIGIN] +1;
   od;
   Print(".\n");
   Print("The incidence matrix of the 1-skeleton of the p-torsion subcomplex is \n");
   return incidenceMatrix;
end;


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


extractPmultipleTorsionCells := function(torsionCells, numberOfTorsionCells, celldata,stabilizerCardinalities, p)
##################################################################
## Make lists of n-cells of each stabilizer cardinality         ##
## in the p-torsion subcomplex.                                 ##
##################################################################
local vcd, pMultipleTorsionCells, numberOfPmultipleTorsionCells, Pmultiples, PmultipleNumero, m, n, cell, j;

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
      Print("There are ",numberOfPmultipleTorsionCells[n+1][PmultipleNumero],
		" orbits of ",n,"-cells in the ",p,"-torsion subcomplex",
		", the stabilizers of which are of cardinality ",m,
		", namely the ones numero "
	);
      for j in [1..Length(pMultipleTorsionCells[n+1][PmultipleNumero])-1] do
	    Print(pMultipleTorsionCells[n+1][PmultipleNumero][j][2],", ");
      od;
      Print(pMultipleTorsionCells[n+1][PmultipleNumero][
	Length(pMultipleTorsionCells[n+1][PmultipleNumero])][2],".\n");
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
	  if Size(commonBoundary) = 1 then
	    Print("The boundary of ",n,"-cell numero ",firstCell[2]," is ",
	      celldata[n+1][firstCell[2]]!.BoundaryImage!.ListIFace,
	      " and that of ",n,"-cell numero ",secondCell[2]," is ",
	      celldata[n+1][secondCell[2]]!.BoundaryImage!.ListIFace,
	      " so they have ",n-1,"-cell numero ",commonBoundary[1],
		" in common.\n");
	    numberOfFusionCandidates := numberOfFusionCandidates +1;
	    fusionCandidates[numberOfFusionCandidates] :=
	      [firstCell, secondCell, commonBoundary[1]];	
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
		Print("The p-torsion subcomplex is bifurcated at the ",
			n-1,"-cell numero ",boundaryCell,".\n");
	  	RemoveSet(remainingFusionCandidates, cellTriple);
	  fi;
	fi;	
      od; 
  od; 
  return remainingFusionCandidates;
end;


IsInnerConjugate := function(G, A, B)
#################################################################
## Check if the subgroups A and B of G are conjugate inside G. ##
## For groups with high cardinality, we use GAP's function     ##
## IsConjugate(G, A, B), whilst for low cardinalities we       ##
## implement an algorithm that is faster for our matrix groups.##
#################################################################
local areConjugate, j, g, card, subcard, a, setOfElements, gAgInv, timeDifference;
 timeDifference := Runtime();
 if IsFinite(G) then
  card := Size(G); 
  if card < 1001 then

    j := 1;
    setOfElements := Set(G);
    areConjugate := false;
    subcard := Size(A);

    while j < card/subcard +1 and areConjugate = false do
	g := setOfElements[1];
	gAgInv := ConjugateGroup(A, g);

	if IsSubset(gAgInv, B) then
	  if IsSubset(B, gAgInv) then
		areConjugate := true;
	  fi;
	fi;
	for a in A do 
		RemoveSet(setOfElements, a*g);
	od;
	j := j+1;
    od;
  else areConjugate := IsConjugate(G, A, B);
  fi;  
 else areConjugate := IsConjugate(G, A, B);
 fi;
 Print("Computing conjugation took ",Runtime()-timeDifference," ms.\n");
 return areConjugate;
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
local firstCell, secondCell, boundaryCell, firstStabilizer, secondStabilizer, boundaryStabilizer, cellTriple, n, g, validated, remainingFusionCandidates;
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

    if IsSubset( boundaryStabilizer, firstStabilizer) then
      if IsSubset( boundaryStabilizer, secondStabilizer) then

	if IsSubset(firstStabilizer, secondStabilizer) then
		validated := true;
	else
		validated := IsConjugate( boundaryStabilizer,
				 firstStabilizer, secondStabilizer);
	fi;
      else Print("****Error: ",n,"-cell numero ",secondCell[2],
		 " not pointwise fixed.****\n");
      fi;
    else Print("****Error: ",n,"-cell numero ",firstCell[2],
		 " not pointwise fixed.****\n");
    fi;
    if validated then 
      if IsSubset(firstStabilizer, boundaryStabilizer) then 
	;  ## in this case,  we can clearly merge the cell triple       ##
	   ## without changing p-primary equivariant Farrell cohomology ##
		Print("The ",n,"-cells numero ",firstCell[2]," and ",
		 secondCell[2]," are subject to a cell fusion.\n"); 
      else
        if IsPnormal( boundaryStabilizer, p) then
	 g := Normalizer(boundaryStabilizer, 
		Center(SylowSubgroup( boundaryStabilizer, p)
	 ));
	 if IsSubset(firstStabilizer, g) and IsSubset(g, firstStabilizer) then
		validated := true;
	   ## if firstStabilizer is isomorphic to g, it is proven       ##
	   ## that we can merge the cell triple                         ##
	   ## without changing p-primary equivariant Farrell cohomology ##
	  else
	  	validated := IsInnerConjugate(
				boundaryStabilizer, firstStabilizer, g);
	  fi; 
	  if validated then
		Print("The ",n,"-cells numero ",firstCell[2]," and ",
		 secondCell[2]," are subject to a non-trivial fusion.\n"); 
          fi;
        else
	 validated := false;
        fi;
      fi;
    fi;
    if validated then ; ## the cell triple is subject ot a cell fusion. ##

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
local firstCell, secondCell, boundaryCell, cellTriple, n, j, reducedTorsionCells, vcd, mergedBoundary;
  mergedBoundary := [];
  reducedTorsionCells := StructuralCopy( torsionCells);
  vcd := Length(reducedTorsionCells) -1;
  for n in [0..vcd] do
	reducedTorsionCells[n+1] := Set(reducedTorsionCells[n+1]);
	mergedBoundary[n+1] := [];
  od;
  for cellTriple in fusionCandidates do 

    firstCell := cellTriple[1];
    secondCell := cellTriple[2];
    n := firstCell[1];
    boundaryCell := cellTriple[3];
    Print("The ",n-1,"-cell numero ",boundaryCell,
	" is amalgamated into the merging ",n,"-cells numero ",
	firstCell[2]," and ",secondCell[2],".\n");
    RemoveSet(reducedTorsionCells[n], [n-1,boundaryCell]);
    mergedBoundary[n+1][firstCell[2]] := Set(Concatenation(
	celldata[n+1][firstCell[2]]!.BoundaryImage!.ListIFace,
	celldata[n+1][secondCell[2]]!.BoundaryImage!.ListIFace
    ));
    RemoveSet(mergedBoundary[n+1][firstCell[2]], boundaryCell);
    Print("The merged boundary of these two cells is ",
		mergedBoundary[n+1][firstCell[2]],".\n");
   # RemoveSet(mergedBoundary[n+1][firstCell[2]], boundaryCell);
    #Print(mergedBoundary[n+1][firstCell[2]]);	 
  od;
  return reducedTorsionCells;
end;


getTorsionSubcomplex := function(groupName, p)
#####################################################################
## Here, p is the prime for which to take the torsion subcomplex.  ##
## We extract the cells the stabilizer of which contains p-torsion.##
#####################################################################
local vcd, stabilizerCardinalities, celldata,
torsionCells, numberOfTorsionCells, n, j;
   Read(Concatenation( DirectoriesPackageLibrary("HAP")[1]![1], 
			"Perturbations/Gcomplexes/",groupName));
   celldata := StructuralCopy(HAP_GCOMPLEX_LIST);
   vcd := Length(celldata) -1;
   Print("Extracting the ",p,"-torsion subcomplex of the ",
		vcd,"-dimensional ",groupName,"-cell complex ... \n");
   stabilizerCardinalities := admissibilityCheck(celldata);
   torsionCells := [];
   numberOfTorsionCells := [];
   for n in [0..vcd] do
	torsionCells[n+1] := [];
	numberOfTorsionCells[n+1] := 0;
	for j in [1..Length(celldata[n+1])] do
	   ## Check if the stabilizer contains p-torsion ##
	   if stabilizerCardinalities[n+1][j] mod p = 0 then
		Print("Extracted ",n,"-cell numero ",j,
			" of stabilizer cardinality ",
			stabilizerCardinalities[n+1][j],".\n");
		numberOfTorsionCells[n+1] 
			:= numberOfTorsionCells[n+1]+1;
	        torsionCells[n+1][numberOfTorsionCells[n+1]]
			:=[n, j];
	   fi;
	od;
   od;
   return
     [torsionCells, numberOfTorsionCells, celldata, stabilizerCardinalities];
end;


InstallGlobalFunction("TorsionSubcomplex", function(groupName, p)
############################################
local torsionCells, numberOfTorsionCells, celldata, sortedData;

	sortedData := getTorsionSubcomplex(groupName, p);
	torsionCells := sortedData[1];
	numberOfTorsionCells := sortedData[2];
	celldata := sortedData[3];
   return
     computeIncidenceMatrix(torsionCells, numberOfTorsionCells, celldata);
end);

InstallGlobalFunction("VisualizeTorsionSkeleton", function(groupName, p)
##################################################
local incidenceMatrix, graphData;
incidenceMatrix := TorsionSubcomplex(groupName, p);
graphData := StructuralCopy(IncidenceMatrixToGraph(incidenceMatrix));
Print("displayed in a separate window on screen now.\n");
GraphDisplay(graphData);
end);

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


getTerminalVertices := function(reducedTorsionCells, celldata)
##########################################################################
## We single out the "terminal" vertices of the p-torsion subcomplex :  ##
## Those with exactly one adjacent edge.                                ##
##########################################################################
local examinedVertex, examinedEdge, examinedEdgeData, numberOfAdjacencies, 
	adjacentEdge, reductionCandidates, numberOfTerminalVertices;
  reductionCandidates:= []; numberOfTerminalVertices := 0;

  for examinedVertex in reducedTorsionCells[1] do 

    examinedVertex := examinedVertex[2];
    numberOfAdjacencies := 0;

    for examinedEdge in reducedTorsionCells[1+1] do
	examinedEdgeData := celldata[1+1][examinedEdge[2]];

	if examinedVertex in examinedEdgeData!.BoundaryImage!.ListIFace then

	  numberOfAdjacencies := numberOfAdjacencies +1;
	  adjacentEdge := examinedEdge[2];
	fi;
    od;
    if numberOfAdjacencies = 1 then 

	Print("The vertex number ",examinedVertex," is terminal.\n");
	numberOfTerminalVertices := numberOfTerminalVertices +1;
	reductionCandidates[numberOfTerminalVertices] := 
				[examinedVertex, adjacentEdge];
    fi;
  od; 
  return reductionCandidates;
end;



checkGruenSwan := function( terminalVertices, celldata, p)
#########################################################################
## Keep only those of the terminal vertices, the stabilizers of which  ##
## satisfy the hypotheses of the Gruen/Swan theorem.   	 	       ##
#########################################################################
local boundaryCell, examinedEdge, EdgeStabilizer, boundaryStabilizer, cellPair, n, g, validated, verticesToReduce;
  verticesToReduce := StructuralCopy( terminalVertices);

  for cellPair in terminalVertices do 

    boundaryCell := cellPair[1];
    examinedEdge := cellPair[2];
    EdgeStabilizer := celldata[1+1][examinedEdge]!.TheMatrixStab;
    boundaryStabilizer := celldata[0+1][boundaryCell]!.TheMatrixStab;
    g := getIdentifier( [1,examinedEdge], boundaryCell, celldata);
    EdgeStabilizer := ConjugateGroup(EdgeStabilizer, g);


    if IsSubset(EdgeStabilizer, boundaryStabilizer) then 
	validated := true;  
	   ## in this case,  we can clearly cut off the edge            ##
	   ## without changing p-primary equivariant Farrell cohomology ##
    else
      if IsPnormal( boundaryStabilizer, p) then
	g := Normalizer(boundaryStabilizer, 
		Center(SylowSubgroup( boundaryStabilizer, p)
	));
	if IsSubset(EdgeStabilizer, g) and IsSubset(g, EdgeStabilizer) then
		validated := true;
	else
	  	validated := IsInnerConjugate(
				boundaryStabilizer, EdgeStabilizer, g);
	fi; 
	   ## if firstStabilizer is isomorphic to g, it is proven       ##
	   ## that we can cut off the edge                              ##
	   ## without changing p-primary equivariant Farrell cohomology ##
      else
	validated := false;
      fi;
    fi;
    if validated then 
		Print("The terminal vertex ",boundaryCell,
		" with edge ",examinedEdge," is subject to a reduction.\n"); 
     else RemoveSet(verticesToReduce, cellPair);
    fi;
  od; 
  return verticesToReduce;
end;


InstallGlobalFunction("ReduceTorsionSubcomplex", function(groupName, p)
###################################################
local torsionCells, numberOfTorsionCells, celldata, sortedData, stabilizerCardinalities, fusionCandidates, reducedTorsionCells,
terminalVertices;

   sortedData := getTorsionSubcomplex(groupName, p);
   torsionCells := sortedData[1];
   numberOfTorsionCells := sortedData[2];
   celldata := sortedData[3];
   stabilizerCardinalities := sortedData[4];

   sortedData := extractPmultipleTorsionCells(torsionCells, 
		numberOfTorsionCells, celldata,stabilizerCardinalities, p);
   fusionCandidates := pairsIntersection(sortedData, celldata);
   fusionCandidates := bifurcationFreeCells(fusionCandidates, torsionCells,
						 celldata);
   fusionCandidates := checkStabilizerConjugacy(fusionCandidates, celldata, 
						p);
   reducedTorsionCells := mergeCells( celldata, fusionCandidates, 
					torsionCells);
   terminalVertices := getTerminalVertices(reducedTorsionCells, celldata);
   terminalVertices := checkGruenSwan( terminalVertices, celldata, p); 
Print("\n At the following terminal vertices, the adjacent edge can be cut off without changing the equivariant Farrell cohomology : ",terminalVertices,"\n");  
   printIsotropyGroups(reducedTorsionCells, celldata);
   # return reducedTorsionCells;
end);
