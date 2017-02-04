InstallGlobalFunction(DavisComplex, function(CoxMat)
local DC, GEN, STAB, BOUNDARY, COEFF, CONJ, DIMENSION, numberOfCells,
irreducible, equalCoxeterMatrices, generateA, generateB, generateD,
H_4, F_4, E_6, E_7, E_8, irreducibleComponents, irr, spherical, extractSpherical, canonicalGenerators, generateStabs, gens, generateBoundary, generateCoeff, 
generateConj, getDimension, getNumberOfSimplices, chains, 
DavisComplexCombinatorial, aux, ExportAsContractibleGcomplex, IsSphericalCoxeterGroup, CreateCoxeterMatrix, MaximalSpherical ;




##########################################################
# DavisComplexes, a HAP sub-package in GAP
# by Ruben J. Sanchez-Garcia and Alexander D. Rahm,
#
# containing Ruben J. Sanchez-Garcia's implementation
# of the Davis complex for Coxeter groups.
#
##########################################################


##########################################################################
# Function IsSphericalCoxeterGroup
# The input is a Coxeter matrix, the output true or false
# Ruben J Sanchez-Garcia & Alexander Rahm
# GAP version 4.3 (TBU)
# May 2015
##########################################################################

#### SPHERICAL ####
## Check whether the Coxeter group represented by a
## Coxeter matrix M is spherical (finite)

## We need firstly the procedures irreducible, irreducibleComponents
## and equalCoxeterMatrices


#### OBTAIN IRREDUCIBLE COMPONENTS OF A COXETER GROUP ####
## if irreducible returns [ true, CoxeterMatrix ] 
## else returns [false, M1, M2] ] where Mi are two 
## submatrices corresponding to two components, 
## maybe not yet irreducible 
irreducible := function(CoxeterMatrix)
    local N, candidates, partition, tested, remove, s, c;
    
    N := Size( CoxeterMatrix );
    if N = 1 then return [true, CoxeterMatrix]; fi;
    
    candidates := [2..N];
    partition := [1];
    tested := [];
    
    while (candidates <> []) and (not IsEqualSet(partition, tested)) do
        s := Difference( partition, tested )[1]; 
        remove := [];
        for c in candidates do
            if CoxeterMatrix[s][c] <> 2 then
                AddSet( partition, c );
                AddSet( remove, c );
            fi;
        od;
        SubtractSet(candidates, remove);
        AddSet(tested, s);
    od;
    
    if candidates = [] then
        return [ true, CoxeterMatrix ];
    else
        return [ false, CoxeterMatrix{partition}{partition}, 
                            CoxeterMatrix{candidates}{candidates} ];
    fi;
    
end;


#### EQUAL AS COXETER MATRICES ####
## Determine whether there is a permutation p s.t.
## M1[p(i)][p(j)] = M2[i][j] all i,j
equalCoxeterMatrices := function( M1, M2)
    local perm, N, ContinueVariable, i, j;
    if DimensionsMat(M1) <> DimensionsMat(M2) then
        return false;
    elif 
        #checks whether entries(M1) = entries(M2)
        (IsEqualSet( Collected( Flat(M1) ), Collected( Flat(M2) ) ) = false) then
        return false;
    else
        N := Size(M1);
        for perm in SymmetricGroup(N) do
            ContinueVariable := true; 
            i := 1;
            while ContinueVariable and (i < N) do
                j := i+1;
                while ContinueVariable and (j <= N) do
                    if M1[i^perm][j^perm] = M2[i][j] then
                        j := j+1;
                    else
                        ContinueVariable := false; 
                    fi;
                od;
                i := i+1;
            od;
            if ContinueVariable = true then return true; fi;
        od;
        return false;
    fi;
end;


#### MATRICES SPHERICAL CASES ####
## Matrices for spherical cases: A_n, B_n, D_n,
generateA := function(N)
    local matrix, i;
    
    matrix := ( NullMat(N,N) + 2 ) - IdentityMat(N); #matrix all 2's except 1's at diagonal
    matrix[1][2] := 3;
    matrix[N][N-1] := 3;
    for i in [2..N-1] do
        matrix[i][i-1] := 3;
        matrix[i][i+1] := 3;
    od;
    return matrix;
end;

generateB := function(N)
    local matrix, i;
    
    matrix := ( NullMat(N,N) + 2 ) - IdentityMat(N); #matrix all 2's except 1's at diagonal
    matrix[1][2] := 4;
    matrix[N][N-1] := 4;
    for i in [2..N-1] do
        matrix[i][i-1] := 3;
        matrix[i][i+1] := 3;
    od;
    return matrix;
end;

generateD := function(N)
    local matrix, i;
    matrix := generateA(N);
    matrix[1][2] := 2;
    matrix[1][3] := 3;
    matrix[2][1] := 2;
    matrix[3][1] := 3;    
    return matrix;
end;

## Matrices for the cases H_4, F_4, E_6, E_7, E_8
H_4 := [ [1,5,2,2], [5,1,3,2], [2,3,1,3], [2,2,3,1] ];
F_4 := [ [1,3,2,2], [3,1,4,2], [2,4,1,3], [2,2,3,1] ];
E_6 := [ [1,3,2,2,2,2], [3,1,3,2,2,2], [2,3,1,3,3,2], [2,2,3,1,2,2], [2,2,3,2,1,3], [2,2,2,2,3,1] ];
E_7 := [ [1,3,2,2,2,2,2], [3,1,3,2,2,2,2], [2,3,1,3,3,2,2], [2,2,3,1,2,2,2], [2,2,3,2,1,3,2], [2,2,2,2,3,1,3], 
            [2,2,2,2,2,3,1] ];
E_8 := [ [1,3,2,2,2,2,2,2], [3,1,3,2,2,2,2,2], [2,3,1,3,3,2,2,2], [2,2,3,1,2,2,2,2], [2,2,3,2,1,3,2,2], 
            [2,2,2,2,3,1,3,2], [2,2,2,2,2,3,1,3], [2,2,2,2,2,2,3,1] ];

##########################################################################
# Auxiliary functions to CoxeterGroups.g
# Ruben J Sanchez
# University of Southampton 2001-04
# GAP version 4.3
# July 2004
##########################################################################



irreducibleComponents := function(CoxeterMatrix)
    local irr, comp1, comp2;
    
    irr := irreducible( CoxeterMatrix );
    if irr[1] = true then
        return [irr[2]];
    else
        comp1 := irreducibleComponents( irr[2] );    
        comp2 := irreducibleComponents( irr[3] );
        return Concatenation(comp1, comp2);
    fi;
end;


#### SPHERICAL FUNCTION ####
## Check whether the Coxeter group represented
## Coxeter matrix M is spherical (finite)
IsSphericalCoxeterGroup:=function(M)
    local N, i, j, sum, CoxeterMatrix, components, irre, 
            comp, fourOrfive, caseA, caseB, caseD;
    
    # finite for size 0 (trivial) or 1 (cyclic order 2)
    N := Size(M);
    if N <= 1 then
        return true;
    fi;
       
    # check whether some order is infinity
    for i in [1..N] do
        for j in [i+1..N] do
            if M[i][j] = 0 then return false;
            fi;
        od;
    od;           
        
    if N = 2 then 
        return true;
    fi;
        
    if N = 3 then
        sum := 1/M[1][2] + 1/M[1][3] + 1/M[2][3];
        if sum > 1 then 
            return true;
        else
            return false;
        fi;
    fi;
    
    # if Size >= 4, get irreducible components
    components := irreducibleComponents( M );
    irre := (Size(components) = 1); #irreducible?
    #if reducible, check recursively IsSphericalCoxeterGroup() on each component
    if not irre then
        for comp in components do
            if not IsSphericalCoxeterGroup(comp) then
                return false;
            fi;
        od;
        return true;
    #if irreducible, check exceptional cases
    else
        #easy checks: there is entries > 5? & #4's or 5's > 1 
        fourOrfive := false;
        for i in [1..N] do
            for j in [i+1..N] do
                if M[i][j] > 5 then
                    return false;
                elif (M[i][j] = 4) or (M[i][j] = 5) then
                    if fourOrfive = true then
                        return false;
                    else fourOrfive := true;
                    fi;
                fi;
            od;
        od;
        #check exceptional cases H_4, F_4, E_6, E_7, E_8
        if N <= 8 then
            if equalCoxeterMatrices(M, H_4) then
               return true;
            fi;
            if equalCoxeterMatrices(M, F_4) then
               return true;
            fi;
            if equalCoxeterMatrices(M, E_6) then
                return true;
            fi;
            if equalCoxeterMatrices(M, E_7) then
                return true;
            fi;
            if equalCoxeterMatrices(M, E_8) then
                return true;
            fi;
        fi;
        #check the cases A_n, B_n and D_n
        caseA := generateA(N); caseB := generateB(N); caseD := generateD(N);
        if equalCoxeterMatrices(M, caseA) then
                return true;
        fi;
        if equalCoxeterMatrices(M, caseB) then
                return true;
        fi;
        if equalCoxeterMatrices(M, caseD) then
                return true;
        fi;
        return false;
    fi;
       
end;






#### SPHERICAL FUNCTION ####
## Check whether the Coxeter group represented
## Coxeter matrix M is spherical (finite)
spherical := function(M)
    local N, i, j, sum, CoxeterMatrix, components, irre, 
            comp, fourOrfive, caseA, caseB, caseD;
    
    # finite for size 0 (trivial) or 1 (cyclic order 2)
    N := Size(M);
    if N <= 1 then
        return true;
    fi;
       
    # check whether some order is infinity
    for i in [1..N] do
        for j in [i+1..N] do
            if M[i][j] = 0 then return false;
            fi;
        od;
    od;           
        
    if N = 2 then 
        return true;
    fi;
        
    if N = 3 then
        sum := 1/M[1][2] + 1/M[1][3] + 1/M[2][3];
        if sum > 1 then 
            return true;
        else
            return false;
        fi;
    fi;
    
    # if Size >= 4, get irreducible components
    components := irreducibleComponents( M );
    irre := (Size(components) = 1); #irreducible?
    #if reducible, check recursively spherical() on each component
    if not irre then
        for comp in components do
            if not spherical(comp) then
                return false;
            fi;
        od;
        return true;
    #if irreducible, check exceptional cases
    else
        #easy checks: there is entries > 5? & #4's or 5's > 1 
        fourOrfive := false;
        for i in [1..N] do
            for j in [i+1..N] do
                if M[i][j] > 5 then
                    return false;
                elif (M[i][j] = 4) or (M[i][j] = 5) then
                    if fourOrfive = true then
                        return false;
                    else fourOrfive := true;
                    fi;
                fi;
            od;
        od;
        #check exceptional cases H_4, F_4, E_6, E_7, E_8
        if N <= 8 then
            if equalCoxeterMatrices(M, H_4) then
               return true;
            fi;
            if equalCoxeterMatrices(M, F_4) then
               return true;
            fi;
            if equalCoxeterMatrices(M, E_6) then
                return true;
            fi;
            if equalCoxeterMatrices(M, E_7) then
                return true;
            fi;
            if equalCoxeterMatrices(M, E_8) then
                return true;
            fi;
        fi;
        #check the cases A_n, B_n and D_n
        caseA := generateA(N); caseB := generateB(N); caseD := generateD(N);
        if equalCoxeterMatrices(M, caseA) then
                return true;
        fi;
        if equalCoxeterMatrices(M, caseB) then
                return true;
        fi;
        if equalCoxeterMatrices(M, caseD) then
                return true;
        fi;
        return false;
    fi;
       
end;



#### EXTRACT SPHERICAL SIMPLICES ####
## Remove from 'allchains' all chains containing a non-spherical subset
extractSpherical := function( allchains, numberOfGenerators, CoxeterMatrix )
    local L, i, j, subset, element, simplices_N, localchains;
    
    localchains := ShallowCopy(allchains);
    L := Size(CoxeterMatrix);
    for i in [L,L-1..2] do
        for subset in Combinations([1..L], i) do
            if spherical( CoxeterMatrix{subset}{subset} ) = false then
                for j in [1..Size(localchains)] do
                    SubtractSet( localchains[j], Filtered( localchains[j], C -> subset in C ) );
                od;
            fi;
        od;
    od;
return localchains;
end;         

#### GENERATING MATRICES (CANONICAL REPRESENTATION) ####
## Each matrix correpond to one generator ##
## Generator s_i correspond to matrix A[i] ##
## Input: Coxeter matrix M ##
canonicalGenerators := function(M)
    local cos, A, i, j, GEN;
    # auxiliary function: cos(n) = cos(Pi/n) = 1/2*(e^(Pi*I/n)-e^(Pi*I/n)^(2*n-1))
    # using GAP function E(n) = e^(2*Pi*I/n)
    cos := function(n)
        if n>0 then
            return (1/2)*( E(2*n) + E(2*n)^(2*n-1) );
        elif n = 0 then   #infinity case
            return 1;
        fi;
    end;

    # generating matrices, stored in array A
    A := [];
    GEN := Size(M);
    for i in [1..GEN] do
        A[i] := IdentityMat(GEN);
        for j in [1..GEN] do
            A[i][i][j] := A[i][i][j] + 2*cos(M[i][j]);
        od;
    od;
    return A;
end;


#### GENERATE STABILIZERS ####
# Generate list of stabilizers out of
# an array of 'chains' (such that chains[i] = chains of length i)
# with respect to the Coxeter matrix 'M'
# Note.- There cannot be a empty list of i-simplices
generateStabs := function( chains, M )
    local L, i, j, groupStabOfChain, simplices, simplex, stabs, gens;
    
    # Returns group corresponding to stabilizer of simplex
    # represented by chain, with generators strored in matrix gens
    groupStabOfChain := function( chain, gens )
        local l, smallest, set;
        if Size(chain[1]) < 1 then return Group( One ( gens[1] ) );
        else
            smallest := chain[1]; l := Size(smallest);
            for set in chain do
                if Size(set) < l then
                    smallest := set;
                    l := Size(set);
                fi;
            od;
            return Group( gens{smallest} );
        fi; 
    end;
    
    L := MaximumList( List( chains, x -> Size(x[1])) ); # length of longest chain
    simplices := []; stabs := [];
    gens := canonicalGenerators(M);
    for i in [1..L] do 
        j := 0; stabs[i] := []; 
        for simplex in chains[i] do
            j := j+1; 
            Add( stabs[i], groupStabOfChain( simplex, gens ) ); #stabs is NOT a set, must keep order
     #       SetName( stabs[i][j], Concatenation("Stab", String(i-1), "_", String(j)));
        od;
    od; 
    return stabs;
end;


#### GENERATE BOUNDARY ####
# Generated list of boundaries of each chain in
# the list of 'chains', where chains[i] = list of i-chains
generateBoundary := function( chains )
    local L, i, boundaryOfChain, simplices, simplex, boundary;
    
   # Returns list of chains = simplices in the boundary of N-simplex
   # given by chain, given by position in N-simplices
    boundaryOfChain := function( chain, simplices )
        local L, i, boundaryComponent, solution, temp, pos;
    
        L := Size(chain); 
        if L = 1 then return []; fi;
        if L < 1 then return fail; fi;
        solution := [];
        for i in Set( List( chain, Size ) ) do
            temp := First( chain, x -> Size(x) = i ); 
            boundaryComponent := Difference( chain, [temp] ); 
            pos := PositionSet( simplices, boundaryComponent ); 
            Add( solution, pos );
        od;
        return solution;
    end; 
    
    L := MaximumList( List( chains, x -> Size(x[1])) ); # length of longest chain
    simplices := []; boundary := [];
    for i in [2..L] do
        simplices[i-1] := ShallowCopy( chains[i] );
        boundary[i-1] := []; 
        for simplex in simplices[i-1] do
            Add( boundary[i-1], boundaryOfChain( simplex, chains[i-1] ) );
        od;
    od; 
    return boundary;
end;


#### GENERATE COEFF AND CONJ (both 'standard') ####
# Generate list COEFF = alternating +-1's and
# list CONJ = constant 1's, out of list BOUNDARY (replacing
# elements by appropiate +-1)
generateCoeff := function( boundary )
    local coeff;
    coeff := List( boundary, x -> List(x, y -> List (y, z -> (-1)^(Position(y,z)+1) )));
    return coeff;
end;
generateConj := function( boundary )
    local conj;
    conj := List( boundary, x -> List(x, y -> List (y, z -> 1)));
    return conj;
end;

#### DIMENSION AND NUMBER OF SIMPLICES (out of collection of chains 'simplices) ####
getDimension := function(simplices)
    return MaximumList( List( simplices, x -> Size(x[1])) ) - 1; # length of longest chain
end;
getNumberOfSimplices := function(simplices)
    return List(simplices, Size);
end;


##########################################################
# DavisComplexes, a HAP sub-package in GAP
# by Ruben J. Sanchez-Garcia and Alexander D. Rahm,
#
# containing Ruben J. Sanchez-Garcia's implementation
# of the Davis complex for Coxeter groups.
#
##########################################################

CreateCoxeterMatrix:= function( numberOfGenerators, CoxeterMatrixEntries)
# CreateCoxeterMatrix := function( numberOfGenerators, CoxeterMatrixEntries)
local GEN, M, m, k, i, j;

   if IsInt(numberOfGenerators) and numberOfGenerators > 0 then
        GEN := numberOfGenerators;
   else
	return fail;
	#break;
   fi;
   if Size(CoxeterMatrixEntries) = GEN*(GEN-1)/2 then

        M := IdentityMat(GEN); #to be completed into a Coxeter matrix
	k := 1;
        for i in [1..GEN-1] do
            for j in [i+1..GEN] do
                m := CoxeterMatrixEntries[k];
		if IsInt(m) and m > -1 then
                    	M[i][j] := Int(m);
		    	k := k+1;
                else
			return fail;
			#break;
		fi;
            od;
        od;    
   else
	return fail;
	#break;
   fi;
#### COMPLETE COXETER MATRIX (symmetric) ####
   for i in [1..GEN] do
        for j in [1..i-1] do
            M[i][j] := M[j][i];
        od;
   od;
   return M;
end;   





##########################################################################
# Program to compute the maximal spherical subgroups of a Coxeter group
# The input is a Coxeter matrix M
# The output is a list of subsets of the Coxeter generators labeled 1..N
# Ruben J Sanchez-Garcia & Alexander Rahm
# GAP version 4.3 (TBU)
# May 2015
##########################################################################
        
MaximalSpherical:= function(M)
    local S, subset, max_spherical, temp_spherical, l, m, found;

    # Coxeter generators
    S := [1..Size(M)];

    # if spherical return all generators
    if IsSphericalCoxeterGroup(M) then
        return [S];
    fi;

    # if not, evaluate each subset recursively
    max_spherical := [];
    for subset in Combinations(S,Length(S)-1) do
        temp_spherical :=  MaximalSpherical(M{subset}{subset});
        # for each returned spherical subset, check is not yet a subset in max_spherical
        for l in temp_spherical do
            found := false; 
            for m in max_spherical do
                if IsSubset(m, subset{l}) then #use subset{l} to keep consistent numbering
                    found := true;
                    #break;
                fi;
            od;
            if not(found) then
                Append(max_spherical, [subset{l}]);
            fi;
        od;
    od;

    return max_spherical;
end;

##########################################################################
# Program to compute the combinatorial structure of the Davis complex
# The input is a Coxeter matrix M
# The output is a list of subsets of the Coxeter generators labeled 1..N
# Ruben J Sanchez-Garcia & Alexander Rahm
# GAP version 4.3 (TBU)
# May 2015
##########################################################################
        

# Subfunction chains
    # computes all chains of subsets of S
    chains := function(S)
        local partialSolution, Hijos, hijo, addChainLink;

        # subfunction - add maximal set to the chains 
        addChainLink := function( chains, set )
            local C, newChain, solution;
            solution := Set(ShallowCopy(chains));
            for C in chains do
                newChain := Set(ShallowCopy(C));
                AddSet( newChain, set );
                AddSet( solution, newChain );
            od;
            AddSet( solution, [set] );
            return solution;
        end;

        # Trivial case
        if Size(S) = 0 then
            return [Combinations(S)];
        fi;

        # Reduction to smaller cases: recursive call
        partialSolution := [];
        Hijos := Combinations( S, Size(S)-1 );
        for hijo in Hijos do
            Append( partialSolution, chains(hijo) );
        od;

        partialSolution := addChainLink(partialSolution, S);

        return partialSolution;
    end;

DavisComplexCombinatorial := function(MaxSpherical)
    local s, i, aux, DC, L;

    # adds all chains of each maximal spherical, deleting duplicates 
    aux := [];
    for s in MaxSpherical do
        aux := Set(Concatenation(aux,chains(s)));
    od;

    # filters chains by length
    DC:= [];
    L := Maximum(List(aux,Length)); # length of longest chain
    for i in [1..L] do
        DC[i] := Filtered( aux, C -> Size(C) = i );
    od;

    return DC;
end;


ExportAsContractibleGcomplex := function(CoxeterMatrix, DC, GEN, STAB, BOUNDARY, COEFF, CONJ, DIMENSION, numberOfCells)
local	j, CoxeterGenerators,
        G, StabilizerGroups, Stabilizer,
        lnth,
        dims,Dimension,
        Boundary,
        boundaryList,
        Elts,
	Rot,Stab,
        RotSubGroups,Action, ActionRecord,
        TransMat,
        x, n,k,s,BI,SGN,tmp, LstEl;

# groupname:= Concatenation("Coxeter",String(GEN),"group");
# for i in [1..GEN-1] do
#      for j in [i+1..GEN] do
#           groupname :=  Concatenation( groupname,"-",String( CoxeterMatrix[i][j]));
#      od;
# od;  
# name:=groupname;

if HAP_GCOMPLEX_SETUP[1] then 
TransMat:=function(x); return x^-1; end;
else
TransMat:=function(x); return x; end;
fi;

lnth:=Size(DC)-1;

dims:=List([1..lnth+1],n->Size(DC[n]));

###################
Dimension:=function(n);
if n>lnth then return 0; fi;
return dims[n+1];
end;
###################

Elts:=[IdentityMat(GEN)];  
StabilizerGroups:=[];
RotSubGroups:=[];
boundaryList:=[];


#######
for n in [1..lnth+1] do
boundaryList[n]:=[];
StabilizerGroups[n]:=[];
RotSubGroups[n]:=[];
  for k in [1..Dimension(n-1)] do
  Append(Elts,Elements(STAB[n][k]));
  Add(StabilizerGroups[n],STAB[n][k]);
  Add(RotSubGroups[n],STAB[n][k]);
  od;
od;

CoxeterGenerators := [];
n := lnth; ## Use the stabilisers of the 
## (top-1 = n)-dimensional cells 
## as generators for the Coxeter group.
  for k in [1..Dimension(n-1)] do

  Append(CoxeterGenerators,Elements(STAB[n][k]));
  od;
####

Elts := SSortedList(Elts);
CoxeterGenerators := SSortedList(CoxeterGenerators);
CoxeterGenerators := Difference(CoxeterGenerators, [IdentityMat(GEN)]);

#######
for n in [1..lnth+1] do ## the cell dimension is n-1.
boundaryList[n]:=[];

  tmp := List([]);
  if n-1 > 0 then
    for j in [1..n] do
  	Append( tmp, [IdentityMat(GEN)]);
    od;		
  fi;
  LstEl:=tmp;

for k in [1..Dimension(n-1)] do
  if n-1 > 0 then   
       BI := BOUNDARY[n-1][k];
       SGN:=COEFF[n-1][k];
  else
       BI:= [];	
       SGN:= [];	
  fi;

   for s in [1..Length(BI)] do
          BI[s]:=[SGN[s]*BI[s], Position(Elts,LstEl[s])];
   od;
Add(boundaryList[n],BI);
od;
od;
####

ActionRecord:=[];
for n in [1..lnth+1] do
ActionRecord[n]:=[];
for k in [1..Dimension(n-1)] do
ActionRecord[n][k]:=[];
od;
od;

G:=Group(CoxeterGenerators);

####################
Boundary:=function(n,k);
if k>0 then
return boundaryList[n+1][k];
else
return NegateWord(boundaryList[n+1][-k]);
fi;
end;
####################

####################
Stabilizer:=function(n,k);
return StabilizerGroups[n+1][k];
end;
####################


####################
Action:=function(n,k,g)
local id,r,u,H,abk,ans;

abk:=AbsInt(k);

if not IsBound(ActionRecord[n+1][abk][g]) then 
H:=StabilizerGroups[n+1][abk];

if Order(H)=infinity then ActionRecord[n+1][abk][g]:=1;
#So we are assuming that any infinite stabilizer group acts trivially!!
else
######
id:=CanonicalRightCosetElement(H,Identity(H));
r:=CanonicalRightCosetElement(H,Elts[g]^-1);
r:=id^-1*r;
u:=r*Elts[g];

if u in RotSubGroups[n+1][abk] then  ans:= 1;
else ans:= -1; fi;

ActionRecord[n+1][abk][g]:=ans;
fi;
######
fi;

return ActionRecord[n+1][abk][g];
end;

##########################
return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:=Elts,
            group:=G,
            stabilizer:=Stabilizer,
            action:=Action,
            properties:=
            [["length",Maximum(1000,lnth)],
             ["characteristic",0],
             ["type","resolution"],
             ["reduced",true]]  ));

end; ## end of function ExportAsContractibleGcomplex
################################################
################################################



#InstallGlobalFunction(DavisComplex, function(CoxMat)
#local DC, GEN, STAB, BOUNDARY, COEFF, CONJ, DIMENSION, numberOfCells;
    if not IsMatrix(CoxMat) then CoxMat:=CoxeterMatrix(CoxMat); fi;
    GEN := Size(CoxMat);
    DC:= DavisComplexCombinatorial(MaximalSpherical(CoxMat));

#### EXTRACT SPHERICAL SIMPLICES ####
    DC := extractSpherical( DC, GEN, CoxMat );
#### REMOVE ANY EMPTY LIST OF CHAINS ####
    DC := Difference( DC, [[]] );
#### GENERATE STABILIZERS ####
    STAB := generateStabs( DC, CoxMat );

#### GENERATE BOUNDARY ####
    BOUNDARY := generateBoundary( DC );

#### GENERATE COEFF AND CONJ (both 'standard') ####
    COEFF := generateCoeff( BOUNDARY );
    CONJ := generateConj( BOUNDARY );

#### DIMENSION AND NUMBER OF SIMPLICES ####
    DIMENSION := getDimension(DC);
    numberOfCells := getNumberOfSimplices(DC);

return ExportAsContractibleGcomplex(CoxeterMatrix, DC, GEN, STAB, BOUNDARY, COEFF, CONJ, DIMENSION, numberOfCells);
# end;
end);


