##########################################################################
##
## Methods for 2x2 congruence subgroups of SL2


##########################################################################
##
## IndexInAmbientGroup( < G > )
##
## index of Gamma0(n) in SL(2,Z)
     InstallMethod(IndexInAmbientGroup,
     "index of Gamma0(n) in SL(2,Z)",
     [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
     function(G)
     local ind, n;
        n:=LevelOfCongruenceSubgroup(G);
        if n<>2 then TryNextMethod(); fi;
        if n=1 then
           ind:=1;
        else
           ind:=n*Product(List(SSortedList(Factors(n)), p->1+1/p));
        fi;
     return ind;
     end);

##########################################################################
##
## ProjectiveSpace( <G> )
InstallMethod(ProjectiveSpace,
     "Projective space",
     [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
     function(G)
     local n;
        if DimensionOfMatrixGroup(G)>2 then TryNextMethod(); fi;

        n := LevelOfCongruenceSubgroup(G);

        return FiniteProjectiveLine(n);
     end);

##########################################################################
##
## AmbientTransversal( <G> )
##
## Right transversal for a congruence subgroup G in its ambient group GG
     InstallMethod(AmbientTransversal,
     "Right transversal for a congruence subgroup G in its ambient group",
     [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
     function(G)
     local n, GG, poscan, cosetPos, transversal;
        if DimensionOfMatrixGroup(G)>2 then TryNextMethod(); fi;
        n:=LevelOfCongruenceSubgroup(G);
        if not IsPrime(n) then TryNextMethod(); fi;
        
     #Now handle case of prime level.
     GG:=AmbientGroupOfCongruenceSubgroup(G);

     cosetPos:=AmbientPosition(G);
     poscan:=function(g) return cosetPos(g^-1); end;
     
     transversal:=List([0..n-1],i->[[1,0],[-i,1]]);
     transversal:=Concatenation(transversal, [  [[0,1],[-1,0]] ]);

        return Objectify( NewType( FamilyObj( GG ),
                    IsHapRightTransversalSLnZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
                    rec( group := GG,
                         subgroup := G,
                         cosets:=transversal,
                         poscan:=poscan 
                    ));
     end);

##########################################################################
##
## AmbientPosition( <G> )
##
## Returns a function cosetPos(g) giving the position of the coset gG in 
## the ambient group. 
     InstallMethod(AmbientPosition,
     "Returns cosetPos(g) function for the congruence subgroup G",
     [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
     function(G)
     local cosetPos, n;
        if DimensionOfMatrixGroup(G)>2 then TryNextMethod(); fi;
        n:=LevelOfCongruenceSubgroup(G);
        if not IsPrime(n) then TryNextMethod(); fi;

        #Now handle case of prime level.
        ###########################################
        cosetPos:=function(g)
           if g[1][1] mod n =0 then return n+1; fi;
           return 1 +((g[2][1]*g[1][1]^-1) mod n);
        end;
        ###########################################
        return cosetPos;

     end);

##########################################################################
##
## AmbientRepresentation( <G> )
##
## Returns a function cosetPos(g) giving a canonical rpresentative of the 
## coset gG in the ambient group. 
     InstallMethod(AmbientRepresentation,
     "Returns cosetPos(g) function for the congruence subgroup G",
     [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
     function(G)
     local cosetRep, n;
        if DimensionOfMatrixGroup(G)>2 then TryNextMethod(); fi;
        n:=LevelOfCongruenceSubgroup(G);
        if not IsPrime(n) then TryNextMethod(); fi;

        #Now handle case of prime level.
        ###########################################
        cosetRep:=function(g)
        if g[1][1] mod n=0 then return [[0,-1],[1,0]]; fi;
        return [[1,0],[(g[2][1]*g[1][1]^-1) mod n,1]];
        end;
        ###########################################
        return cosetRep;

     end);


