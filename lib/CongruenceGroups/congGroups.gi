#############################################################################
##
## This file contains constructors of congruence subgroups and methods for 
## general congruence subgroups. Faster methods for specific groups of interest
## are installed elsewhere.

#############################################################################
##
## CongruenceSubgroupGamma0( < GG >, < n >)
InstallMethod( CongruenceSubgroupGamma0, "for integer matrix group and positive integer",
    [ IsIntegerMatrixGroup, IsPosInt],
    function(GG,n)
    local type, G, membership, ambientMembership, S,T,U;
    type := NewType( FamilyObj([One(GG)]),
                     IsGroup and
                     IsAttributeStoringRep and
                     IsFinitelyGeneratedGroup and
                     IsMatrixGroup and
                     IsHAPCongruenceSubgroup);
    G := rec();
    ObjectifyWithAttributes( G, type,
      DimensionOfMatrixGroup, DimensionOfMatrixGroup(GG),
      AmbientGroupOfCongruenceSubgroup, GG,
      OneImmutable, One(GG),
      IsIntegerMatrixGroup, IsIntegerMatrixGroup(GG),
      IsFinite, IsFinite(GG),
      StabilizerSubgroup, Group(One(GG)),
      LevelOfCongruenceSubgroup, n,
      IsHAPPrincipalCongruenceSubgroup, false,
      IsHAPCongruenceSubgroupGamma0, true,  
      IsHAPCongruenceSubgroupGamma1, false,
      IsHAPConjugatedCongruenceSubgroup, false 
      );
     if GG=SL(2,Integers) then
             S:=[[0,-1],[1,0]];; T:=[[1,1],[0,1]]; U:=S*T;
             GG:=AmbientGroupOfCongruenceSubgroup(G);
             GG!.GeneratorsOfMagmaWithInverses := [S*U,S*U^2];
     fi;


    if DimensionOfMatrixGroup(GG)=2 then
    ###################################################
    membership:=function(g)
    if not g in GG then return false; fi;
    if not g[2][1] mod n = 0  then return false;
    else return true; fi;
    end;
    ###################################################
    ###################################################
    ambientMembership:=function(g)
    if not g[2][1] mod n = 0  then return false;
    else  return true; fi;
    end;
    ###################################################
    fi;
    if DimensionOfMatrixGroup(GG)=3 then
    ###################################################
    membership:=function(g)
    if not g in GG then return false; fi;
    if not g[2][1] mod n = 0  then return false; fi;
    if not g[3][1] mod n = 0  then return false; fi;
    return true;
    end;
    ###################################################
    ###################################################
    ambientMembership:=function(g)
    if not g[2][1] mod n = 0  then return false; fi;
    if not g[3][1] mod n = 0  then return false; fi;
    return true; 
    end;
    ###################################################
    fi;
    if DimensionOfMatrixGroup(GG)=4 then
    ###################################################
    membership:=function(g)
    if not g in GG then return false; fi;
    if not g[2][1] mod n = 0  then return false; fi;
    if not g[3][1] mod n = 0  then return false; fi;
    if not g[4][1] mod n = 0  then return false; fi;
    return true;
    end;
    ###################################################
    ###################################################
    ambientMembership:=function(g)
    if not g[2][1] mod n = 0  then return false; fi;
    if not g[3][1] mod n = 0  then return false; fi;
    if not g[4][1] mod n = 0  then return false; fi;
    return true;
    end;
    ###################################################
    fi;

    if DimensionOfMatrixGroup(GG)>4 then
    Print("Gamma0 not yet implemented for matrices of dimension >4.\n");
    return fail;
    fi;

    G!.membership:=membership;
    G!.ambientMembership:=ambientMembership;

    ##################################################
    ##
    ## Remaining components to be computed in other functions.
    if DimensionOfMatrixGroup(GG)<4 then ProjectiveSpace(G);fi;
    AmbientTransversal(G);
    AmbientPosition(G); #A generic method will be used to construct the
    AmbientRepresentation(G); #these two functions except for cases with a
                         #bespoke implementation.
                         #cosetPos and poscan are used interchangably.

    return G;
        end);

#############################################################################
##
## PrincipalCongruenceSubgroup( < GG >, < n >)
InstallMethod( PrincipalCongruenceSubgroup, 
    "for integer matrix group and positive integer",
    [ IsIntegerMatrixGroup, IsPosInt],
    function(GG,n)
    local type, G, membership, ambientMembership, S,T,U;
    type := NewType( FamilyObj([One(GG)]),
                     IsGroup and
                     IsAttributeStoringRep and
                     IsFinitelyGeneratedGroup and
                     IsMatrixGroup and
                     IsHAPCongruenceSubgroup);
    G := rec();
    ObjectifyWithAttributes( G, type,
      DimensionOfMatrixGroup, DimensionOfMatrixGroup(GG),
      AmbientGroupOfCongruenceSubgroup, GG,
      OneImmutable, One(GG),
      IsIntegerMatrixGroup, IsIntegerMatrixGroup(GG),
      IsFinite, IsFinite(GG),
      LevelOfCongruenceSubgroup, n,
      IsHAPPrincipalCongruenceSubgroup, true,
      IsHAPCongruenceSubgroupGamma0, false,
      IsHAPCongruenceSubgroupGamma1, false,
      IsHAPConjugatedCongruenceSubgroup, false
      );
     if GG=SL(2,Integers) then
             S:=[[0,-1],[1,0]];; T:=[[1,1],[0,1]]; U:=S*T;
             GG:=AmbientGroupOfCongruenceSubgroup(G);
             GG!.GeneratorsOfMagmaWithInverses := [S*U,S*U^2];
     fi;

    ###################################################
    membership:=function(g)
    local i,j;
    if not g in GG then return false; fi;
    for i in [1..DimensionOfMatrixGroup(G)] do
    for j in [i+1..DimensionOfMatrixGroup(G)] do
        if not g[i][j] mod n = 0  then return false; fi;
        if not g[j][i] mod n = 0  then return false; fi;
    od;od;
    for i in [1..DimensionOfMatrixGroup(G)] do
        if not g[i][i] mod n = 1 then return false; fi;
    od;
    return true; 
    end;
    ###################################################
    ###################################################
    ambientMembership:=function(g)
    local i,j;
    for i in [1..DimensionOfMatrixGroup(G)] do
    for j in [i+1..DimensionOfMatrixGroup(G)] do
        if not g[i][j] mod n = 0  then return false; fi;
        if not g[j][i] mod n = 0  then return false; fi;
    od;od;
    for i in [1..DimensionOfMatrixGroup(G)] do
        if not g[i][i] mod n = 1 then return false; fi;
    od;
    return true; 
    end;
    ###################################################
    G!.membership:=membership;
    G!.ambientMembership:=ambientMembership;
    
    ###################################################
    ##
    ## Remaining components to be computer elsewhere.
    AmbientTransversal(G);
    AmbientPosition(G); #A generic method will be used to construct the
    AmbientRepresentation(G); #these two functions except for cases with a
                         #bespoke implementation.
                         #cosetPos and poscan are used interchangably.

    return G;
        end);

##########################################################################
##
## \IN
##
InstallOtherMethod(IN,
"membership for congruence subgroup",
[ IsMatrix, IsHAPCongruenceSubgroup ],
function(x, G)
return G!.membership(x);
end);

##########################################################################
##
## IsSubset(SL(n,Z),G)
##
## for an integer congruence subgroup G
InstallMethod( IsSubset,
    "for a natural SL_n(Z) and a congruence subgroup",
    [ IsNaturalSL, IsHAPCongruenceSubgroup and IsIntegerMatrixGroup ],
    0,
    function( SL, G )
    local x;
    if IsSubset(SL, AmbientGroupOfCongruenceSubgroup(G)) then
       return true;
    fi;
    for x in GeneratorsOfGroup(G) do
       if Determinant(x)<>1 then return false; fi;
    od;
    return true;
    end);

##########################################################################
##
## StabilizerSubgroup( <G> )
##
## Stabilizer subgroup for Gamma0 is chosen to be the identity group
     InstallMethod(StabilizerSubgroup,
     "Stabilizer subgroup of Gamma0(n) in SL(n,Z)",
     [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
     function(G)
        return Group(One(G));
     end);

##########################################################################
##
## StabilizerSubgroup( <G> )
##
## Stabilizer subgroup for PrincipalGamma is chosen to be the identity group
## except in the 2x2 case.
     InstallMethod(StabilizerSubgroup,
     "Stabilizer subgroup of Gamma(n) in SL(n,Z)",
     [ IsIntegerMatrixGroup and IsHAPPrincipalCongruenceSubgroup ],
     function(G)
     local S, T, U, GG;
        if DimensionOfMatrixGroup(G)=2 then
           if Size(LevelOfCongruenceSubgroup(G))<=2 then
              return Group(One(G));
           else
              return Group( [[0,-1],[1,0]] * [[1,1],[0,1]] );
           fi;
        fi;
        if DimensionOfMatrixGroup(G)>2 then
        return Group(One(G));
        fi;

     end);


##########################################################################
##
## AmbientTree( <G> )
##
## Returns a tree representing the cosets of G in the ambient group GG. This
## method assumes very little of G. For specific cases of interest a faster 
## method will be installed elsewhere. For the principal congruence group in
## SL(2,Z) we work with a non-trivial Stabilizer in order to produce a free
## generating set (and the tree is thus not quite a coset tree). This is 
## currently the only case where we use non-trivial stabilizer, though other 
## SL(2,Z) cases should really be "upgraded" to non-trivial stabilizer. Due
## to potentially non-trivial stabilizers we have to distinguish between 
## ambient trees and ambient transversals. 
     InstallMethod(AmbientTree,
     "Constructs a coset tree for a generic congruence subgroup (slow method)",
     [ IsHAPCongruenceSubgroup ],
     function(G)
     local GG,tree,InGmodU,Ugrp,v,p,g,s,n,q,vv,gens,
      nodes, leaves, ambientGenerators, InLowDegreeNodesModG,
      one, genTriples, vertex2word, triple2word, csts,U;

     GG:=AmbientGroupOfCongruenceSubgroup(G);
     ambientGenerators:=GeneratorsOfGroup(GG);
     one:=One(GG);

     Ugrp:=StabilizerSubgroup(G);
     Ugrp:=Elements(Ugrp);
     tree:=[1 ];
     genTriples:=[];
     cnt:=1;
     leaves:=NewDictionary(one,true,GG);
     nodes:=[one];
     AddDictionary(leaves,one,1);

     if Size(Ugrp)>1 then
        ###########################################
        InGmodU:=function(g)
        local u;
        for u in Ugrp do
           if G!.ambientMembership(u*g) then return true; fi;;
        od;
        return false;
        end;
        ###########################################
     else
        InGmodU:=G!.ambientMembership;
     fi;

     ###########################################
     InLowDegreeNodesModG:=function(g)
     local x,gg,B1,B2;
     gg:=g^-1;

     for x in nodes do
     if InGmodU(x*gg) then return x; fi;
     od;

     return false;
     end;
     ###########################################


     ############Tree Construction########################
     while Size(leaves)>0 do
     vv:=leaves!.entries[1];
     v:=vv[1];
         for s in [1..Length(ambientGenerators)] do
             g:=v*ambientGenerators[s];
             q:=InLowDegreeNodesModG(g);
             if q=false then
               Add(nodes,g);
               AddDictionary(leaves,g,Length(nodes));
               p:=LookupDictionary(leaves,v);
               Add(tree,[p, s]);
             else Add(genTriples,[v,g,q]);
             fi;
         od;
     RemoveDictionary(leaves,v);
     od;
     #####################################################

     G!.nodes:=nodes;

     #####################################################
     vertex2word:=function(v)
     local word;
     word:=one;
     while not v=1 do
        word:=ambientGenerators[tree[v][2]]*word;
        v:=tree[v][1];
     od;
     return word;
     end;
     #####################################################

     #####################################################
     triple2word:=function(x)
     local u,uu,g,q,c,ans;

     ans:=[];
     for u in Reversed(Ugrp) do
        c:=x[2]*u*vertex2word(  Position(nodes,x[3])  )^-1;
        if G!.ambientMembership(c) then Add(ans,c); fi;
     od;

     return ans;
     end;
     #####################################################

     Unbind(tree[1]);
     #MUST ADJUST THE FOLLOWING FOR SIZE(STAB)>1
     genTriples:=List(genTriples,x->triple2word(x));
     genTriples:=Concatenation(genTriples);
     genTriples:=List(genTriples,x->Minimum(x,x^-1));
     genTriples:=SSortedList(genTriples);
     SetGeneratorsOfMagmaWithInverses(G,genTriples);
     return tree;
     end);

############################################################################
##
## GeneratorsOfGroup( <G> )
##
## returns the generators of a congruence subgroup
     InstallMethod(GeneratorsOfGroup,
     "Generators of a congruence subgroup",
     [ IsHAPCongruenceSubgroup ],
     function(G)
     AmbientTree(G);
     return GeneratorsOfMagmaWithInverses(G);
     end);

############################################################################
##
## AmbientTransversal( <G> )
##
## returns a transversal for a generic congruence subgroup G in the ambient 
## group GG. The method will invoke AmbientTree if StabilizerSubgroup(G)=1.
     InstallMethod(AmbientTransversal,
     "returns a transversal for G in the ambient group (slow method)",
     [ IsHAPCongruenceSubgroup ],
     function(G)
     local GG,tree,v,p,g,s,n,q,vv,gens,
      nodes, nodesinv, leaves, ambientGenerators, InLowDegreeNodesModG,
      one, poscan;

     GG:=AmbientGroupOfCongruenceSubgroup(G);

     if Size(StabilizerSubgroup(G))=1 then
         AmbientTree(G);
         nodes:=G!.nodes;
         poscan:=AmbientPosition(G); 

     else
     ################################################################
     ################################################################
     ambientGenerators:=GeneratorsOfGroup(GG);
     one:=One(GG);
     tree:=[1 ];
     cnt:=1;
     leaves:=NewDictionary(one,true,GG);
     nodes:=[one];
     AddDictionary(leaves,one,1);

     ###########################################
     InLowDegreeNodesModG:=function(g)
     local x,gg,B1,B2;
     gg:=g^-1;
     for x in nodes do
        if G!.ambientMembership(x*gg) then return x; fi;
     od;
     return false;
     end;
     ###########################################

     ############Tree Construction########################
     while Size(leaves)>0 do
     vv:=leaves!.entries[1];
     v:=vv[1];
     for s in [1..Length(ambientGenerators)] do
        g:=v*ambientGenerators[s];
        q:=InLowDegreeNodesModG(g);
        if q=false then
         Add(nodes,g);
         AddDictionary(leaves,g,Length(nodes));
            p:=LookupDictionary(leaves,v);
            Add(tree,[p, s]);
        fi;
     od;
     RemoveDictionary(leaves,v);
     od;
     #####################################################
     
     nodesinv:=List(nodes,g->g^-1);

     ####################################################
     poscan:=function(x)
     local i;
     for i in [1..Length(nodesinv)] do
        if G!.ambientMembership(  x*nodesinv[i]  ) then return i; fi;
     od;
     end;
     ####################################################
     fi;
     ################################################################
     ################################################################
     

     return Objectify( NewType( FamilyObj( GG ),
                    IsHAPRightTransversalCongruenceSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := GG,
               subgroup := G,
               cosets:=nodes,
               poscan:=poscan ));
     end);

############################################################################
##
## RightTransversal( <GG> , <G> )
##
## returns a right transversal of a generic congruence subgroup in its ambient
## group
     InstallMethod(RightTransversal,
     "right transversal of a congruence subgroup in its ambient group",
     [ IsMatrixGroup , IsHAPCongruenceSubgroup ],
     function(GG,G)
     local T, poscan, Tinv;
     if GG=AmbientGroupOfCongruenceSubgroup(G) then
        return AmbientTransversal(G);
     fi;
 
     if not IsHAPCongruenceSubgroup(GG) then
         TryNextMethod();
     fi;
     
     T:=AmbientTransversal(G);
     T:=Filtered(T!.cosets, x->GG!.ambientMembership(x));
     Tinv:=List(T, x->x^-1);

     ########################
     poscan:=function(g)
     local i;
     #gg:=g^-1;     
     for i in [1..Length(T)] do
         #if T[i]*gg in G then return i; fi;
         if g*Tinv[i] in G then return i; fi;
     od;
     end;
     ########################

     return 
     Objectify( NewType( FamilyObj( AmbientGroupOfCongruenceSubgroup(GG) ),
                IsHAPRightTransversalCongruenceSubgroup and IsList and
                IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := GG,
               subgroup := G,
               cosets:=T,
               poscan:=poscan ));
     
     end);

##########################################################################
##
## AmbientTree( <G> )
##
## Returns a tree representing the cosets of G in the ambient group GG.
     InstallMethod(AmbientTree,
     "Constructs a coset tree for a congruence subgroup in the ambient group, and uses it to obtain a generating set (faster method)",
     [ IsHAPCongruenceSubgroup and HasAmbientPosition and HasAmbientRepresentation ],
     function(G)
     local GG, ambientGenerators, tree,InGmodU,Ugrp,U,v,p,g,s,n,q,vv,
      leaves,genTriples,generators, InLowDegreeNodesModG, csts, cnt,
      vertex2word, one, cosetPos, cosetRep, triple2word,
      i,j,u,c,a,b;

        if DimensionOfMatrixGroup(G)>2 then TryNextMethod(); fi;
        n:=LevelOfCongruenceSubgroup(G);

     GG:=AmbientGroupOfCongruenceSubgroup(G);
     ambientGenerators:=GeneratorsOfGroup(GG);
     one:=One(GG);
     cosetPos:=AmbientPosition(G);
     cosetRep:=AmbientRepresentation(G);

     Ugrp:=StabilizerSubgroup(G);
     Ugrp:=Elements(Ugrp);

     tree:=[ ];
     genTriples:=[];
     cnt:=1;
     leaves:=NewDictionary(one,true,GG);
     csts:=[];
     csts[cosetPos(one)]:=1;
     AddDictionary(leaves,one,cosetPos(one));

     ###########################################
     InLowDegreeNodesModG:=function(g)
     local pos;
     pos:=cosetPos(g);
     if not IsBound(csts[pos]) then return false; else;
     return pos; fi;
     end;
     ###########################################


     #########The work horse##############################
     while Size(leaves)>0 do
         vv:=1*leaves!.entries[1];
         v:=vv[1];
         p:=cosetPos(v);
         for s in [1..Length(ambientGenerators)] do
            g:=ambientGenerators[s]*v;
            q:=InLowDegreeNodesModG(g);
            if q=false then
                 q:=cosetPos(g);
                 AddDictionary(leaves,g,q);
                 csts[q]:=1;
                 tree[q]:=[p,s];
            else Add(genTriples,[p,s,v,g]);  #Usp=Uq mod G
            fi;
         od;
         RemoveDictionary(leaves,v);
     od;
     #####################################################

     #####################################################
     vertex2word:=function(v)
     local word;
     word:=one;
     while IsBound(tree[v]) do
         word:=word*ambientGenerators[tree[v][2]];
         v:=tree[v][1];
     od;
     return word;
     end;
     #####################################################

     #####################################################
     triple2word:=function(x)
     local u,uu,g,q,c;
     for u in Ugrp do
         c:=x[4]^-1*u*vertex2word(cosetPos(x[4]));
         if G!.ambientMembership(c) then return c; fi;
     od;
     return fail;  #This should never happen
     end;
     #####################################################

     #MUST ADJUST FOR SIZE(STAB)>1
     genTriples:=List(genTriples,x->triple2word(x));
     genTriples:=List(genTriples,x->Minimum(x,x^-1));
     genTriples:=SSortedList(genTriples);
     SetGeneratorsOfMagmaWithInverses(G,genTriples);

     return tree;
     end);

##########################################################################
##
## AmbientPosition( <G> )
##
## Returns a function cosetPos(g) giving the position of the coset gG in 
## the ambient group. 
     InstallMethod(AmbientPosition,
     "Returns cosetPos(g) function for the congruence subgroup G",
     [ IsIntegerMatrixGroup and IsHAPPrincipalCongruenceSubgroup ],
     function(G)
     local cosetPos, n, sln, Uelts, RR, R, g, one;
     n:=LevelOfCongruenceSubgroup(G);
     sln:=SL(DimensionOfMatrixGroup(G) , Integers mod n );
     one:=One(sln);
     Uelts:=StabilizerSubgroup(G);
     Uelts:=Elements(Uelts)*one;
     RR:=Enumerator(sln);;
     R:=[];;
     for g in RR do
         Add(R,Minimum(List(Uelts,u->u*g)));
     od;
     R:=SSortedList(R);

     ###################################################
     cosetPos:=function(g)
     local gg;
     gg:=g*one;
     gg:=Minimum(List(Uelts,u->u*gg));
     return Position(R,gg);
     end;
     ###################################################

     return cosetPos;

     end);

##########################################################################
##
## AmbientRepresentation( <G> )
##
## Returns a function cosetRep(g) giving a representative of the coset gG in 
## the ambient group. 
     InstallMethod(AmbientRepresentation,
     "Returns cosetPos(g) function for the congruence subgroup G",
     [ IsIntegerMatrixGroup and IsHAPPrincipalCongruenceSubgroup ],
     function(G)
     local cosetRep, n, sln, Uelts,  one;
     n:=LevelOfCongruenceSubgroup(G);
     sln:=SL(DimensionOfMatrixGroup(G) , Integers mod n );
     one:=One(sln);
     Uelts:=StabilizerSubgroup(G);
     Uelts:=Elements(Uelts)*one;

     ###################################################
     cosetRep:=function(g)
     local gg;
     gg:=g*one;
     gg:=Minimum(List(Uelts,u->u*gg));
     return gg;
     end;
     ###################################################

     return cosetRep;

     end);

##########################################################################
##
## IndexInAmbientGroup( < G > )
##
## index of congruence subgroup in ambient group 
     InstallMethod(IndexInAmbientGroup,
     "index of Gamma0(n) in ambient",
     [ IsHAPCongruenceSubgroup ],
     function(G)
     return Length(AmbientTransversal(G));
     end);

##########################################################################
##
## Length( T )
##
## length of a right transversal of congruence subgroups
     InstallOtherMethod(Length, 
     "length of a transversal",
     [IsHAPRightTransversalCongruenceSubgroup],
     function(T)
     return Length(T!.cosets);
     end);

##########################################################################
##
     InstallOtherMethod(PositionCanonical,
     "for HapRightTransversals of subrougs ",
     [IsHAPRightTransversalCongruenceSubgroup,IsObject],
     function(R,x)
     return R!.poscan(x);
     end);

##########################################################################
##
     InstallOtherMethod(ELM_LIST,
     "for HAPRightTransversalCongruenceSubgroup",
     [IsHAPRightTransversalCongruenceSubgroup,IsInt],
     function(R,i)
     return R!.cosets[i];
     end);

##########################################################################
##
## AmbientPosition( <G> )
##
     InstallOtherMethod(AmbientPosition,
     "for a generic congruence subgroup (slow method via ambient tree)",
     [ IsHAPCongruenceSubgroup ],
     function(G)
     local cosetPos, nodesinv;

     nodesinv:=List(G!.nodes,g->g^-1);
     ###################################################
     cosetPos:=function(x) 
     local i;
     for i in [1..Length(nodesinv)] do
        if G!.ambientMembership(  x*nodesinv[i]  ) then return i; fi;
     od;
     end;
     ####################################################
     return cosetPos;
     end);

##########################################################################
##
## AmbientRepresentation( <G> )
##
     InstallOtherMethod(AmbientRepresentation,
     "for a generic congruence subgroup (slow method via ambient tree)",
     [ IsHAPCongruenceSubgroup ],
     function(G)
     local cosetRep,T;
     T:=AmbientTransversal(G);
     cosetRep:=function(g) return T[T!.poscan(g)]  ; end;
     return cosetRep;
     end);

############################################################################
##
## AmbientTreeDisplay( G )
##
## Display the ambient tree.
     InstallGlobalFunction(AmbientTreeDisplay,
     function(G)
     local A, T, i, j,L;
     T:=AmbientTree(G);
     A:=0*IdentityMat(Length(T));;
     for i in [1..Length(T)] do
         if IsBound(T[i]) then
             A[i][T[i][1]]:=1;
             A[T[i][1]][i]:=1;
         fi;
     od;
     L:=Filtered([1..Length(A)],i-> not IsZero(A[i]));
     A:=A{L};
     A:=TransposedMat(A)*1;
     A:=A{L};

     A:=IncidenceMatrixToGraph(A);
     GraphDisplay(A);
     end);


############################################################################
##
## IntersectionWithCongugate( <H> , g )
##
InstallMethod(IntersectionWithConjugate,
" returns the intersection of H and H^g ",
[ IsHAPCongruenceSubgroup, IsMatrix ],
function(H,g)
local  type, G, gg, membership, HH, S, U, T;

    type := NewType( FamilyObj([One(H)]),
                     IsGroup and
                     IsAttributeStoringRep and
                     IsFinitelyGeneratedGroup and
                     IsMatrixGroup and
                     IsHAPCongruenceSubgroup);
    G := rec();
    ObjectifyWithAttributes( G, type,
      DimensionOfMatrixGroup, DimensionOfMatrixGroup(H),
      AmbientGroupOfCongruenceSubgroup, AmbientGroupOfCongruenceSubgroup(H),
      OneImmutable, One(H),
      IsIntegerMatrixGroup, false,  #it could be integer but mostly it won't be
      IsFinite, IsFinite(H),
      StabilizerSubgroup, Group(One(H)),
      LevelOfCongruenceSubgroup, infinity,  #Think!
      IsHAPPrincipalCongruenceSubgroup, false,
      IsHAPCongruenceSubgroupGamma0, false,
      IsHAPCongruenceSubgroupGamma1, false,
      IsHAPConjugatedCongruenceSubgroup,true 
      );

gg:=g^-1;

###################################################
membership:=function(x)
if not x in H then return false; fi;
return x^gg in H;
end;
###################################################

HH:=Intersection(H,H^g);
G!.membership:=membership;
G!.ambientMembership:=membership;

SetGeneratorsOfMagmaWithInverses(G,GeneratorsOfGroup(HH));

G!.AmbientGroupOfCongruenceSubgroup:= AmbientGroupOfCongruenceSubgroup(H);
return G;

end);

############################################################################
##
## IntersectionW( <H> , <K> )
##
InstallOtherMethod(Intersection2,
" returns the intersection of two congruence subgroups H and K ",
[ IsHAPCongruenceSubgroup, IsHAPCongruenceSubgroup ],
function(H,K)
local  type, G, membership, ambientMembership;

    if DimensionOfMatrixGroup(H)<>DimensionOfMatrixGroup(K) then
        Print("The matrix groups be be of the same dimension.");
        return fail;
    fi; 

    type := NewType( FamilyObj([One(H)]),
                     IsGroup and
                     IsAttributeStoringRep and
                     IsFinitelyGeneratedGroup and
                     IsMatrixGroup and
                     IsHAPCongruenceSubgroup);
    G := rec();
    ObjectifyWithAttributes( G, type,
      DimensionOfMatrixGroup, DimensionOfMatrixGroup(H),
      AmbientGroupOfCongruenceSubgroup, AmbientGroupOfCongruenceSubgroup(H),
      OneImmutable, One(H),
      IsIntegerMatrixGroup, IsIntegerMatrixGroup(H) and IsIntegerMatrixGroup(K),  
      IsFinite, IsFinite(H) or IsFinite(K),
      StabilizerSubgroup, Group(One(H)),
      LevelOfCongruenceSubgroup, Lcm(LevelOfCongruenceSubgroup(K),LevelOfCongruenceSubgroup(H)),  
      IsHAPPrincipalCongruenceSubgroup, false,
      IsHAPCongruenceSubgroupGamma0, false,
      IsHAPCongruenceSubgroupGamma1, false,
      IsHAPConjugatedCongruenceSubgroup,false
      );

###################################################
membership:=function(x)
if not H!.membership(x) then return false; fi;
if not K!.membership(x) then return false; fi;
return true;
end;
###################################################

###################################################
ambientMembership:=function(x)
if not H!.ambientMembership(x) then return false; fi;
if not K!.ambientMembership(x) then return false; fi;
return true;
end;
###################################################

G!.membership:=membership;
G!.ambientMembership:=ambientMembership;
G!.AmbientGroupOfCongruenceSubgroup:= AmbientGroupOfCongruenceSubgroup(H);

AmbientTransversal(G);
AmbientPosition(G);
AmbientRepresentation(G);

return G;

end);

