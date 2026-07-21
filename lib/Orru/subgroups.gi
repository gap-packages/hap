
###################################################################
###################################################################
InstallGlobalFunction(HAP_GenericCongruenceSubgroup,
function(family,integer,ring,ideal)
    local R, n, one, mat, type, G, I, fam;

# This function returns a template for a congruence subgroup G in some 
# classical group such as GL(n,R), SL(n,R), Sp_n(R) over a ring of integers R.
# The subgroup G is defined using some condition involving an ideal I<=R. 
# The input "family" is a string such as "SL", "GL", Sp". In this code we
# denote any such family by "SL".
# 
# Various components such as "membership, "tree", "cosetPos" will ineed to be 
# implemented using a specific secondary function that depends on the specific 
# family and ring.
#
# The *fast* method will use a formula to first implement "cosetPos" and then
# compute the coset "tree", The *slow* method will use a general brute force 
# algorithm to first comput the "tree" and then implement "cosetPos".

    R := ring;
    n := integer;

    if IsInt(ideal) then
        I := IdealByGenerators(R,[ideal]); 
    else
        I := ideal; 
    fi;
    if not IsBound(I) then
        return fail; 
    fi;

    one := One(R);
    mat := one*IdentityMat(n);
    fam := Concatenation(family,"(",String(n),", ",Name(R)," )");

    type := NewType( FamilyObj([ mat ]),
                     IsGroup and
                     IsAttributeStoringRep and
                     IsFinitelyGeneratedGroup and
                     IsMatrixGroup and
                     IsHapCongruenceSubgroup);

    G:=rec(
        ringOfIntegers:=R,	#Matrices in the subgroup G lie in the ring R.
        level:=I,              	#This ideal I<=R is used to define G. 
        fam:=fam,		#A string like "SL(n,R)".
        dimension:=n,		#A group of nxn matrices
        ambientGroup:=fail,     #A supergroup like SL(n,Integers)
        ambientGenerators:=fail, #Generators for the ambient group
        membership:= fail,      #true if a matrix g lies in G.
        membershipLight:=fail,  #true if a matrix g in SL(n,R) lies in G.
        gens:=fail,		#"Nice" generating set for SL(n,R).
        tree:= fail,   		#Coset tree of G with respect to gens.
        generators:= fail,	#Generating set for G. 
        index:=fail,   		#Index of G in SL(n,R).
        cosetRep:= fail,        #CosetRep(g) represents g*G for g in SL(n,R).
        cosetPos:= fail,       	#cosetPos(g) is the position of the coset g*G. 
        ugrp:= Group(mat),	#The trivial (or vertex stabilizer) group.
        name:="Congruence subgroup");

        ObjectifyWithAttributes(G, type,
				DimensionOfMatrixGroup, n
							);

    return G;
end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallMethod( HAPCongruenceSubgroupGamma0_alt,
"for integer n and integer m",
[IsInt, IsInt],
function(n,m)
    local G, sl, membership, membershipLight, CosetRep, CanonicalRep, CosetPos, ProjLine,
    CosetOfInt, S, T, U, UU;
    if not (n=2 and m>0) then TryNextMethod(); fi;

    #The following implements G=Gamm0(m) < SL(2,Z)

    sl := SL(2,Integers);
    G := HAP_GenericCongruenceSubgroup("SL",2,Integers,m);
    UU := Units(Integers mod m);  #Graham
    UU := List(UU, u->Int(u));  #Graham
    G!.Units:=UU; #Graham

    ###################################################
    membership:=function(g)
        if not g in sl then
            return false;
        fi;
        if not g[2][1] mod m = 0  then
            return false;
        else
            return true;
        fi;
    end;
    ###################################################
    ###################################################
    membershipLight:=function(g)
        if not g[2][1] mod m = 0 then
            return false;
        else
            return true;
        fi;
    end;
    ###################################################
    ####################
    S:=[[0,-1],[1,0]];;
    T:=[[1,1],[0,1]];
    U:=S*T;
    ####################


    G!.membership := membership;
    G!.membershipLight := membershipLight;
    G!.level := m;

    G!.ugrp := Group([[1,0],[0,1]]);
    G!.name := "CongruenceSubgroupGamma0";
    if m = 1 then
        G!.index := m;
    else
        G!.index := m*Product(List(SSortedList(Factors(m)), p->1+1/p));
    fi;

    ProjLine := FiniteProjectiveLine_alt(m);
    CosetPos := function(g)
        local v, vv, U, u, w;
        v := [g[1][1], g[2][1]];
        vv := List(v, x -> x mod m);
        #U := Units(Integers mod m);
        for u in G!.Units do
            #w := List(vv, x -> (Int(u)*x) mod m);
            w := List(vv, x -> (u*x) mod m);
            if w in ProjLine.Reps then                  #Takes all the time
                return Position(ProjLine.Reps,w);       #
            fi;
        od;
    end;

#    CosetRep := function(i)   #Should input a group element

    CosetOfInt := function(i)
        local a, c, b, d, gg;
        a := ProjLine.Reps[i][1];
        c := ProjLine.Reps[i][2];
        if a = 0 then
            return [[0,-1],[1,0]];
        fi;
        gg := Gcdex(a,c);
        b := -gg.coeff2;
        d :=  gg.coeff1;
        return [[a,b],[c,d]];
    end;

    CosetRep:=function(g);
        return CosetOfInt(CosetPos(g));
    end;

    G!.cosetRep := CosetRep;   
    G!.cosetPos := CosetPos;   
    G!.ambientGenerators:=[S,S*U];
    G!.transversal:=List([1..Length(ProjLine.Reps)],i->CosetOfInt(i)^-1);

    G := ObjectifyWithAttributes(G, TypeObj(G),
        IsIntegerMatrixGroup, true,
        IsFinite, false);

    return G;
end);
###################################################################
#NEW
InstallMethod( HAPCongruenceSubgroupGamma0,
"for integer n and integer m",
[IsInt, IsInt],
function(n,m)
    local G, sl, membership, membershipLight, CosetRep, CanonicalRep, CosetPos, ProjLine, CosetOfInt, S, T, U, UU;
    if not (n=2 and m>0) then TryNextMethod(); fi;

    #The following implements G=Gamm0(m) < SL(2,Z)

    sl := SL(2,Integers);
    G := HAP_GenericCongruenceSubgroup("SL",2,Integers,m);

    UU := Units(Integers mod m);  #Graham
    UU := List(UU, u->Int(u));  #Graham
    G!.Units:=UU; #Graham


    ###################################################
    membership:=function(g)
        if not g in sl then
            return false;
        fi;
        if not g[2][1] mod m = 0  then
            return false;
        else
            return true;
        fi;
    end;
    ###################################################
    ###################################################
    membershipLight:=function(g)
        if not g[2][1] mod m = 0 then
            return false;
        else
            return true;
        fi;
    end;
    ###################################################
    ####################
    S:=[[0,-1],[1,0]];;
    T:=[[1,1],[0,1]];
    U:=S*T;
    ####################


    G!.membership := membership;
    G!.membershipLight := membershipLight;
    G!.level := m;

    G!.ugrp := Group([[1,0],[0,1]]);
    G!.name := "CongruenceSubgroupGamma0";
    if m = 1 then
        G!.index := m;
    else
        G!.index := m*Product(List(SSortedList(Factors(m)), p->1+1/p));
    fi;

    ProjLine := FiniteProjectiveLine(m);

    CanonicalRep := function(g)
        local v, vv, U, d, dd, x, y;
        v := [g[1][1], g[2][1]];
        vv := List(v, x -> x mod m);
        #U := Units(Integers mod m);
        if vv[1] mod m = 0 then
            return [0,1];
        #elif ZmodnZObj(vv[1],m) in U then
        elif vv[1] in G!.Units then
            return [1,(Inverse(vv[1]) mod m)*vv[2] mod m];
        else
            d := Gcd(vv[1],m);
            dd := m/d;
            x := vv[1]/d;
            y := vv[2]/x mod dd;
            while not Gcd(d,y) = 1 do
                y := y + dd;
            od;
            return [d, y];
        fi;
    end;

    CosetPos := function(g)
        local w;
        w := CanonicalRep(g);
        return Position(ProjLine,w);
    end;

#    CosetRep := function(i)   #Should input a group element

    CosetOfInt := function(i)
    local a, c, b, d, gg;
        a := ProjLine[i][1];
        c := ProjLine[i][2];
        if a = 0 then
            return [[0,-1],[1,0]];
        fi;
        gg := Gcdex(a,c);
        b := -gg.coeff2;
        d :=  gg.coeff1;
        return [[a,b],[c,d]];
    end;

    CosetRep:=function(g);
        return CosetOfInt(CosetPos(g));
    end;

    G!.cosetRep := CosetRep;   
    G!.cosetPos := CosetPos;   
    G!.ambientGenerators:=[S,S*U];
    G!.transversal:=List([1..Length(ProjLine)],i->CosetOfInt(i)^-1);

    G := ObjectifyWithAttributes(G, TypeObj(G),
        IsIntegerMatrixGroup, true,
        IsFinite, false);

    return G;
end);
###################################################################

###################################################################
###################################################################
#THERE SHOULD BE JUST ONE IMPLEMENTATION FOR ALL CASES. 
#NEED TO ADJUST THIS.
InstallMethod(HAPCongruenceSubgroupTree,
"Coset tree for congruence subgroup",
[IsHapCongruenceSubgroup],
function(G)
    if not (G!.dimension = 2 and G!.ringOfIntegers = Integers) then
        TryNextMethod();
    fi;
    HAP_SL2ZSubgroupTree_fast(G);
end);
###################################################################
###################################################################

###################################################################
###################################################################

InstallOtherMethod( RightTransversal,
"Right transversal of congruence subgroup G in SL(2,Z)",
[IsMatrixGroup, IsHapCongruenceSubgroup],
100000,
function(G,H)
    if not (H!.dimension = 2 and Name(G) = "SL(2,Integers)") then
        TryNextMethod();
    fi;

    #HAPCongruenceSubgroupTree(H);

    return HAP_TransversalCongruenceSubgroupInAmbientGroup(G,H);
end);
###################################################################
###################################################################

#SL3

InstallMethod( HAPCongruenceSubgroupGamma0,
"for SL(3,Z)",
[IsInt, IsInt],
function(n,m)
    local G,sl,membership,membershipLight, ProjPlane, CosetRep, CosetPos, MatrixInSL3_Hermite, S, T, U, UU, CosetOfInt;
    
    if not (n = 3 and m > 0) then
        TryNextMethod();
    fi;

    sl := SL(3, Integers);
    G  := HAP_GenericCongruenceSubgroup("SL", 3, Integers, m);

    UU := Units(Integers mod m);  #Graham
    UU := List(UU, u->Int(u));  #Graham
    G!.Units:=UU; #Graham


    membership := function(g)
        if not g in sl then
            return false;
        fi;
        if g[2][1] mod m <> 0 then
            return false;
        fi;
        if g[3][1] mod m <> 0 then
            return false;
        fi;
        return true;
    end;

    membershipLight := function(g)
        if g[2][1] mod m <> 0 then
            return false;
        fi;
        if g[3][1] mod m <> 0 then
            return false;
        fi;
        return true;
    end;

    G!.membership := membership;
    G!.membershipLight := membershipLight;
    G!.level := m;
    G!.name := "CongruenceSubgroupGamma0";

    MatrixInSL3_Hermite := function(v)
        local Herm;
        Herm := HermiteNormalFormIntegerMatTransform([[v[1]],[v[2]],[v[3]]]);
        return Inverse(Herm!.rowtrans);
    end;

    ProjPlane := FiniteProjectivePlane(m);
    CosetPos := function(g)
        local v, vv, u, w;
        v := [g[1][1], g[2][1], g[3][1]];
        vv := List(v, x -> x mod m);
        #U := Units(Integers mod m);
        for u in G!.Units do
            #w := List(vv, x -> (Int(u)*x) mod m);
            w := List(vv, x -> (u*x) mod m);
            if w in ProjPlane.Reps then
                return Position(ProjPlane.Reps,w);
            fi;
        od;
    end;
#    CosetRep := function(i)    #Should input a group element
    CosetOfInt:=function(i)
    local x,y,z;
        x := ProjPlane.Reps[i][1];
        y := ProjPlane.Reps[i][2];
        z := ProjPlane.Reps[i][3];

        return MatrixInSL3_Hermite([x,y,z]);
    end;

    CosetRep:=function(g);
        return CosetOfInt(CosetPos(g));
    end;

    S := [ [1,0,1], [0,-1,-1], [0,1,0] ];
    T := [ [0,1,0], [0,0,1], [1,0,0]];
    U := [[0,1,0], [1,0,0], [-1,-1,-1]];


    G!.cosetRep := CosetRep;  
    G!.cosetPos := CosetPos;  
    G!.ambientGenerators:=[S,T,U];
    G!.transversal:=List([1..Length(ProjPlane.Reps)],i->CosetOfInt(i)^-1);

    G := ObjectifyWithAttributes(G, TypeObj(G),
        IsIntegerMatrixGroup, true,
        IsFinite, false);

    return G;
end);

InstallMethod(HAPCongruenceSubgroupTree,
"Coset tree for congruence subgroup",
[IsHapCongruenceSubgroup],
function(G)
if not (G!.dimension=3 and G!.ringOfIntegers=Integers) then TryNextMethod(); fi;
HAP_SL3ZSubgroupTree_fast(G);
end);

InstallOtherMethod(\in,
"membership test for HapCongruenceSubgroups",
[IsMatrix, IsHapCongruenceSubgroup and IsGroup],
1000000,  #There must be a better way to ensure this method is used!
function(g,G)
return  G!.membership(g);
end);

InstallMethod(GeneratorsOfGroup,
"Generating set for HapCongruenceSubgroups",
[IsHapCongruenceSubgroup and IsGroup],
1000000,  #There must be a better way to ensure this method is used!
function(G)
HAPCongruenceSubgroupTree(G);
return  G!.GeneratorsOfMagmaWithInverses;
end);

InstallOtherMethod( RightTransversal,
"Right transversal of congruence subgroup G in SL(3,Z)",
[IsMatrixGroup, IsHapCongruenceSubgroup],
1000000,
function(G,H)
    if not (H!.dimension = 3 and Name(G) = "SL(3,Integers)") then
        TryNextMethod();
    fi;

    #HAPCongruenceSubgroupTree(H);

    return HAP_TransversalCongruenceSubgroupInAmbientGroup(G,H);
end);

