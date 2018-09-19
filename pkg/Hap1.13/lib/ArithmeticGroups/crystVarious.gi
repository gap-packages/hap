
#######################################################################
#0
#F  IsIntList
##  Input: A list L 
##         
##  Output: True if L is a list of integers
##          False otherwise   
##

InstallGlobalFunction(IsIntList,
function(list)
local i;

    for i in list do
        if not IsInt(i) then 
            return false;
        fi;
    od;
    return true;
end);
#######################################################################

#######################################################################
#0
#F  VectorToCrystMatrix
##  Input: A n-dimensional vector v 
##          
##  Output: A pure translation crystallographic matrix 
##             
##
InstallGlobalFunction(VectorToCrystMatrix,      
function(v)
local M,n;

    v:=Flat(v);
    n:=Length(v);
    M:=IdentityMat(n+1);
    Add(v,1);
    Remove(M);
    Add(M,v);
    return M;
end);
################### end of VectorToCrystMatrix #######################

######################################################################
#0
#F  CrystTranslationMatrixToVector
##  Input: A pure translation crystallographic matrix g
##         
##  Output: An n-dimensional vector v  
##             
##
InstallGlobalFunction(CrystTranslationMatrixToVector,
function(g)
local n,v;

    n:=Length(g);
    v:=g[n];
    v:=Flat(v);
    Remove(v);
    return v;
end);
################### end of CrystTranslationMatrixToVector ############

######################################################################
#0
#F  TranslationSubGroup
##  Input: A crystallographic group G 
##         
##  Output: Translation subgroup of G 
##             
##
InstallGlobalFunction(TranslationSubGroup,
function(G)
local B,SbGrp,trsltbasis;

    B:=TranslationBasis(G);
    if not IsBound(G!.TranslationBasis) then return false;fi;
    B:=G!.TranslationBasis;
    trsltbasis:=List(B,w->VectorToCrystMatrix(Flat(w)));
    SbGrp:=Group(trsltbasis);
    SbGrp!.TranslationBasis:=B;
    SetIsCrystTranslationSubGroup(SbGrp,true);
    return SbGrp;
end);

################### end of TranslationSubGroup #######################




######################################################################
#0
#F  Method in
##  Input: A matrix g and a translation subgroup G of 
##                                           a crystallographic group  
##         a positive integer p 
##  Output: True if g in G 
##          False if g not in G   
##
InstallOtherMethod(\in, 
        "for TranslationSubGroup of a CrystGroup",
          [IsMatrix,IsCrystTranslationSubGroup],
function(g,G)
local B,v,n;

    n:=DimensionSquareMat(g)-1;
    if not LinearPartOfAffineMatOnRight(g)=IdentityMat(n) then 
         return false;
    fi;
    B:=G!.TranslationBasis;
    v:=CrystTranslationMatrixToVector(g);

    return IsIntList(v*TransposedMat(B)^-1);

end);

################### end of Method g in G ###########################

####################################################################
#0
#F  IsCrystSameOrbit
##  Input: Two points u, v in R^n and a group G.
##          
##  Output: return True if u, v in the same orbit. 
##          Otherwise returns False. 
##             
##
InstallGlobalFunction(IsCrystSameOrbit,
function(arg)
    local G,T,H,u,v,B,x,w;
    
    G:=arg[1];
    if Length(arg)=3 then
        H:=TranslationSubGroup(G);
        T:=RightTransversal(G,H);
        B:=H!.TranslationBasis;
        u:=arg[2];
        v:=arg[3];
    else
        B:=arg[2];
        T:=arg[3];
        u:=arg[4];
        v:=arg[5];
    fi;
    u:=Flat(u);
    v:=Flat(v);
    Add(u,1);
    Add(v,1);
    for x in T do
        w:=u*x-v;
        w:=Flat(w);
        Remove(w);
        if IsIntList(w*TransposedMat(B)^-1) then 
            return x*VectorToCrystMatrix(w)^-1;
        fi;
    od;
    return false;
end);

################### end of IsCrystSameOrbit ##############################

##########################################################################
#0
#F  CombinationDisjointSets
##  Input: A list of k positive integers $(a_i)$. 
##
##  Output: A 2-dimensional array L such that 0 <= L[i][j] < a_j. 
##             
##

InstallGlobalFunction(CombinationDisjointSets,
function(arg)
local b,list,n1,i,g,h;

    g:=arg[1];
    if g=[] then return [[]];fi;
    n1:=g[1];
    h:=g{[2..Length(g)]};
    list:=[];
    b:=CombinationDisjointSets(h);
    for i in [0..(n1-1)] do
        Append(list,List(b,w->AddFirst(w,i)));
    od;
    return list;
end);

################### end of CombinationDisjointSets #######################

##########################################################################
#0
#F  AddFirst
##  Input: A list w and an element g.  
##          
##  Output: List with g in the first position 
##             
##
InstallGlobalFunction(AddFirst,
function(list,g)                 # add g in the first position in list
local w;
    w:=[g];
    Append(w,list);
    return w;
end);
##
################### end of AddFirst ######################################

##########################################################################
#0
#F  IsCrystSufficientLattice
##  Input: A lattice's basis B and a transversal S of the translation
##         subgroup in crystallographic group G. 
##          
##  Output: True if G acts on the lattice, otherwise return False. 
##             
##
InstallGlobalFunction(IsCrystSufficientLattice,
function(B,SS,T)
local 
    v,x,w,B1,c,i,A,Origin,S1,S2,S;
    S:=StructuralCopy(List(SS));
    Append(S,GeneratorsOfGroup(T));
    Origin:=0*B[1];
    Origin:=Flat(Origin);
    Add(Origin,1);
    A:=StructuralCopy(B);
    v:=Sum(A)/2;
    v:=Flat(v);
    Add(v,1);
    for x in S do
        w:=Flat(v*x-v);
        if not IsIntList(w*A^-1) then 
            return false;
        fi;
    od;
    for i in [1..Length(A)] do
        A[i]:=Flat(A[i]);
        Add(A[i],1);
    od;
    for x in S do
        B1:=List(A,w->w*x);
        c:=Flat(Origin*x-Origin);
        B1:=List(B1,w->w-c);
        for x in B1 do
            Remove(x);
        od;
        S2:=Set(B1);
        for x in B do
            if not ((x in S2) or (-x in S2)) then 
            return false;
            fi;
        od;
    od;
    return true;
end);
################### end of IsCrystSufficientLattice ######################

##########################################################################
#0
#F  CrystFinitePartOfMatrix
##  Input: A crystallographic matrix g 
##          
##  Output: Finite part of g. 
##             
##
InstallGlobalFunction(CrystFinitePartOfMatrix,
function(g)
local 
    x,w,i;
    w:=[];
    for i in [1..(Length(g)-1)] do
        x:=Flat(g[i]);
        Remove(x);
        Add(w,x);
    od;
    return w;
end);
################### end of CrystFinitePartOfMatrix #######################

##########################################################################
#0
#F  ResolutionBoundaryOfWordOnRight
##  Input: A free resolution R, degree n, a list of word w 
##          
##  Output: The boundary of w respects to the right action. 
##             
##
InstallGlobalFunction(ResolutionBoundaryOfWordOnRight,
function(R,n,W)
local 
    x, DW, Boundary, Dimension,Elts,pos, ans,H;

    Dimension:=R!.dimension;
    Boundary:=R!.boundary;
    Elts:=R!.elts;
    DW:=[];

    for x in W do
        ans:=Boundary(n,x[1]);
        ans:=List(ans, a->[a[1],Elts[a[2]]]);
        ans:=List(ans, a->[a[1],a[2]*Elts[x[2]]]);
        Append(DW,ans);
    od;

    DW:= AlgebraicReduction(DW);
    for x in DW do
        if not x[2] in Elts then 
            Add(Elts,x[2]);
        fi;
    od;
    DW:=List(DW,x->[x[1],Position(Elts,x[2])]);
    DW:=List(DW,x->[R!.Sign(n-1,x[1],x[2])*x[1],x[2]]);
    for x in DW do
        H:=R!.stabilizer(n-1,AbsInt(x[1]));
        pos:=Position(R!.elts,CanonicalRightCountableCosetElement(
                                               H,R!.elts[x[2]]));
        if pos=fail then 
            Add(R!.elts,CanonicalRightCountableCosetElement(
                                               H,R!.elts[x[2]])); 
            x[2]:=Length(R!.elts);
        else 
            x[2]:=pos;
        fi;
    od;
    DW:=List(DW,x->[R!.Sign(n-1,x[1],x[2])*x[1],x[2]]);
    DW:= AlgebraicReduction(DW);
    return DW;
end);
################### end of ResolutionBoundaryOfWordOnRight ###############

##########################################################################
#0
#F  CrystCubicalTiling
##  Input: Dimension n 
##          
##  Output: A list of some cubical tiling in n-dimensional space 
##             
InstallGlobalFunction(CrystCubicalTiling,
function(n)
local 
    combin,x,w,Til,i;
    combin:=Combinations([1..n],2);
    Til:=[];
    Add(Til,IdentityMat(n));
    for i in [1..Length(combin)] do
        w:=combin[i];
        x:=IdentityMat(n);
        x[w[1]][w[1]]:=-1/2;
        x[w[1]][w[2]]:=-1/2;
        x[w[2]][w[1]]:=1;
        x[w[2]][w[2]]:=-1;
        Add(Til,x);
    od;
    return Til;
end);
################### end of CrystCubicalTiling ############################

##########################################################################
#0
#F  AverageInnerProduct
##  Input: An affine group G and 2 vector u,v   
##          
##  Output: the avarage inner product of u and v. 
##             
InstallGlobalFunction(AverageInnerProduct,
function(G,u,v)
local 
    i,sum,n,Elts;
    n:=Order(G);
    Elts:=Elements(G);
    sum:=0;
    for i in [1..n] do
        sum:=sum+(u*Elts[i])*(v*Elts[i]);
    od;
    sum:=sum/n;
    return sum;
end);

################### end of AverageInnerProduct ##########################

#########################################################################
#0
#F  OrthogonalizeBasisByAverageInnerProduct
##  Input: A basis B and group G 
##          
##  Output: An orthogonal basis B' and a a matrix change of basis   
##             
##
InstallGlobalFunction(OrthogonalizeBasisByAverageInnerProduct,
function(B,G)
local 
    Project,i,j,A,n;
    A:=StructuralCopy(B);
    n:=Length(B);
    if not RankMat(B)=Length(B) then 
        Print("Input is not a basis");
        return fail;
    fi;

    Project:=function(u,v)    
                        #This operator projects the vector v orthogonally 
                          #onto the line spanned by vector u
        return (AverageInnerProduct(G,u,v)/AverageInnerProduct(G,u,u))*u;
    end;

    for i in [2..n] do
        for j in [1..(i-1)] do
            B[i]:=B[i]-Project(B[j],A[i]);
        od;
    od;
    for i in [1..n] do
        B[i]:=(1/Sqrt(AverageInnerProduct(G,B[i],B[i])))*B[i];
    od;
    return B;
end);

########## end of OrthogonalizeBasisByAverageInnerProduct ################

##########################################################################
#0
#F  CrystMatrix
##  Input: nxn matrix M 
##          
##  Output: the crystallographic form of M. 
##             
##
InstallGlobalFunction(CrystMatrix,
function(M)
local 
    i,n,x,N;
    N:=StructuralCopy(M);
    if IsMatrix(N) then
        n:=Length(N);
        x:=0*N[1];
        Add(x,1);
        for i in [1..n] do
            Add(N[i],0);
        od;
        Add(N,x);
    else
        N:=VectorToCrystMatrix(N);
    fi;
return N;
end);

################### end of CrystMatrix ###################################

##########################################################################
#0
#F  IsRigid
##  Input: A HAP G-complex
##  
##  Output: Either ``true'' or ``false''.
##
##
InstallGlobalFunction(IsRigid,
function(C)
local 
    i,j,bdr, w,s;

    for i in [1..5000] do #SLOPPY!! 	        
	for j in [1..C!.dimension(i)] do
            bdr:=C!.boundary(i,j);
            for w in bdr do
            s:=ConjugateGroup(C!.stabilizer(i-1,AbsInt(w[1])),
                                      C!.elts[w[2]]^-1); 
            if not IsSubgroup(s,C!.stabilizer(i,j)) then Print([i,j],"\n");
 return false; fi;
            od;
        od;
        i:=i+1;
    od; 
    return true;
end);
#######################End of IsRigid  ###################################


##########################################################################
#0
#F  IsRigidOnRight
##  Input: A HAP G-complex
##  
##  Output: Either ``true'' or ``false''.
##
##
InstallGlobalFunction(IsRigidOnRight,
function(C)
local 
    i,j,L,bdr,intst;

    i:=1;    
    while C!.dimension(i)>0 do
        for j in [1..C!.dimension(i)] do
            bdr:=C!.boundary(i,j);
            L:=List(bdr,w->Elements(ConjugateGroup(C!.stabilizer(i-1,AbsInt(w[1])),C!.elts[w[2]])));
            intst:=Intersection(L);
            if not Elements(C!.stabilizer(i,j))=Elements(intst) then
Print([i,j]);
                return false;
            fi;
        od;
        i:=i+1;
    od; 
    return true;
end);
#######################End of IsRigidOnRight  ###################################
