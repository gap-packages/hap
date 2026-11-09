
###########################################################
###########################################################
InstallGlobalFunction(Gamma0_SL3ZTopRationalHomology,
function(n)
    local FirstCondition, SecondCondition, ThirdCondition_1, ThirdCondition_2, A, B, C, i, j, M, SolutionSpace,FiniteProjectivePlane,ProjPlane;

#Returns the dimension of the degree 3 rational homology of the 
#congruence subgroup Gamma_0(n) in SL(3,Z). 

##########################################
##########################################

    ProjPlane := FiniteProjectiveSpace(2,n);

    FirstCondition := function(p)
        local v;
        v := [-p[2],p[1],p[3]] mod n;
        return ProjPlane.RepOf(v);
    end;

    SecondCondition := function(p)
        local v;
        v := [p[3],p[1],p[2]] mod n;
        return ProjPlane.RepOf(v);
    end;

    ThirdCondition_1 := function(p)
        local v;
        v := [-p[2],p[1]-p[2],p[3]] mod n;
        return ProjPlane.RepOf(v);
    end;

    ThirdCondition_2 := function(p)
        local v;
        v := [p[2]-p[1],p[1],p[3]] mod n;
        return ProjPlane.RepOf(v);
    end;

    A := SparseIdentityMat(Length(ProjPlane.Reps));
    B := SparseIdentityMat(Length(ProjPlane.Reps));
    C := SparseIdentityMat(Length(ProjPlane.Reps));

    for i in [1..Length(ProjPlane.Reps)] do
        j := Position(ProjPlane.Reps,FirstCondition(ProjPlane.Reps[i]));
        SparseMatAddToEntry(A,i,j,1);
        j := Position(ProjPlane.Reps,SecondCondition(ProjPlane.Reps[i]));
        SparseMatAddToEntry(B,i,j,-1);
        j := Position(ProjPlane.Reps,ThirdCondition_1(ProjPlane.Reps[i]));
        SparseMatAddToEntry(C,i,j,1);
        j := Position(ProjPlane.Reps,ThirdCondition_2(ProjPlane.Reps[i]));
        SparseMatAddToEntry(C,i,j,1);
    od;
    M:=SparseMatConcatenation(A,B);Unbind(A);Unbind(B);
    M:=SparseMatConcatenation(M,C);Unbind(C);
    return Length(ProjPlane.Reps) - RankMatDestructive(M);
end);
###########################################################
###########################################################

