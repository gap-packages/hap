#DeclareGlobalFunction("SpunKnotCompliment");
#DeclareGlobalFunction("SpunAboutInitialHyperPlane");
###############################################################################
######################## Spinning Knot Complements ############################
###############################################################################
################# Input: an integer pair (n,p) corresponding ##################
######################## to the pth prime knot on n crossings, ################
######################## denoted by K. ########################################
################ Output: a 5-dimensional regular CW-complex S(K*) #############
######################## homotopy equivalent to the complement ################
######################## of the spinning of K* (a knotted arc formed ##########
######################## from K by removing an unknotted segment) #############
######################## about a plane. #######################################
###############################################################################
InstallGlobalFunction(SpunKnotComplement,
function(k)

    local
          d, K, DeletedPair, pair, C, U, C2, omicron;

    d:=[0,0,1,1,2,3,7,21,49,165,552];

    if IsList(k)
        then
        if not (IsInt(k[1]) and IsInt(k[2]) and SignInt(k[1])=1 and SignInt(k[2])=1)
            then
            Error("the input must be a pair of positive integers.\n");
        elif k[1]>11
            then
            Error("only knots with less than 12 crossings are stored in HAP.\n");
        elif k[2]>d[k[1]]
            then
            Error("no such prime knot exists.\n");
        fi;
        K:=PureCubicalKnot(k[1],k[2]);
    elif IsPureComplex(k)
        then
        K:=ShallowCopy(k);
    else
        Error("input must be an integer pair or a knot.\n");
    fi;

    DeletedPair:=function(K)
        # inputs a cubical knot
        # outputs a list of cubical complexes [C,U]
        # C being the complement of K with a certain line removed
        # (except for its end-points)
        # U is a plane intersecting C at just these endpoints
        local
              C, array, 0row, i, pos, a, b, subarray, j;

        C:=PureComplexComplement(K);
        array:=ShallowCopy(C!.binaryArray);
        0row:=0; # location of the line to be altered

        while 0row=0
            do
            for i in [1..Length(array[2])]
                do
                if 0 in array[2][i]
                    then
                    0row:=0row+i;
                    break;
                fi;
            od;
        od;

        pos:=Positions(array[2][0row],0);
        a:=pos[1]; # the endpoints of 0row
        b:=pos[Length(pos)];

        array[2][0row]:=0*ShallowCopy(array[2][0row])+1;
        array[2][0row][a]:=0;
        array[2][0row][b]:=0;
        array[1][0row]:=ShallowCopy(array[2][0row]);

        subarray:=0*ShallowCopy(array); # the plane about which the knot
        subarray[1]:=ShallowCopy(array[1]); # complement will be spun

        return [PureCubicalComplex(array), PureCubicalComplex(subarray)];

    end;

    pair:=DeletedPair(K);
    C:=ShallowCopy(pair[1]);
    U:=ShallowCopy(pair[2]);
    C2:=ContractedComplex(C,U);

    omicron:=RegularCWMap(C2,U);

    return Spin(omicron);

end);
###############################################################################
###############################################################################

###############################################################################
######################## Spinning About Hyperplane ##############################
###############################################################################
################# Input: Pure cubical Complex K ###############################
################ Output: a regular CW-complex S(K) ############################
######################## homotopy equivalent to the space got  ################
######################## by spinning K about the 'initial' hyperplane #########
###############################################################################
InstallGlobalFunction(SpunAboutHyperplane,
function(K)
local U,i,inc;

U:=1*K!.binaryArray;
for i in [2..Length(U)] do    ##CHANGED
U[i]:=0*U[i];
od;
U:=PureCubicalComplex(U);
inc:=RegularCWMap(K,U);

return ContractedComplex(Spin(inc));
end);
###############################################################################
###############################################################################

###############################################################################
###############################################################################
InstallGlobalFunction(SpunLinkComplement,
function(K)
local C,Y;
Y:=PureCubicalKnot(K);
Y:=PureComplexComplement(Y);
Y:=SpunAboutHyperplane(Y);
Y:=ContractedComplex(Y);
return Y;
end);
###############################################################################
###############################################################################

