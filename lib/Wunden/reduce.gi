#Iteratively simplifies a G-complex L by repeatedly finding a cell e in some
#dimension n whose k-star is
#   (i)  G-balanced      -- IsGBalanced(L,e,n,k) = true, and
#   (ii) a disk           -- HomeomorphicToDisk(L,e,n,k) = true,
#and replacing that star by a single new cell via Replace2(L,e,n,k).
#
#The search restarts from dimension 0 after every successful replacement,
#since Replace2 renumbers cells across several dimensions at once and it is
#not safe to keep old indices around. This continues until a full pass over
#the whole complex finds nothing left to replace.
#

InstallGlobalFunction(ReduceGComplex, function(L)
local CurrentTopDim, Progress, n, e, k, dim, kmax;
 
# L:=arg[1];
# if Length(arg)>1 and arg[2]>0 then
#     MaxK:=arg[2];
# else
#     MaxK:=infinity;
# fi;
 
#Highest dimension in which L currently has any cells.
CurrentTopDim:=function()
local d;
    d:=0;
    while L!.dimension(d+1) > 0 do
        d:=d+1;
    od;
    return d;
end;
 
repeat
    Progress:=false;
    n:=0;
    while n <= CurrentTopDim() and not Progress do
        dim:=L!.dimension(n);
        e:=1;
        while e <= dim and not Progress do
            kmax:=CurrentTopDim()-n;
            k:=1;
            while k <= kmax and not Progress do
                if IsGBalanced(L,e,n,k) and HomeomorphicToDisk(L,e,n,k) then
                    Info(InfoWunden, 1, "ReduceGComplex: replacing the ",k,
                          "-star of cell ",e," in dimension ",n);
                    L:=Replace2(L,e,n,k);
                    Progress:=true;
                fi;
                k:=k+1;
            od;
            e:=e+1;
        od;
        n:=n+1;
    od;
until not Progress;
 
return L;
end);
