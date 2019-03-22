#You must set PKGDIR equal to the directory in which GAP packages are stored
# on your computer.

PKGDIR=/home/graham/pkg;

#You must set GACDIR equal to the directory in which the GAP compiler gac is
#stored on your computer.

GACDIR=/usr/local/lib/gap4r4/bin/i686-pc-linux-gnu-gcc;



#####################################################################
#DON'T CHANGE ANYTHING BELOW
#####################################################################


LIB=$PKGDIR/Hap1.19/lib;

rm $PKGDIR/Hap1.19/boolean;
echo "COMPILED:=true;" > $PKGDIR/Hap1.19/boolean;

$GACDIR/gac -d $LIB/CompiledGAP/*.c;
mkdir  $LIB/CompiledGAP/Compiled;
mv *.so $LIB/CompiledGAP/Compiled/;


$GACDIR/gac -d $LIB/ArtinCoxeter/*.gi; 
mkdir  $LIB/ArtinCoxeter/Compiled;
mv *.so $LIB/ArtinCoxeter/Compiled/;

$GACDIR/gac -d $LIB/FreeGmodules/*.gi;
mkdir $LIB/FreeGmodules/Compiled;
mv *.so $LIB/FreeGmodules/Compiled/;

$GACDIR/gac -d $LIB/Functors/*.gi;
mkdir $LIB/Functors/Compiled;
mv *.so $LIB/Functors/Compiled/;

$GACDIR/gac -d $LIB/Homology/*.gi;
mkdir $LIB/Homology/Compiled;
mv *.so $LIB/Homology/Compiled/;

$GACDIR/gac -d $LIB/NonabelianTensor/*.gi;
mkdir $LIB/NonabelianTensor/Compiled;
mv *.so $LIB/NonabelianTensor/Compiled/;

$GACDIR/gac -d $LIB/Perturbations/*.gi;
mkdir $LIB/Perturbations/*.gi;
mkdir $LIB/Perturbations/Compiled/;
mv *.so $LIB/Perturbations/Compiled/;

$GACDIR/gac -d $LIB/Polycyclic/*.gi;
mkdir $LIB/Polycyclic/Compiled/;
mv *.so $LIB/Polycyclic/Compiled/;

$GACDIR/gac -d $LIB/Polymake/*.gi;
mkdir  $LIB/Polymake/Compiled/;
mv *.so $LIB/Polymake/Compiled/;

$GACDIR/gac -d $LIB/Resolutions/*.gi;
mkdir $LIB/Resolutions/Compiled/;
mv *.so $LIB/Resolutions/Compiled/;

$GACDIR/gac -d $LIB/ResolutionsModP/*.gi;
mkdir $LIB/ResolutionsModP/Compiled/;
mv *.so $LIB/ResolutionsModP/Compiled/;

$GACDIR/gac -d $LIB/Rings/*.gi;
mkdir $LIB/Rings/Compiled/;
mv *.so $LIB/Rings/Compiled/;

$GACDIR/gac -d $LIB/ModPRings/*.gi;
mkdir $LIB/ModPRings/Compiled/;
mv *.so $LIB/ModPRings/Compiled/;

$GACDIR/gac -d $LIB/FpGmodules/*.gi;
mkdir $LIB/FpGmodules/Compiled/;
mv *.so $LIB/FpGmodules/Compiled/;

$GACDIR/gac -d $LIB/PolyComplexes/*.gi;
mkdir $LIB/PolyComplexes/Compiled/;
mv *.so $LIB/PolyComplexes/Compiled/;

$GACDIR/gac -d $LIB/SimplicialGroups/*.gi;
mkdir $LIB/SimplicialGroups/Compiled/;
mv *.so $LIB/SimplicialGroups/Compiled/;

$GACDIR/gac -d $LIB/RegularCWComplexes/*.gi;
mkdir $LIB/RegularCWComplexes/Compiled/;
mv *.so $LIB/RegularCWComplexes/Compiled/;


$GACDIR/gac -d $LIB/GraphsOfGroups/*.gi;
mkdir $LIB/GraphsOfGroups/Compiled/;
mv *.so $LIB/GraphsOfGroups/Compiled/;
