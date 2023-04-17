#You must set PKGDIR equal to the directory in which GAP packages are stored
# on your computer.

PKGDIR=$(dirname $(pwd));

#You must set GACDIR equal to the directory in which the GAP compiler gac is
#stored on your computer.

GACDIR=$(dirname $(which gac))

#####################################################################
#DON'T CHANGE ANYTHING BELOW
#####################################################################
read version < version;

LIB=$PKGDIR/hap-$version/lib;

rm $PKGDIR/hap-$version/boolean;
echo "COMPILED:=false;" > $PKGDIR/hap-$version/boolean;

rm -rf $PKGDIR/hap-$version/lib/*/Compiled








