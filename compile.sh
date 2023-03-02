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
echo "COMPILED:=true;" > $PKGDIR/hap-$version/boolean;

#$GACDIR/gac -d $LIB/CompiledGAP/*.c;
#mkdir  $LIB/CompiledGAP/Compiled;
#mv *.so $LIB/CompiledGAP/Compiled/;

for dir in $(ls -d $LIB/*/); do
mkdir $dir/Compiled;
for file in $(ls $dir*.gi); do
$GACDIR/gac -d $file;
mv *.so $dir/Compiled;
rm *.la.la;
done;
done;




