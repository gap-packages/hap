#You must set PKGDIR equal to the directory in which GAP packages are stored
# on your computer.

PKGDIR=/home/graham/pkg;

#You must set GACDIR equal to the directory in which the GAP compiler gac is
#stored on your computer.

GACDIR=/usr/share/gap-4.10.2/;



#####################################################################
#DON'T CHANGE ANYTHING BELOW
#####################################################################
read version < version;

LIB=$PKGDIR/Hap$version/lib;

rm $PKGDIR/Hap$version/boolean;
echo "COMPILED:=true;" > $PKGDIR/Hap$version/boolean;

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




