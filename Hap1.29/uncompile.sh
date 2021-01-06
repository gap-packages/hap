#The directory PKGDIR should be the directory in which GAP packages are stored on your computer. 

read version < version;
PKGDIR=/home/graham/pkg/Hap$version;

rm $PKGDIR/boolean;

rm -rf $PKGDIR/lib/*/Compiled
echo "COMPILED:=false;" > $PKGDIR/boolean;




