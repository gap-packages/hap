#The directory PKGDIR should be the directory in which GAP packages are stored on your computer. 

PKGDIR=/home/graham/pkg;

rm $PKGDIR/Hap1.19/boolean;

rm -rf $PKGDIR/Hap1.19/lib/*/Compiled
echo "COMPILED:=false;" > $PKGDIR/Hap1.19/boolean;




