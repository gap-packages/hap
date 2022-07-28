LoadPackage( "hap" );

if not GAPInfo.Architecture="x86_64-pc-linux-gnu-default32-kv8" then

if not EvalString(GAPInfo.Version{[3,4]})>=11 then

TestDirectory(
[DirectoriesPackageLibrary( "hap", "tst/testall" ),
DirectoriesPackageLibrary( "hap", "tst/testall2" ),
DirectoriesPackageLibrary( "hap", "tst/testextra" ),
DirectoriesPackageLibrary( "hap", "tst/testextra2" )],
#DirectoriesPackageLibrary( "hap", "tst/testall3" )],
rec(exitGAP     := true,
    testOptions := rec(compareFunction := "uptowhitespace") ) );

else

TestDirectory(
[DirectoriesPackageLibrary( "hap", "tst/testall" ),
DirectoriesPackageLibrary( "hap", "tst/testall2" ),
DirectoriesPackageLibrary( "hap", "tst/testextra" ),
DirectoriesPackageLibrary( "hap", "tst/testextra2" ),
DirectoriesPackageLibrary( "hap", "tst/testall3" ),
DirectoriesPackageLibrary( "hap", "tst/testallV11" )],
rec(exitGAP     := true,
    testOptions := rec(compareFunction := "uptowhitespace") ) );
fi;
else

TestDirectory(
[DirectoriesPackageLibrary( "hap", "tst/testall" ),
DirectoriesPackageLibrary( "hap", "tst/testall2" ),
DirectoriesPackageLibrary( "hap", "tst/testextra" ),
DirectoriesPackageLibrary( "hap", "tst/testextra2" )],
#DirectoriesPackageLibrary( "hap", "tst/testall3" )],
rec(exitGAP     := true,
    testOptions := rec(compareFunction := "uptowhitespace") ) );
fi;

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
