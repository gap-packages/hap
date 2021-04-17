LoadPackage( "hap" );

TestDirectory(
DirectoriesPackageLibrary( "hap", "tst/testextra" ),
DirectoriesPackageLibrary( "hap", "tst/testall" ),
DirectoriesPackageLibrary( "hap", "tst/testall2" ),
DirectoriesPackageLibrary( "hap", "tst/testextraextra" )
testOptions := rec(exitGAP:=true,
                   compareFunction := "uptowhitespace") ) );

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
