LoadPackage( "hap" );

TestDirectory(
[DirectoriesPackageLibrary( "hap", "tst/testall3" ),
DirectoriesPackageLibrary( "hap", "tst/testallV11" ),
DirectoriesPackageLibrary( "hap", "tst/testextra" ),
DirectoriesPackageLibrary( "hap", "tst/testextra" ),
DirectoriesPackageLibrary( "hap", "tst/testall" ),
DirectoriesPackageLibrary( "hap", "tst/testall2" ),
DirectoriesPackageLibrary( "hap", "tst/testextraextra" )],
rec(exitGAP:=true,
                   compareFunction := "uptowhitespace") ) );

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
