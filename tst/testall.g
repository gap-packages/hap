LoadPackage( "hap" );

dirs := Concatenation(
  DirectoriesPackageLibrary( "hap", "tst/testall" ),
  DirectoriesPackageLibrary( "hap", "tst/testall2" ),
  DirectoriesPackageLibrary( "hap", "tst/testextra" ),
  DirectoriesPackageLibrary( "hap", "tst/testextra2" )
);

if not GAPInfo.Architecture="x86_64-pc-linux-gnu-default32-kv8" then
    Append(dirs, DirectoriesPackageLibrary( "hap", "tst/testall3" ));
fi;

TestDirectory(dirs,
  rec(exitGAP     := true,
      testOptions := rec(compareFunction := "uptowhitespace") ) );

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
