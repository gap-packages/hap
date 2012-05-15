#############################################################################
##
##  PackageInfo.g                HAP Package                Graham Ellis 
##
#############################################################################

SetPackageInfo( rec(

  PackageName := "HAP",
  Subtitle  := "Homological Algebra Programming",
  Version := "1.10.1",
  Date    := "16/03/2012",
  ArchiveURL 
          := "http://hamilton.nuigalway.ie/Hap/hap1.10.1",
  ArchiveFormats 
          := ".tar.gz",


  Persons := [ 
    rec( 
      LastName      := "Ellis",
      FirstNames    := "Graham",
      IsAuthor      := true,
      IsMaintainer  := true,
      Email         := "graham.ellis@nuigalway.ie",
      WWWHome       := "http://hamilton.nuigalway.ie",
      PostalAddress := Concatenation( [
                         "Graham Ellis\n",
                         "Mathematics Department\n",
                         "NUI Galway\n",
                         "Galway\n",
                         "Ireland" ] ),
      Place         := "Galway",
      Institution   := "National University of Ireland, Galway"
    )
  ],  

  Status  := "accepted",
  CommunicatedBy 
          := "Derek Holt (Warwick)",
  AcceptDate 
          := "03/2006",

  README_URL := "http://hamilton.nuigalway.ie/Hap/README.HAP",
  PackageInfoURL := "http://hamilton.nuigalway.ie/Hap/PackageInfo.g",

  AbstractHTML := 
    "This package provides some functions for group cohomology. ",

  PackageWWWHome := "http://hamilton.nuigalway.ie/Hap/www",
                  
  PackageDoc := rec(
    BookName  := "HAP",
    ArchiveURLSubset := ["doc", "www"],
    HTMLStart := "www/index.html",
    PDFFile   := "doc/manual.pdf",
    SixFile   := "doc/manual.six",
    LongTitle := "Homological Algebra Programming Package",
    Autoload := true 
  ),


  Dependencies := rec(
    GAP := ">= 4.4.7",
    NeededOtherPackages := [],
    SuggestedOtherPackages := [[ "polycyclic", ">=1.1" ],
    		 	       [ "nq",         ">=1.1" ],
                               ["nql",         ">=0.0" ],
                               ["homology",    ">=0.0"   ],
			       [ "aclib",      ">=1.1" ],
			       [ "edim",      ">=1.2.2" ],
			       [ "gapdoc",      ">=0.0" ],
			       [ "singular", ">=06.07.23" ],
                               [ "congruence", ">=0,0" ]
			      ],

    ExternalConditions := [["Some optional functions require Polymake software",
    "http://www.math.tu-berlin.de/polymake/"],
    ["Some optional functions require Graphviz software",
    "http://www.graphviz.org/"],
     ["One optional function requires the Simplicial Homology GAP package",
         "http://www.cis.udel.edu/~dumas"]
    ]
  ),

AvailabilityTest := ReturnTrue,

BannerString     := Concatenation( "Loading HAP ",
                            String( ~.Version ), " ...\n" ),

Autoload := false,

TestFile := "test/hap.tst",

Keywords := [ "homology", "cohomology", "resolution", "homotopy group", 
"module of identities" ]

));

