#Run this script whenever the package version number has increased. It 
#first updates the .xml manual files, and then changes the version number 
#in several www files and in the README.HAP file.

gap -l '/home/graham/;' -A -r ./doc/updateUndocumented.gi;
gap -l '/home/graham/;' -r ./doc/updateIndex.gi;
./doc/rd.sh
