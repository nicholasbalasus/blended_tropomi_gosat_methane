#!/bin/csh -f
# @(#) dindex

# Get the storage directory from the first argument
set StorageDir=$1

#-----------------------------------------------------------------
# Set path
#-----------------------------------------------------------------
set Input_Dir   = "${StorageDir}/oversample/input/"
set Output_Dir  = "${StorageDir}/oversample/output/"

#-----------------------------------------------------------------
# Output resolution you want
#-----------------------------------------------------------------
set Res = 0.01

set filenames = (`find $Input_Dir -name "*.csv" -printf "%f\n"`)

foreach Input_Filename ($filenames)
#-----------------------------------------------------------------
# Set file names
#-----------------------------------------------------------------
set Input_Filename = $Input_Filename
set addition = "_oversampled"
set Output_Filename = "`echo $Input_Filename | sed 's/.csv//1'`$addition$Res.csv"
set inin = $Input_Dir$Input_Filename
echo $inin
echo $Output_Filename
#-----------------------------------------------------------------
# Call RegridPixel.x, and pass user inputs
#-----------------------------------------------------------------
./RegridPixels.x<<EOF
$Output_Dir
$inin
$Output_Filename
$Res
EOF

quit:
end

exit
