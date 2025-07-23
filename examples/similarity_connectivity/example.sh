#!/bin/bash

#####################################
#									#
#	This variables you can change	#
#									#
#####################################

# Directory where assets are stored
export DNATCO_ASSETS_PATH="$HOME/gitlab-projects/library/libllka/assets"
# Name of the folder where the examples will be created
FOLDER="./Examples"
# Just prexif on console
PREFIX="[Example Script]:"

#####################################
#									#
#		Do not touch this!			#
#									#
#####################################

print(){
	local MESSAGE="$1"
	echo -e "$PREFIX $MESSAGE"
}

sourceFileExists() {
	local FILE="$1"
	if [ -f "$FILE" ]; then
		return 1
	else
		print "ERROR: $FILE not exists"
		exit
	fi
}

fileExists(){
    local FILE="$1"
    if [ -f "$FILE" ]; then
        return 0
    else
        return 1
    fi
}

if [ ! -d "$FOLDER" ]; then
	mkdir -p "$FOLDER"
	echo "$PREFIX Creating folder."
else
	echo "$PREFIX Folder already exists."
fi

copy(){
	local NAME="$1"
	if fileExists ./$FOLDER/$NAME; then
		print " "
	else
		cp ./$NAME ./$FOLDER/$NAME
	fi
	
}

cd $FOLDER
print " "
print "Creating all connectivity file"
print "Using: './similarity_connectivity -i 1bna.cif -c all_connectivity'"
../similarity_connectivity -i ../1bna.cif -c all_connectivity
print " "
print "Creating all similarity file"
print "Using: './similarity_connectivity -i 1bna.cif -s all_similarity'"
rm all_similarity.json 2>/dev/null
../similarity_connectivity -i ../1bna.cif -s all_similarity
print " "
print "Creating Connectivity first step"
print "Using: './similarity_connectivity -i 1bna.cif -c connectivity_first_step -x 1'"
rm connectivity_first_step.json 2>/dev/null
../similarity_connectivity -i ../1bna.cif -c connectivity_first_step -x 1
print " "
print "Creating Connectivity first step with NtC AA01"
print "Using: './similarity_connectivity -i 1bna.cif -c connectivity_first_step_ntc_AA01 -x 1 -n AA01'"
rm connectivity_first_step_ntc_AA01.json 2>/dev/null
../similarity_connectivity -i ../1bna.cif -c connectivity_first_step_ntc_AA01 -x 1 -n AA01
print " "
print "Creating Connectivity second step with NtC AA01 and prev NtC AA02"
print "Using: './similarity_connectivity -i 1bna.cif -c connectivity_second_step_ntc_AA01_prev_AA02 -x 2 -n AA01 -b AA02'"
rm connectivity_second_step_ntc_AA01_prev_AA02.json 2>/dev/null
../similarity_connectivity -i ../1bna.cif -c connectivity_second_step_ntc_AA01_prev_AA02 -x 2 -n AA01 -b AA02
print " "
print "Creating Connectivity second step with NtC AA01 and next NtC AA02 and distance cut"
print "Using: './similarity_connectivity -i 1bna.cif -c connectivity_second_step_ntc_AA01_next_AB02_distance_cut -x 2 -n AA01 -a AB02 -d 0.5 '"
rm connectivity_second_step_ntc_AA01_next_AB02_distance_cut.json 2>/dev/null
../similarity_connectivity -i ../1bna.cif -c connectivity_second_step_ntc_AA01_next_AB02_distance_cut -x 2 -n AA01 -a AB02 -d 0.5 
print " "
print "Creating Extended mfcif file"
print "Using './similarity_connectivity -i 1bna.cif -o 1bna_Extended'"
rm 1bna_Extended.cif 2>/dev/null
../similarity_connectivity -i ../1bna.cif -o 1bna_Extended
print " "
print "Creating Similarity for second step"
print "Using: './similarity_connectivity -i 1bna.cif -s similarity_second_step -x 2'"
rm similarity_second_step.json 2>/dev/null
../similarity_connectivity -i ../1bna.cif -s similarity_second_step -x 2
print " "
print "Creating Similarity for second step with RMSD cut 0.5"
print "Using: './similarity_connectivity -i 1bna.cif -s similarity_second_step_rmsd_cut -x 2 -r 0.5'"
rm similarity_second_step_rmsd_cut.json 2>/dev/null
../similarity_connectivity -i ../1bna.cif -s similarity_second_step_rmsd_cut -x 2 -r 0.5 
print " "
print "All done!"
