#!/bin/bash

#####################################
#									#
#	This variables you can change	#
#									#
#####################################

# Directory where assets are stored
export DNATCO_ASSETS_PATH="$HOME/gitlab-projects/library/libllka/assets"
# Just prexif on console
PREFIX="[Testing Script]:"

#####################################
#									#
#		Do not touch this!			#
#									#
#####################################

RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
AQUA="\e[36m"
CRESET="\e[0m"
SUCCESS=0
WARNS=0
FAILS=0
SUCCESC=0
WARNC=0
FAILC=0
TOTALS=0
TOTALC=0

print(){
	local MESSAGE="$1"
	echo -e "$PREFIX $MESSAGE"
}

if [ ! -d "$DNATCO_ASSETS_PATH" ]; then
  print "Path '$DNATCO_ASSETS_PATH' not configured!"
  exit 1
fi

assetsExists(){
	local NAME="$1"
	if [ ! -f "$DNATCO_ASSETS_PATH/$NAME" ]; then
		print ""
		print "${RED}ERROR${CRESET}: File $NAME in $DNATCO_ASSETS_PATH not found!"
		print ""
		exit 1
	fi
}

assetsExists "clusters.csv"
assetsExists "confal_percentiles.csv"
assetsExists "confals.csv"
assetsExists "golden_steps.csv"
assetsExists "nu_angles.csv"

fileExists() {
	local FILE="$1"
	if [ -f "$FILE" ]; then
		return 1
	else
		print ""
		print "${RED}ERROR${CRESET}: $FILE not exists"
		print ""
		exit
	fi
}

fileExists "similarity_connectivity"
fileExists "1bna.cif"

testSim(){
	((TOTALS++))
	((TOTALS++))
	local refRMSD="$1"
	local RMSD="$2"
	local refED="$3"
	local ED="$4"
	local refRMSDTrim="${refRMSD:0:-1}"
	local RMSDPart="${refRMSDTrim#*.}"
	local RMSDLen="${#RMSDPart}"
	local RMSD_int="${RMSD%%.*}"
	local RMSD_dec="${RMSD#*.}"
	RMSD_dec="${RMSD_dec:0:RMSDLen}"
	local RMSDTrim="${RMSD_int}.${RMSD_dec}"
	local refEDTrim="${refED:0:-1}"
	local EDPart="${refEDTrim#*.}"
	local EDLen="${#EDPart}"
	local ED_int="${ED%%.*}"
	local ED_dec="${RMSD#*.}"
	ED_dec="${ED_dec:0:EDLen}"
	local EDTrim="${RMSD_int}.${RMSD_dec}"
    if [ "$refRMSD" = "$RMSD" ]; then
        print "refRMSD: ${AQUA}$refRMSD${CRESET}   RMSD: ${AQUA}$RMSD${CRESET}    Result: ${GREEN}Correct${CRESET}"
		((SUCCESS++))
	elif [ "$refRMSDTrim" = "$RMSDTrim" ]; then
		print "refRMSD: ${YELLOW}$refRMSD${CRESET}   RMSD: ${YELLOW}$RMSD${CRESET}    Result: ${YELLOW}Warning${CRESET}"
		((WARNS++))
    else
        print "refRMSD: ${RED}$refRMSD${CRESET}   RMSD: ${RED}$RMSD${CRESET}    Result: ${RED}Failed${CRESET}"
		((FAILS++))
    fi
	
	if [ "$refED" = "$ED" ]; then
        print "refED: ${AQUA}$refED${CRESET}   ED: ${AQUA}$ED${CRESET}    Result: ${GREEN}Correct${CRESET}"
		((SUCCESS++))
	elif [ "$refEDTrim" = "$EDTrim" ]; then
		 print "refED: ${YELLOW}$refED${CRESET}   ED: ${YELLOW}$ED${CRESET}    Result: ${YELLOW}Warning${CRESET}"
		 ((WARNS++))
    else
        print "refED: ${RED}$refED${CRESET}   ED: ${RED}$ED${CRESET}    Result: ${RED}Failed${CRESET}"
		((FAILS++))
    fi
}

testOutputSimilarity(){
	local output="$1"
	local refRMSD="$2"
	local refED="$3"
	local struct="$4"
	local step="$5"
	local ntc="$6"
    local rmsd=$(echo "$output" | awk '{print $2}')
    local ed=$(echo "$output" | awk '{print $4}')
	print " "
	print "Testing structure: ${AQUA}$struct${CRESET} for step ${AQUA}$step${CRESET} for ntc ${AQUA}$ntc${CRESET}."
	testSim $refRMSD $rmsd $refED $ed
	print " "
}

testCon(){
	((TOTALC++))
	((TOTALC++))
	local refC5="$1"
	local C5="$2"
	local refO3="$3"
	local O3="$4"
	local refO3Trim="${refO3:0:-1}"
	local O3Part="${refO3Trim#*.}"
	local O3Len="${#O3Part}"
	local O3_int="${O3%%.*}"
	local O3_dec="${O3#*.}"
	O3_dec="${O3_dec:0:${O3Len}}"
	local O3Trim="${O3_int}.${O3_dec}"
	local refC5Trim="${refC5:0:-1}"
	local C5Part="${refC5Trim#*.}"
	local C5Len="${#C5Part}"
	local C5_int="${C5%%.*}"
	local C5_dec="${C5#*.}"
	C5_dec="${C5_dec:0:${C5Len}}"
	local C5Trim="${C5_int}.${C5_dec}"
    if [ "$refC5" = "$C5" ]; then
        print "refC5: ${AQUA}$refC5${CRESET}   C5: ${AQUA}$C5${CRESET}    Result: ${GREEN}Correct${CRESET}"
		((SUCCESC++))
	elif [ "$refC5Trim" = "$C5Trim" ]; then
		print "refC5: ${YELLOW}$refC5${CRESET}   C5: ${YELLOW}$C5${CRESET}    Result: ${YELLOW}Warning${CRESET}"
		((WARNC++))
    else
        print "refC5: ${RED}$refC5${CRESET}   C5: ${RED}$C5${CRESET}    Result: ${RED}Failed${CRESET}"
		((FAILC++))
    fi
	
	if [ "$refO3" = "$O3" ]; then
        print "refO3: ${AQUA}$refO3${CRESET}   O3: ${AQUA}$O3${CRESET}    Result: ${GREEN}Correct${CRESET}"
		((SUCCESC++))
	elif [ "$refO3Trim" = "$O3Trim" ]; then
		 print "refO3: ${YELLOW}$refO3${CRESET}   O3: ${YELLOW}$O3${CRESET}    Result: ${YELLOW}Warning${CRESET}"
		 ((WARNC++))
    else
        print "refO3: ${RED}$refO3${CRESET}   O3: ${RED}$O3${CRESET}    Result: ${RED}Failed${CRESET}"
		((FAILC++))
    fi
}

testOutputConnectivity() {
	local output="$1"
	#Prev

	local line1=$(echo "$output" | sed -n '1p')
	local refC5Prev="$2"
	local refO3Prev="$3"

	local C5Prev=$(echo "$line1" | awk '{print $2}')
	local O3Prev=$(echo "$line1" | awk '{print $4}')
	#Next
	local line2=$(echo "$output" | sed -n '2p')
	local refC5Next="$4"
	local refO3Next="$5"
	local C5Next=$(echo "$line2" | awk '{print $2}')
	local O3Next=$(echo "$line2" | awk '{print $4}')
	local struct="$6"
	local step="$7"
	local ntc="$8"
	local prevNtC="$9"
	local nextNtC="${10}"
	
	print "Testing prev. connectivity structure: ${AQUA}$struct${CRESET} for step ${AQUA}$step${CRESET} for ntc ${AQUA}$ntc${CRESET} for prevNtC ${AQUA}$prevNtC${CRESET}"
	testCon $refC5Prev $C5Prev $refO3Prev $O3Prev
	print " "
	print "Testing next connectivity structure: ${AQUA}$struct${CRESET} for step ${AQUA}$step${CRESET} for ntc ${AQUA}$ntc${CRESET} for prevNtC ${AQUA}$nextNtC${CRESET}"
	testCon $refC5Next $C5Next $refO3Next $O3Next
	print " "
	}
print " "
print "Testing similarity"
testOutputSimilarity "$(./similarity_connectivity -i 1bna.cif -t -s -x 1 -n BA01)" "0.390329" "77.5079" "1bna" "1" "BA01"
testOutputSimilarity "$(./similarity_connectivity -i 1bna.cif -t -s -x 10 -n AA01)" "1.12085" "273.687" "1bna" "10" "AA01"
testOutputSimilarity "$(./similarity_connectivity -i 1bna.cif -t -s -x 14 -n OP03)" "2.33127" "265.795" "1bna" "14" "OP03"
testOutputSimilarity "$(./similarity_connectivity -i 1bna.cif -t -s -x 17 -n AB1S)" "0.865051" "185.206" "1bna" "17" "AB1S"
testOutputSimilarity "$(./similarity_connectivity -i 1bna.cif -t -s -x 21 -n ZZS2)" "1.61601" "255.602" "1bna" "21" "ZZS2"
print " Testing similarity complete"
print " Result: "
print " $SUCCESS/$TOTALS Correct"
print " $WARNS/$TOTALS Warning"
print " $FAILS/$TOTALS Failed"
print ""
print "Testing connectivity"
print " "
testOutputConnectivity "$(./similarity_connectivity -t -c -x 10 -a BA13 -b BBS1 -i 1bna.cif -n BB07)" "0.720698" "0.0623711" "0.813698" "0.0879142" "1bna" "10" "BB07" "BBS1" "BA13"
testOutputConnectivity "$(./similarity_connectivity -t -c -x 2 -a AA01 -b ZZS1 -i 1bna.cif -n ZZS2)" "3.15064" "2.68763" "1.30261" "3.00979" "1bna" "2" "ZZS2" "ZZS1" "AA01"
testOutputConnectivity "$(./similarity_connectivity -t -c -x 4 -a OP02 -b AB1S -i 1bna.cif -n AA03)" "1.04119" "0.721056" "2.7636" "2.55017" "1bna" "4" "AA03" "AB1S" "OP02"
testOutputConnectivity "$(./similarity_connectivity -t -c -x 14 -a AB1S -b OP02 -i 1bna.cif -n AB2S)" "0.962" "1.77422" "1.32313" "0.649352" "1bna" "14" "AB2S" "OP02" "AB1S"
testOutputConnectivity "$(./similarity_connectivity -t -c -x 19 -a ZZS2 -b AA02 -i 1bna.cif -n OP04)" "0.777052" "2.50414" "1.47854" "4.15402" "1bna" "10" "OP04" "AA02" "ZZS2"

print " Testing connectivity complete"
print " Result: "
print " $SUCCESC/$TOTALC Correct"
print " $WARNC/$TOTALC Warning"
print " $FAILC/$TOTALC Failed"
print ""
print " Testing complete"
print " Total result: "
print " $((SUCCESC+SUCCESS))/$((TOTALC+TOTALS)) Correct"
print " $((WARNS+WARNC))/$((TOTALC+TOTALS)) Warning "
print " $((FAILS+FAILC))/$((TOTALC+TOTALS)) Failed"
print ""
