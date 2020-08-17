LIB_NAME="Boomer"

function StartTime() { date +%s >&1; }

function Timestamp() { local PASSED=$(echo $(date +%s) - $1 | bc); echo "["$(date -d@$PASSED -u +%H:%M:%S)"]" >&1; }

function RandomName() { cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w $1 | head -n 1 >&1; }

function Killer() {
	local MODULE_NAME="Killer"
	lst=$(jobs -p)
	if [[ "$lst" = "" ]]; then :
	else {
		kill -s KILL $lst
		echo ""$LIB_NAME"."$MODULE_NAME": Killed jobs: "$lst""
		echo
	} 
	fi
}

function Logo() { echo; echo "------ "$1" ------"; echo; echo "Now: "$(date +'%Y-%m-%d %H:%M:%S')""; echo; Killer; }

function Seal() {
	local MODULE_NAME="Seal"
	start_time=$(StartTime);
	md5deep -lr $1/$2 > $1/all.md5
	chmod -R 555 $1
	echo ""$LIB_NAME"."$MODULE_NAME": Folder "$1" is sealed "$(Timestamp $start_time)""
}

function Read() { ( zcat $1 || bzcat $1 || cat $1 ) >&1 2> /dev/null; }

function CountLines() { Read $1 | wc -l >&1; }

function BamIndex() {
	local MODULE_NAME="BamIndex"
	local start_time=$(StartTime);
	local files=""
	for INFILE in "$@"; do
		samtools index $INFILE
		files+=""$'\t'""$INFILE""$'\n'""
	done
	echo ""$LIB_NAME"."$MODULE_NAME": Indexed files "$(Timestamp $start_time)":"
	echo "$files"
}

function CreateFIFO() { 
	local NEW_UID="/tmp/"$(RandomName 32)".pipe"
	mkfifo $NEW_UID
	echo $NEW_UID >&1
}


function FileDir() { local full=$(realpath $1); dirname $full >&1; }

function FileExt() { local full=$(realpath $1); local basename=$(basename $full); echo ${basename#*.} >&1; }

function FileBase() { local full=$(realpath $1); basename $full .$(FileExt $1) >&1; }

function Row() { awk "NR == "$1" { print } NR > "$1" { exit }" $2 >&1; }

function RealList() { for file in $(ls -1); do echo $(realpath $file) >&1; done }
