## Boomer Library, v0.03
#Tiny *BASH* library for bioinformatics.
# 
#Install:
# 
#```bash
#echo alias boomer=\"source $path_to_boomer/boomer.sh\" >> ~/.bashrc
#boomer
#```
# 
#Compile this doc:
# 
#```bash
#boomer; CompileDoc
#```

LIB_NAME="Boomer"

### Timestamp
#Make timestamp in format [hh:mm:ss].
# 
#```bash
#start_time=$(StartTime) 
#echo "$(Timestamp $start_time)"
#```

function StartTime() { date +%s >&1; }
function Timestamp() { local PASSED=$(echo $(date +%s) - $1 | bc); echo "["$(date -d@$PASSED -u +%H:%M:%S)"]" >&1; }

### Random Name
#Make random string `a-zA-Z0-9` of certain length.
# 
#```bash
#random_filename="$(RandomName $string_length)".txt
#```

function RandomName() { cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w $1 | head -n 1 >&1; }

### Process Killer
#Kill all stopped processes and clear the memory.
# 
#```bash
#Killer
#```

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

### Logo
#Make simple but pretty logo.
# 
#```bash
#Logo "$program_name"
#```

function Logo() { echo; echo "------ "$1" ------"; echo; echo "Now: "$(date +'%Y-%m-%d %H:%M:%S')""; echo; Killer; }

### Seal Folder
#Count checksum of the folder, write it to `all.md5`, then make folder read-only.
# 
#```bash
#Seal "$folder_path"
#```

function Seal() {
	local MODULE_NAME="Seal"
	start_time=$(StartTime);
	md5deep -lr $1/$2 > $1/all.md5
	chmod -R 555 $1
	echo ""$LIB_NAME"."$MODULE_NAME": Folder "$1" is sealed "$(Timestamp $start_time)""
}

### Read File
#Cat any file, even bzipped and gzipped.
# 
#```bash
#Read "$filename" | command1
#```

function Read() { ( zcat $1 || bzcat $1 || cat $1 ) >&1 2> /dev/null; }

### Count Lines
#Count lines in any file, even bzipped and gzipped.
# 
#```bash
#lines_num=$(CountLines $filename)
#```

function CountLines() { Read $1 | wc -l >&1; }

### Bam Index
#Make index of *BAM* files.
# 
#```bash
#BamIndex $folder_name/*.bam
#```

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

### Create FIFO pipe
#Create FIFO pipe with random name.
# 
#```bash
#pipe_name=$(CreateFIFO)
#cat $pipe_name | command1 & process_with_stdout | tee $pipe_name | command2
#rm $pipe_name
#```

function CreateFIFO() { 
	
	local NEW_UID="/tmp/"$(RandomName 32)".pipe"
	mkfifo $NEW_UID
	echo $NEW_UID >&1
}

### Filename Procedures
#Return absolute dir path, extension and basename of file.
# 
#```bash
#dir=$(FileDir $filename)
#extension=$(FileExt $filename)
#basename=$(FileBase $filename)
#```

function FileDir() { local full=$(realpath $1); dirname $full >&1; }
function FileExt() { local full=$(realpath $1); local basename=$(basename $full); echo ${basename#*.} >&1; }
function FileBase() { local full=$(realpath $1); basename $full .$(FileExt $1) >&1; }

function CompileDoc() ( cat boomer.sh | grep -oP "(?<=^#).*" > boomer.md; echo ""$LIB_NAME": Doc is ready." )
