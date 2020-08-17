# Boomer Library, v0.03
Tiny *BASH* library for bioinformatics.
 
Install:
 
```bash
echo alias boomer=\"source $path_to_boomer/boomer.sh\" >> ~/.bashrc
boomer
```
 
Compile this doc:
 
```bash
boomer; CompileDoc
```
## Timestamp
Make timestamp in format [hh:mm:ss].
 
```bash
start_time=$(StartTime) 
echo "$(Timestamp $start_time)"
```
## Random Name
Make random string `a-zA-Z0-9` of certain length.
 
```bash
random_filename="$(RandomName $string_length)".txt
```
## Process Killer
Kill all stopped processes and clear the memory.
 
```bash
Killer
```
## Logo
Make simple but pretty logo.
 
```bash
Logo "$program_name"
```
## Seal Folder
Count checksum of the folder, write it to `all.md5`, then make folder read-only.
 
```bash
Seal "$folder_path"
```
## Read File
Cat any file, even bzipped and gzipped.
 
```bash
Read "$filename" | command1
```
## Count Lines
Count lines in any file, even bzipped and gzipped.
 
```bash
lines_num=$(CountLines $filename)
```
## Bam Index
Make index of *BAM* files.
 
```bash
BamIndex $folder_name/*.bam
```
## Create FIFO pipe
Create FIFO pipe with random name.
 
```bash
pipe_name=$(CreateFIFO)
cat $pipe_name | command1 & process_with_stdout | tee $pipe_name | command2
rm $pipe_name
```
## Filename Procedures
Return absolute dir path, extension and basename of file.
 
```bash
dir=$(FileDir $filename)
extension=$(FileExt $filename)
basename=$(FileBase $filename)
```
## Row from file
Return row by number.
```bash
Row $number $file
```
