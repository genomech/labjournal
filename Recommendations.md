# Мудрость биоинформатика

* [Архивация](#zip)
* [Pickle](#pickle)
* [Парсинг](#parsing)
* [Контрольные суммы](#checksum)
* [Таймстампы](#timestamp)
* [Blister IO](#blister_io)
* [Графики](#graph)

<a name="zip"></a>
## Архивация

Не архивировать промежуточные файлы (e.g., сохранять не в .fastq.gz, а просто в .fastq).
Скорость работы bash-скрипта возрастает в разы.

<a name="pickle"></a>
## Pickle

Всегда снимать [pickle](https://docs.python.org/3/library/pickle.html) при работе с большими библиотеками.
Поверь, это лучше, чем отрезать куски от файла, и быстрее, чем парсить файл каждый раз.
Но только pickles не снимать, обязательно дублировать в формате БД.

```python
import pickle

# dump
with open(dump_filename, 'wb') as f:
    pickle.dump(data, f)

# load
with open(dump_filename, 'rb') as f:
    data_new = pickle.load(f)
```

<a name="parsing"></a>
## Парсинг

Использовать для парсинга уже имеющиеся библиотеки.
Это круто по двум причинам: 1) не нужно изобретать велосипед; 2) Не нужно делать мусорные файлы head'ом для теста -- при правильной настройке функции отрезается нужное количество записей.

1. VCF: [PyVCF](https://pyvcf.readthedocs.io/en/latest/)
    
```python
import vcf

vcf_reader = vcf.Reader(open(vcf_filename, 'r'))

for record in vcf_reader:
    # do stuff
```

2. Sequence форматы, а также SAM/BAM: [BioPython](https://biopython.org/wiki/Documentation)

```python
from Bio import SeqIO

# input
with open(input_filename, "rU") as input_handle:
    for record in SeqIO.parse(input_handle, input_format):
        # do stuff

# output
with open(output_filename, "w") as output_handle:
    SeqIO.write(sequences, output_handle, output_format)
    
# convert
count = SeqIO.convert(input_filename, input_format, output_filename, output_format)
print(f"Converted {count} records from {input_format} to {output_format}", end='\n')
```

3. CSV: [pandas](https://pandas.pydata.org/pandas-docs/stable/)

```python
import pandas as pd

data = pd.read_csv(input_filename, sep=',', nrows=None)
data.to_csv(output_filename, sep='\t', index=False, mode='w')
```

<a name="checksum"></a>
## Контрольные суммы

Сохранять контрольные суммы больших файлов и всегда их проверять перед работой.

```
$ md5sum [file] > [file].md5
$ md5sum -c [file].md5
```

<a name="timestamp"></a>
## Таймстампы

```python
import time

start_time = time.time()

# do stuff

print(f"Stuff is done [%.2f sec]" % (time.time() - start_time), end="\n")
```

<a name="blister_io"></a>
## Blister IO

В bash это делается довольно просто.

```bash
#!/bin/bash

# dir
mkdir -p /dir/path

# files
for file in `ls -d *.*`
do
# do stuff
done
```

В питоне посложнее.
Вот блистерная библиотечка для работы с файлами.

```python
from PyQt5.QtCore import QFileInfo, QFile, QDir
import glob
import gzip
import bzip2

# Input file(s) handle. Can do file masks.

def blister_input(filenames):

    file_list = []
    fileinfo_list = []
    fileinfo_unreadable = []
    
    for filename in filenames:
        file_list += glob.glob(QFileInfo(filename).absoluteFilePath())
        
    file_list = list(set(file_list))
    file_list.sort()
    
    for _file in file_list:
        fileinfo = QFileInfo(_file)
        if fileinfo.isFile():
            if fileinfo.isReadable():
                fileinfo_list += [fileinfo.absoluteFilePath()]
            else:
                fileinfo_unreadable += [fileinfo.absoluteFilePath()]
    
    if fileinfo_unreadable:
        print(f"Blister: List of unreadable files (will not be processed):", end="\n")
        for fileinfo in fileinfo_unreadable: print(f"\t{fileinfo}", end="\n")
    
    if not fileinfo_list:
        print(f"Blister: No input files exist.", end="\n")
        return False
    
    print(f"Blister: List of input files:", end="\n")
    for fileinfo in fileinfo_list: print(f"\t{fileinfo}", end="\n")
    return fileinfo_list

# Output file handle.

def blister_output(filename, output_dir, mod, suffix, rewrite=True):
    
    if mod != "": mod = "_[" + mod + "]"
    
    if (suffix != ""): suffix = "." + suffix
        
    fileinfo_old = QFileInfo(filename)
    fileinfo = QFileInfo(QDir(output_dir), fileinfo_old.baseName() + mod + suffix)
    
    if (fileinfo.exists() and (not fileinfo.isFile())):
        print(f"Blister: This path is a dir:\n\t{fileinfo.absoluteFilePath()}", end="\n")
        return False
    if ((fileinfo.exists() and (not fileinfo.isWritable())) or ((not fileinfo.exists()) and (not QFileInfo(fileinfo.absolutePath()).permission(QFile.WriteUser)))):
        print(f"Blister: Writing this file is not permitted:\n\t{fileinfo.absoluteFilePath()}", end="\n")
        return False
    if (fileinfo.exists() and (rewrite == False)):
        fileinfo = QFileInfo(QDir(output_dir), fileinfo_old.baseName()+ "_" + str(int(time.time()) % 100000) + suffix)
        print(f"Blister: File to write already exists [rewriting is forbidden]. It will be renamed:\n\t{fileinfo_old.absoluteFilePath()} --> {fileinfo.absoluteFilePath()}", end="\n")
    
    return fileinfo.absoluteFilePath()

# Output dir handle.

def blister_dir(_dir, create=True):
    
    dir_info = QFileInfo(_dir)
    
    if (dir_info.exists() and (not dir_info.permission(QFile.WriteUser))):
        print(f"Blister: Writing in this dir is not permitted:\n\t{dir_info.absoluteFilePath()}", end="\n")
        return False
    if ((not dir_info.exists()) and (create == False)):
        print(f"Blister: This dir does not exist [creating new is forbidden]:\n\t{dir_info.absoluteFilePath()}", end="\n")
        return False
    if ((not dir_info.exists()) and (create == True)):
        result = QDir.mkpath(QDir(), dir_info.absoluteFilePath())
        if not result:
            print(f"Blister: Creating new dir was failed:\n\t{dir_info.absoluteFilePath()}", end="\n")
            return False
        else:
            print(f"Blister: New dir was created:\n\t{dir_info.absoluteFilePath()}", end="\n")
    
    return dir_info.absoluteFilePath()
    
# Gzip and bzip2 check

def blister_gzip_check(filename):
    GZIP_MAGIC_NUMBER = "1f8b"
    with open(filename, 'rb') as file_check:
        return file_check.read(2).hex() == GZIP_MAGIC_NUMBER

def blister_bzip2_check(filename):
    BZIP2_MAGIC_NUMBER = "425a68"
    with open(filename, 'rb') as file_check:
        return file_check.read(3).hex() == BZIP2_MAGIC_NUMBER
```

Usage:

```python
# input files
filenames = blister_input(input_file)
if not filenames: return 1

# output_dir
_dir = blister_dir(namespace.output, create=True)
if not _dir: return 1

# output file
new_filename = blister_output(filename, _dir, mod, suffix)
if not new_filename: return 1

# input_handle
input_handle = gzip.open(filename, "rt") if blister_gzip_check(filename) else (bz2.open(filename, "rt") if blister_bzip2_check(filename) else open(filename, "rU"))
```

<a name="graph"></a>
## Графики

[MatPlotLib](https://matplotlib.org/3.1.1/index.html) в помощь.
Хорош тем, что может писать вывод в SVG.
Все мы знаем, что графики в растре -- это б-гомерзость.
