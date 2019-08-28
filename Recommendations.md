# Мудрость биоинформатика

* [Архивация](#zip)
* [subprocess](#subprocess)
* [Pickle](#pickle)
* [Парсинг](#parsing)
* [Контрольные суммы](#checksum)
* [Таймстампы](#timestamp)
* [IO](#io)
* [Графики](#graph)

<a name="zip"></a>
## Архивация

Не архивировать промежуточные файлы (e.g., сохранять не в .fastq.gz, а просто в .fastq).
Скорость работы bash-скрипта возрастает в разы.

<a name="subprocess"></a>
## subprocess

Использовать [subprocess](https://docs.python.org/3/library/subprocess.html) вместо bash-скриптов.
Лучше отчётность, проще следить.
Также помогает использовать уже готовые утилиты в скрипте.

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

<a name="io"></a>
## IO

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
В скриптах есть блистерная библиотечка для работы с файлами.

<a name="graph"></a>
## Графики

[MatPlotLib](https://matplotlib.org/3.1.1/index.html) в помощь.
Хорош тем, что может писать вывод в SVG.
Все мы знаем, что графики в растре -- это б-гомерзость.
