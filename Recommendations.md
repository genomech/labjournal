# Основные рекомендации 

1. Всегда снимать [pickle](https://docs.python.org/3/library/pickle.html) при работе с большими библиотеками.
Поверь, это лучше, чем отрезать куски от файла, и быстрее, чем парсить файл каждый раз.
    * VCF: [PyVCF](https://pyvcf.readthedocs.io/en/latest/)
    * Sequence форматы: [BioPython](https://biopython.org/wiki/Documentation)
