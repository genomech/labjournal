# Конфигурация машины

Ubuntu, разумеется.

## Репозиторий Ubuntu

| Package | Description |
|:--------|:-----|
| ffmpeg | Нарезка видео (см. ниже) |
| gimp | Редактор растровых изображений |
| inkscape | Редактор векторных изображений |
| kile | GUI для TeX |
| libreoffice | Набор офисных программ |
| mc | Консольный файловый менеджер |
| python3 | Python 3 |
| qbittorrent | Торрент-клиент |
| texlive | Набор пакетов TeX |
| vim | Консольный текстовый редактор |

Биоинформатика:

| Package | Description |
|:--------|:-----|
| bcftools ||
| bowtie2 ||
| cutadapt ||
| fastqc ||
| hisat2 ||
| samtools ||
| vcftools ||

## Сторонние приложения

| Package | Description | Link |
|:--------|:-----|:-----|
| Arduino IDE | Программирование микроконтроллеров | [здесь](https://www.arduino.cc/en/Main/Software) |
| Qt | Фреймворк C++ | [здесь](https://www.qt.io/download-qt-installer)|
| SQLite Studio | GUI для редактирования СУБД | [здесь](https://sqlitestudio.pl/index.rvt?act=download) |
| Telegram | Мессенджер ||

## R

Установка:

```bash
ver=$(lsb_release -cs); echo "deb https://cloud.r-project.org/bin/linux/ubuntu "$ver"-cran35/" | sudo tee -a /etc/apt/sources.list

gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | sudo apt-key add -

sudo apt update

sudo apt-get install r-base r-base-dev
```

Установку пакетов *R* лучше всего производить через **sudo**.

## Прокси и Tor

**privoxy + tor**: локальный прокси-сервер с выходом в луковую сеть.

Установка:

```bash
sudo apt install tor tor-geoipdb privoxy
```

Если не устанавливается *privoxy*, нужно следовать указаниям на [официальном сайте](https://www.privoxy.org/).

Если не устанавливается *tor*:

```bash
ver=$(lsb_release -cs); echo "deb https://deb.torproject.org/torproject.org "$ver" main" | sudo tee -a /etc/apt/sources.list

curl https://deb.torproject.org/torproject.org/A3C4F0F979CAA22CDBA8F512EE8CBC9E886DDD89.asc | gpg --import

gpg --export A3C4F0F979CAA22CDBA8F512EE8CBC9E886DDD89 | sudo apt-key add -

sudo apt update

sudo apt install tor tor-geoipdb
```

Настройка:

```bash
sudo mv /etc/privoxy/config /etc/privoxy/config.backup

echo -e "forward-socks5 / 127.0.0.1:9050 .\nconfdir /etc/privoxy\nlogdir /var/log/privoxy\nactionsfile default.action\nactionsfile user.action\nfilterfile default.filter\nlogfile logfile\ndebug 4096\ndebug 8192\nuser-manual /usr/share/doc/privoxy/user-manual\nlisten-address 127.0.0.1:8118\ntoggle 1\nenable-remote-toggle 0\nenable-edit-actions 0\nenable-remote-http-toggle 0\nbuffer-limit 4096" | sudo tee /etc/privoxy/config

sudo service privoxy restart

sudo service tor restart
```

Далее можно использовать как прокси-сервер `127.0.0.1:9050`.
Работает как с торрент-клиентами, так и с *Telegram*.

**Дополнительно:** установка Tor Browser:

```bash
sudo apt install torbrowser-launcher
```

## ffmpeg

Скрипт для упрощения нарезки видео:

```bash
#!/bin/bash

IN_FILE="/media/avicenna/Dopamine_Python/GoT/GOT/Season 08/GOT.[S08E02].2xRu.En.[qqss44].mkv"
OUT_FILE="/dev/my/MyDocs/temp/jenny.mkv"
TIME_START="00:51:05"
TIME_END="00:53:01"
V_TRACK=0
A_TRACK=2

ffmpeg -i "$IN_FILE" -qscale 0 -map 0:v:$V_TRACK -map 0:a:$A_TRACK -ss $TIME_START -to $TIME_END "$OUT_FILE"
```
