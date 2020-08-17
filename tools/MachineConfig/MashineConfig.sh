# Big data folder
bigdata=~/.scepter

# Folder to move annovar/humandb to
annovar_humandb=$bigdata/annovar_humandb

# Create big data folders

mkdir -p $bigdata;
mkdir -p $annovar_humandb;

# Ubuntu repo packages

sudo apt install -y \
bcftools			# BCF \
bowtie2				# Aligner \
bwa					# Aligner \
curl				# FTP Downloader \
cutadapt			# Adapter Cut & Split \
fastqc				# FastQ Quality \
ffmpeg				# Video \
freebayes			# Variant Caller \
gimp				# Raster Images \
git					# Work with Repos \
git-lfs				# GATK Dependency \
hisat2				# RNA-Seq Aligner \
inkscape			# Vector \
kile				# LaTeX \
libreoffice			# Office \
libvcflib-tools		# VCF \
libxml-libxml-perl	# SRA-Tools Dependency \
mc					# File Manager \
picard-tools		# Bioinformatics Tools \
python3				# Python 3 \
python3-pip			# Python 3 Package Manager \
qbittorrent			# Torrent \
samtools			# SAM/BAM/CRAM \
texlive				# TeX \
texlive-xetex		# XeLaTeX \
texlive-fonts-recommended \
texlive-fonts-extra \
vcftools			# VCF \
vim;				# Text Processor

# python3

sudo pip3 install --user --upgrade pip;
sudo python3 -m pip install \
biopython			# FastQ, Fasta, etc. \
matplotlib			# Plots \
numpy				# Number Arrays\
pandas				# Tables \
psutil				# System Management \
pybedtools			# BED \
PyQt5				# Qt \
pysam				# SAM/BAM/CRAM \
pyvcf				# VCF \
tabulate			# Tables to MD or TeX or that stuff \
termcolor;			# Colorize console

# tor + privoxy

sudo apt install tor tor-geoipdb privoxy torbrowser-launcher;

sudo mv /etc/privoxy/config /etc/privoxy/config.backup;
echo -e "forward-socks5 / 127.0.0.1:9050 .\nconfdir /etc/privoxy\nlogdir /var/log/privoxy\nactionsfile default.action\nactionsfile user.action\nfilterfile default.filter\nlogfile logfile\ndebug 4096\ndebug 8192\nuser-manual /usr/share/doc/privoxy/user-manual\nlisten-address 127.0.0.1:8118\ntoggle 1\nenable-remote-toggle 0\nenable-edit-actions 0\nenable-remote-http-toggle 0\nbuffer-limit 4096" | sudo tee /etc/privoxy/config;
sudo service privoxy restart;
sudo service tor restart;

# R

ver=$(lsb_release -cs); echo "deb https://cloud.r-project.org/bin/linux/ubuntu "$ver"-cran35/" | sudo tee -a /etc/apt/sources.list > /dev/null;
gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 &&
gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | sudo apt-key add - &&
sudo apt update &&
sudo apt install r-base r-base-dev;

# PREBUILDS

# sra-tools

curl https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz | tar -xz -C ./prebuild/;

# IGV

target="~/IGV_2.8.2.zip";
curl https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_2.8.2.zip > $target &&
unzip -q -d ./prebuild/ $target &&
rm $target;

# Annovar

curl www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz | tar -xz -C ./prebuild/ &&
mv ./prebuild/annovar/humandb $annovar_humandb;

# Juicebox

mkdir -p ./prebuild/Juicebox/;
curl https://s3.amazonaws.com/hicfiles.tc4ga.com/public/Juicebox/Juicebox_1.11.08.jar > ./prebuild/Juicebox/Juicebox_1.11.08.jar;

# MODULES

# git submodule update --init --recursive;

# GATK

# https://github.com/broadinstitute/gatk
# TODO Через conda
