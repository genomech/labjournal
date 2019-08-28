from PyQt5.QtCore import QFileInfo, QFile, QDir
import glob
import gzip
import bzip2

# --------- BLISTER IO ---------

# blister_input(filenames): return list of full filepaths to input (can handle mask) of False
# blister_output(filename, output_dir, mod, suffix, rewrite=True): return full filepath to output (output_dir/basename_[mod].suffix) or False
# blister_dir(_dir, create=True): return full path to output or False
# blister_gzip_check(filename): return True if file is gzipped, False otherwise
# blister_bzip2_check(filename): return True if file is bzipped, False otherwise
# blister_input_open(filename, mode): return input handle whether file is zipped or not

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
    
def blister_gzip_check(filename):
    GZIP_MAGIC_NUMBER = "1f8b"
    with open(filename, 'rb') as file_check:
        return file_check.read(2).hex() == GZIP_MAGIC_NUMBER

def blister_bzip2_check(filename):
    BZIP2_MAGIC_NUMBER = "425a68"
    with open(filename, 'rb') as file_check:
        return file_check.read(3).hex() == BZIP2_MAGIC_NUMBER

def blister_input_open(filename, mode):
    return gzip.open(filename, mode) if blister_gzip_check(filename) else (bz2.open(filename, mode) if blister_bzip2_check(filename) else open(filename, mode))
