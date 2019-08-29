import subprocess

print(subprocess.run(["md5sum", "./blister.py"], capture_output=True).stdout)
