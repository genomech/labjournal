from blister import *

from PyQt5.QtCore import QCryptographicHash, QByteArray

def blister_thread_id(item):
	_hash = QCryptographicHash(QCryptographicHash.Sha256);
	_hash.addData(QByteArray(bytes(item, 'utf-8')))
	result = _hash.result().toHex().__str__()[-5:-1]
	
	print(result)

input_filenames = blister_input(["./*.*"])

for item in input_filenames:
	blister_thread_id(item)
