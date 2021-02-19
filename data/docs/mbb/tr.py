from google_trans_new import google_translator  
import time

translator = google_translator()

new_text = []

for item in open("./textryba", 'rt').readlines():
	item = item[:-1]
	print(f"'{item}'")
	line = translator.translate(item, lang_src='en', lang_tgt='ru') if item != "" else item
	print(line)
	new_text += [ line + '\n' ]
	time.sleep(1)
	#break

open("./textryba.translate", 'wt').writelines(new_text)
