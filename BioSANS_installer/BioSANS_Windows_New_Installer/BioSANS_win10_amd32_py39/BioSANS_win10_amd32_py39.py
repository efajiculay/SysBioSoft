import sys
from os import getcwd, path, system
from subprocess import Popen, PIPE
from wget import download
import requests

try:
	gcwd = sys._MEIPASS
except:	
	gcwd = getcwd()


	
session = requests.Session()
#url = "https://drive.google.com/file/d/12f9XBttlVJsctNrquJNRNhfE7iRXZsYT/view?usp=sharing"
#url = "https://drive.google.com/uc?export=download&id=12f9XBttlVJsctNrquJNRNhfE7iRXZsYT"


def download_file_from_google_drive(id, destination):
	URL = "https://docs.google.com/uc?export=download"

	session = requests.Session()

	response = session.get(URL, params = { 'id' : id }, stream = True)
	token = get_confirm_token(response)

	if token:
		params = { 'id' : id, 'confirm' : token }
		response = session.get(URL, params = params, stream = True)

	save_response_content(response, destination)	

def get_confirm_token(response):
	for key, value in response.cookies.items():
		if key.startswith('download_warning'):
			return value

	return None

def save_response_content(response, destination):
	CHUNK_SIZE = 32768
	with open(destination, "wb") as f:
		count = 0
		for chunk in response.iter_content(CHUNK_SIZE):
			if chunk: # filter out keep-alive new chunks
				f.write(chunk)
				count = count + 1
				print(" Download progress :",round(100*CHUNK_SIZE*count/134512640,3),"%",end="\r")

if __name__ == "__main__":
	file_id = '12f9XBttlVJsctNrquJNRNhfE7iRXZsYT'
	destination = path.join(gcwd,"BioSANS_py39_amd64_pre_installed.exe")
	
	print("\nDownloading BioSANS Installer\n")
	download_file_from_google_drive(file_id, destination) 
	print("\n\nInstallation will start in a few seconds\n")
	system("start /B "+str(destination).replace("\\","/"))
