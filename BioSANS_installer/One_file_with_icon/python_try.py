import sys

wd = sys._MEIPASS
file_path = os.path.join(wd,"read.txt")

with open(file_path,"r") as f:
	for x in f:
		print(x)
		
x = 1
while x!="Y":
	x = input("Do you want to exit? (Y/N)").strip().upper()