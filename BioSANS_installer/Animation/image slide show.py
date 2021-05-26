import tkinter as tk
from itertools import cycle
from PIL import ImageTk, Image

images = ["ffa1c83.jpg", "ffa4600.jpg", "faa4149.jpg", "f099e64.png"]
photos = cycle(ImageTk.PhotoImage(Image.open(image)) for image in images)