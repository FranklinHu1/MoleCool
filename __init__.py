#Main file for running the application

from tkinter import *
from LewisStructureGen import *
import os
import os.path

def mainFunc():
    windowPopup = Tk()
    windowPopup.wm_title("Welcome!")
    windowPopup.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    windowPopup.resizable(width = False, height = False) #Cannot resize init screen
    label = Label(windowPopup, text = 'Please select your preferred screen size, width x height', wraplength = 400, font = 'Arial 12')
    label.pack(side = 'top', fill = X, pady = 10)
    B1 = Button(windowPopup, text = '500 x 500', command = lambda:[windowPopup.destroy(), run(500, 500)])
    B2 = Button(windowPopup, text = '700 x 700', command = lambda: [windowPopup.destroy(), run(700, 700)])
    B3 = Button(windowPopup, text = '1000 x 1000', command = lambda: [windowPopup.destroy(), run(1000, 1000)])
    B1.pack(fill = X), B2.pack(fill = X), B3.pack(fill = X)
    windowPopup.mainloop()

mainFunc()
