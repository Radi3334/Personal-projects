import matplotlib.pyplot as plt
import numpy as np
import pygame
import pandas as pd

#graphing boundries
x_data = np.linspace(0, 1, 100)
gamma=1.4
t_end=10

#creating object for each cell
class volume:  #assuming a volume
       def __init__(self, x, y, column, row, inside):  #1 is left, or top, 2 is right or bottom
        self.x_position=x
        self.y_position=y
        self.row=row
        self.column=column
        self.inside=inside

count=50
tolerance=1/count
filtered=[]
times=[]
#set initial conditions
for i in range(1,count):
    for j in range(1,count):
        if i>5*count/6:
            object=volume(x=i/count, y=j/count, row=j, column=i, inside=True)
            filtered.append(object)
        else:
            object=volume(x=i/count, y=j/count, row=j, column=i, inside=True)
            filtered.append(object)

pygame.init()
# Set up the display
width, height = 800, 800
screen = pygame.display.set_mode((width, height))
pygame.display.set_caption("Plot Example")
WHITE = (255, 255, 255)
RED = (255, 0, 0)

def map_to_screen(x, y, x_min=0, x_max=1, y_min=0, y_max=1):
    """Convert data coordinates (x, y) to screen coordinates (px, py)"""
    px = int((x - x_min) / (x_max - x_min) * width)
    py = int(height - (y - y_min) / (y_max - y_min) * height)  # y=0 at bottom
    return px, py

dx,dy=1/count, 1/count
running=True
grid = {(v.row, v.column): v for v in filtered}
while running:
    screen.fill(WHITE)
    rect_width = width/count
    
    temporary_storage={}
    for i in filtered:
            cell_east  = grid.get((i.row, i.column+1),None)
            cell_west  = grid.get((i.row, i.column-1),None)
            cell_north = grid.get((i.row-1, i.column),None)
            cell_south = grid.get((i.row+1, i.column),None)

            #boundry conditions
            if cell_west is None:
                hello=0
                
                
            #elif  cell_east is None:
  
            #elif cell_north is None:

            #elif cell_south is None:

            #else: 

            
    pygame.display.flip()

    # Event handling
for event in pygame.event.get():
    if event.type == pygame.QUIT:
        running=False












