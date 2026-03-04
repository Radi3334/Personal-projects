import matplotlib.pyplot as plt
import numpy as np
import pygame
import pandas as pd

#graphing boundries
x_data = np.linspace(0.1, 1, 100)
y_data_1= x_data/x_data
y_data_2=0*x_data
t_end=10
steps_per_frame=50


#goal temperature
goal=300
#creating object for each particle
class volume:  #assuming a particle 
       def __init__(self,temperature, x, y, column, row, inside):  #1 is left, or top, 2 is right or bottom
        self.x_position=x
        self.y_position=y
        self.temperature=temperature # kelvein
        self.row=row
        self.column=column
        self.inside=inside
        self.alpha=2.2*10**-5
        #incompressible so momentum just sum of forces

count=50
tolerance=1/count

filtered=[]
times=[]
#set initial conditions
for i in range(1,count):
    for j in range(1,count):
        object=volume(temperature=0, x=i/count, y=j/count, row=j, column=i, inside=True)
        filtered.append(object)

#for object in filtered:
    #position=np.abs((x_data - object.x_position)).argmin()
    #if  not(object.y_position-tolerance >= y_data_1[position] and object.y_position-tolerance <= y_data_2[position]):
     #   object.inside=False
int_sum=0
def controller(current_temperature_controller, temperature, goal, dt, int_sum):
    proportional_gain= 0   #1
    derivative_gain= 0   # 1
    integral_gain=0 #0.01 #dont make the integral term to big
    #temperature is a list of all the temperatures up to that point
    if len(temperature)>1:
        error=goal-temperature[-1]
        der_2=(temperature[-2]-temperature[-1])/dt
        if np.abs(der_2)>50:der_2*=0.1
        derivative=derivative_gain*der_2
    else:
        error=0
        derivative=0
    if np.abs(int_sum)>10**6:
         integral_gain*=0.001
    prop=error*proportional_gain
    integral=integral_gain*int_sum
    output=current_temperature_controller+prop+derivative+integral
    return 500, 10
    if output<200:
         return 200, error
    elif output>1000:
         return 1000, error
    else:
        return output, error

temperatures_list=[]
forcing=[]
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
t=0
dt= dx**2/4/filtered[0].alpha
running=True
grid = {(v.row, v.column): v for v in filtered}
while running:
    screen.fill(WHITE)
    rect_width = width/count
    
    temperary_storage={}
    for i in filtered:
            cell_east  = grid.get((i.row, i.column+1),None)
            cell_west  = grid.get((i.row, i.column-1),None)
            cell_north = grid.get((i.row-1, i.column),None)
            cell_south = grid.get((i.row+1, i.column),None)
            #boundry conditions
            #left wall is controller
            if cell_west is None:
                output1, output2=controller(i.temperature,temperatures_list, goal, dt, int_sum)
                temperary_storage[(i.row,i.column)]=output1
                
                
            elif  cell_east is None:
                temperary_storage[(i.row,i.column)]=0
            elif cell_north is None:
                temperary_storage[(i.row,i.column)]=0
            elif cell_south is None:
                temperary_storage[(i.row,i.column)]=0
            else: 
                dT = i.alpha * ((cell_east.temperature - 2*i.temperature + cell_west.temperature)/dx**2 +(cell_north.temperature - 2*i.temperature + cell_south.temperature)/dy**2)
                temperary_storage[(i.row,i.column)]=(i.temperature+dt*dT)
    output1, output2=controller(i.temperature,temperatures_list, goal, dt, int_sum)
    int_sum+=output2*dt
    forcing.append(output1)

    for i in filtered:
            i.temperature= temperary_storage[(i.row,i.column)]

    avg_temperature=0
    for i in filtered:
        avg_temperature+=i.temperature/count**2
        if i.inside:
                point_x, point_y = map_to_screen(i.x_position, i.y_position)
                if i.temperature>=0:
                    green=0
                    red=min(i.temperature*0.5, 255)
                else: 
                    green=255
                    red=0
                
                pygame.draw.rect(screen, (red, green, 0),((point_x-rect_width/2,point_y-rect_width/2,int(rect_width),int(rect_width))),0)
    temperatures_list.append(avg_temperature)
    times.append(t)
    points_1 = [map_to_screen(x_data[i], y_data_1[i]) for i in range(len(x_data))]
    pygame.draw.lines(screen, RED, False, points_1, 2)

                    # Draw curve y_data_2
    points_2 = [map_to_screen(x_data[i], y_data_2[i]) for i in range(len(x_data))]
    pygame.draw.lines(screen, RED, False, points_2, 2)
    pygame.display.flip()
    t+=dt
    # Event handling
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            t = t_end
            running=False

plt.subplot(1, 2, 1)
plt.plot(times,temperatures_list)
plt.ylabel("temperatures")
plt.subplot(1, 2, 2)
plt.plot(times,forcing)
plt.ylabel("controller")
plt.show()

df = pd.DataFrame({
    "forcing_function": forcing,
    "response_function": temperatures_list
})

# Write to CSV
df.to_csv("forcing_vs_response.csv", index=False)










