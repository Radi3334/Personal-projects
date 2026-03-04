import matplotlib.pyplot as plt
import numpy as np
import pygame

#graphing tube
x_data = np.linspace(0, 1, 100)
y_data_1= 0.4*np.exp(-(x_data - 0.5)**2)
y_data_2=1-y_data_1

gamma=1.4
t_end=10

#creating object for each particle
class volume:  #assuming a particle 
       def __init__(self, rho, u1, u2, v1, v2, p, x, y, column, row, inside):  #1 is left, or top, 2 is right or bottom
        self.rho = rho
        u = u2 - u1
        v = v2 - v1
        self.E = p/(gamma-1) + 0.5 * rho * (u*u + v*v)
        self.x_position=x
        self.y_position=y
        self.p=p
        self.row=row
        self.column=column
        self.inside=inside
        self.u1=u1
        self.u2=u2
        self.v1=v1
        self.v2=v2
        self.p_old=0
        self.rho_old=0
        self.E_old=0
        #incompressible so momentum just sum of forces

count=50
tolerance=1/count

filtered=[]
#set initial conditions
for i in range(1,count):
    for j in range(1,count):
        if i>5*count/6:
            object=volume(rho=100, u1=0,u2=0,  v1=0,v2=0, p=0.5*10**6, x=i/count, y=j/count, row=j, column=i, inside=True)
            filtered.append(object)
        else:
            object=volume(rho=100, u1=0,u2=0,  v1=0,v2=0, p=i*40*10**6, x=i/count, y=j/count, row=j, column=i, inside=True)
            filtered.append(object)

for object in filtered:
    position=np.abs((x_data - object.x_position)).argmin()
    if  not(object.y_position-tolerance >= y_data_1[position] and object.y_position-tolerance <= y_data_2[position]):
        object.inside=False



faces_dict_horizontal={}
faces_dict_vertical={} #horizontal faces are always given by the coordinate of the cell below it, vertical faces by the cell left of it
#value of item in dictionary is the value of the flow speed moving

#still have to create faces for the bottom
for object in filtered:
    faces_dict_horizontal[(object.row,object.column)]=0
    faces_dict_vertical[(object.row,object.column)]=0


#make dictionaries to find the faces, then you can use the coordinates of the cells to find the face you want

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
# Main game loop
dx,dy=1/count, 1/count
t=0
dt= 0.1#min(dx/(np.abs(i.u1-i.u2)+a),dy/(np.abs(i.v1-i.v2)+a))
grid = {(v.row, v.column): v for v in filtered}
running=True


#you have to make it that it computes the pressures first, and then it recomputes them a few times
#change the code to faces
while running:
    screen.fill(WHITE)
    rect_width = width/count
    
    for i in filtered:
        cell_east  = grid.get((i.row, i.column+1),None)
        cell_west  = grid.get((i.row, i.column-1),None)
        cell_north = grid.get((i.row-1, i.column),None)
        cell_south = grid.get((i.row+1, i.column),None)
        
        #boundry conditions
        if cell_west is None and i.inside:
            i.u1 =0
            pressure_gradient_east=cell_east.p-i.p
            pressure_gradient_north=cell_north.p-i.p
            pressure_gradient_south=cell_south.p-i.p

            au2=-pressure_gradient_east/(i.rho*dx)
            faces_dict_vertical(i.row,i.column)+=au2
            av1=-pressure_gradient_north/(i.rho*dx)
            faces_dict_horizontal(i.row,i.column)+=av1

        elif  cell_east is None and i.inside:
            i.u2 = 100
            pressure_gradient_west=cell_west.p-i.p
            pressure_gradient_north=cell_north.p-i.p
            pressure_gradient_south=cell_south.p-i.p
            au1=-pressure_gradient_west/(i.rho*dx)
            i.u1+=au1*dt
            av1=-pressure_gradient_north/(i.rho*dx)
            i.v1+=av1*dt
            av2=-pressure_gradient_south/(i.rho*dx)
            i.v2+=av2*dt
        elif cell_north is None and i.inside:
            i.v1 = 0
            pressure_gradient_east=cell_east.p-i.p
            pressure_gradient_west=cell_west.p-i.p
            pressure_gradient_south=cell_south.p-i.p
            au1=-pressure_gradient_west/(i.rho*dx)
            i.u1+=au1*dt
            au2=-pressure_gradient_east/(i.rho*dx)
            i.u2+=au2*dt
            av2=-pressure_gradient_south/(i.rho*dx)
            i.v2+=av2*dt
        elif cell_south is None and i.inside:
            i.v2 = 0
            pressure_gradient_east=cell_east.p-i.p
            pressure_gradient_west=cell_west.p-i.p
            pressure_gradient_north=cell_north.p-i.p

            au1=-pressure_gradient_west/(i.rho*dx)
            i.u1+=au1*dt
            au2=-pressure_gradient_east/(i.rho*dx)
            i.u2+=au2*dt
            av1=-pressure_gradient_north/(i.rho*dx)
            i.v1+=av1*dt
        else: 
            if i.inside:
                pressure_gradient_east=cell_east.p-i.p
                pressure_gradient_west=cell_west.p-i.p
                pressure_gradient_north=cell_north.p-i.p
                pressure_gradient_south=cell_south.p-i.p
                #positive flow is outwards
                au1=-pressure_gradient_west/(i.rho*dx)
                i.u1+=au1*dt
                au2=-pressure_gradient_east/(i.rho*dx)
                i.u2+=au2*dt
                av1=-pressure_gradient_north/(i.rho*dy)
                i.v1+=av1*dt
                av2=-pressure_gradient_south/(i.rho*dy)
                i.v2+=av2*dt
                
                #i.u1=np.sign(pressure_gradient_west)*np.sqrt(abs(pressure_gradient_east/(i.rho/2+cell_east.rho/2)))
                #i.u2=-1*np.sign(pressure_gradient_east)*np.sqrt(abs(pressure_gradient_west/(i.rho/2+cell_west.rho/2)))
                #i.v1=np.sign(pressure_gradient_north)*np.sqrt(abs(pressure_gradient_north/(i.rho/2+cell_north.rho/2)))
                #i.v2=-1*np.sign(pressure_gradient_south)*np.sqrt(abs(pressure_gradient_south/(i.rho/2+cell_south.rho/2)))
        
        i.p_old=i.p
        i.E_old=i.E
        i.rho_old=i.rho
    storage_density={}
    storage_pressure={}
    for i in filtered:
        if i.inside:
            cell_east  = grid.get((i.row, i.column+1),None)
            cell_west  = grid.get((i.row, i.column-1),None)
            cell_north = grid.get((i.row-1, i.column),None)
            cell_south = grid.get((i.row+1, i.column),None)
            #continuity equation: flux out is positive
            if cell_west is None and i.inside:
                pressure=(cell_east.p+cell_north.p+cell_south.p)/3-dx*i.rho*(i.u2-i.u1+i.v1-i.v2)/4/dt
            elif  cell_east is None and i.inside:
                pressure=(cell_north.p+cell_south.p+cell_west.p)/3-dx*i.rho*(i.u2-i.u1+i.v1-i.v2)/4/dt
            elif cell_north is None and i.inside:
                pressure=(cell_east.p+cell_south.p+cell_west.p)/3-dx*i.rho*(i.u2-i.u1+i.v1-i.v2)/4/dt
            elif cell_south is None and i.inside:
                pressure=(cell_east.p+cell_north.p+cell_west.p)/3-dx*i.rho*(i.u2-i.u1+i.v1-i.v2)/4/dt
            else: 
                if i.inside:
                    pressure=(cell_east.p+cell_north.p+cell_south.p+cell_west.p)/4-dx*i.rho*(i.u2-i.u1+i.v1-i.v2)/4/dt
            storage_pressure[(i.row,i.column)]=pressure
            if (cell_east is not None):
                uE = i.u2
                rhoE = i.rho_old if uE > 0 else cell_east.rho_old
                flux_E = rhoE * uE
            else:
                flux_E=0
            if (cell_west is not None):
                uW = i.u1
                rhoW = cell_west.rho_old if uW > 0 else i.rho_old
                flux_W = rhoW * uW
            else:
                flux_W=0
            if (cell_north is not None):
                vN = i.v1
                rhoN = cell_north.rho_old if vN > 0 else i.rho_old
                flux_N = rhoN * vN
            else:
                flux_N=0
            if (cell_south is not None):
                vS = i.v2
                rhoS = i.rho_old if vS > 0 else cell_south.rho_old
                flux_S = rhoS * vS
            else:
                flux_S=0
            density = max((i.rho_old - dt/dx * (flux_E - flux_W) - dt/dy * (flux_N - flux_S)),0.01)
            storage_density[(i.row,i.column)]=density
    for i in filtered:
        if i.inside:
            i.rho=storage_density[(i.row,i.column)]
            i.p=storage_pressure[(i.row,i.column)]
    for i in filtered:
        if i.inside:
            point_x, point_y = map_to_screen(i.x_position, i.y_position)
            if i.u2>0:
                green=max(0,min(255,i.u2*100))
                red=0
            else: 
                green=0
                red=max(0,min(np.abs(i.u2*100),255))
            pygame.draw.rect(screen, (red, green, 0),((point_x-rect_width/2,point_y-rect_width/2,int(rect_width),int(rect_width))),0)
    
    # Draw curve y_data_1
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













