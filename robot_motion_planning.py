import numpy as np
import pygame

robot_count=10
dt=0.1

goal_x=0.5
goal_y=0.2
max_speed=0.0001

class robot:
    def __init__(self, x, y, vy, vx, thrust, throttle, mass, radius):
        self.x=x
        self.y=y
        self.vy=vy
        self.vx=vx
        self.thrust=thrust  #how fast the robot moves
        self.throttle=throttle #percentage
        self.radius=radius
        self.mass=mass
    def move_horizontal(self,throttle): #if direction is positive moves to positive x
        a_x=throttle/100*self.thrust/self.mass
        self.vx+=a_x*dt
        if np.abs(self.vx)>max_speed:
            self.vx=np.sign(self.vx)* max_speed
        self.x+=self.vx*dt
    def move_vertical(self, throttle): #if direction is positive moves to positive x
        a_y=throttle/100*self.thrust/self.mass
        self.vy+=a_y*dt
        if self.vy>max_speed:
            self.vy=np.sign(self.vy)* max_speed
        self.y+=self.vy*dt

    def distance_to_goal(self, goal_x, goal_y):
        distance=np.sqrt((goal_x-self.x)**2+(goal_y-self.y)**2)
        return distance
    
    def direction_finder(self, goal_x, goal_y):
        x_error= self.x-goal_x
        y_error= self.y-goal_y
        if x_error<0:
            angle= np.atan2(-x_error, -y_error)
        elif x_error>0:angle= np.atan2(x_error, y_error)
        else: angle=np.sign(-y_error)*np.pi/2
        return angle
    
    def throttle_control(self, direction):
        throttle_x=100*np.cos(direction)
        throttle_y=100*np.sin(direction)
        return throttle_x, throttle_y


class wall:
    def __init__(self, x, y, orientation):
        self.x=x
        self.y=y
        self.orientation=orientation #true if horizontal, false if vertical


pygame.init()
font = pygame.font.SysFont(None, 24)
width, height = 800, 800
screen = pygame.display.set_mode((width, height))
pygame.display.set_caption("Robot Motion")
WHITE = (255, 255, 255)
RED = (255, 0, 0)
def map_to_screen(x, y, x_min=0, x_max=1, y_min=0, y_max=1):
    px = int((x - x_min) / (x_max - x_min) * width)
    py = int(height - (y - y_min) / (y_max - y_min) * height)  # y=0 at bottom
    return px, py


#creates robots
robots=[]
for i in range(robot_count):
    robot_int=robot(x=i/robot_count,y=0.5, vy=0,vx=0, thrust=0.00001, throttle=0, radius=10, mass=10)
    robots.append(robot_int)


t=0
running=True
while running:
    screen.fill(WHITE)
    for i in robots:
        throttle_x, throttle_y=i.throttle_control(i.direction_finder(goal_x, goal_y))
        i.move_horizontal(throttle_x)
        i.move_vertical(throttle_y)

    for i in robots:
        point_x, point_y=map_to_screen(i.x,i.y)
        pygame.draw.circle(screen, (255,0,0), (point_x,point_y), i.radius)
        text = font.render(str(np.round(np.sqrt(height**2+width**2)*i.distance_to_goal(goal_x,goal_y),0)), True, (0, 0, 0))
        text_rect = text.get_rect(center=(point_x, point_y - 15))
        screen.blit(text, text_rect)

    point_x_g, point_y_g=map_to_screen(goal_x,goal_y)
    pygame.draw.circle(screen, (0,255,0), (point_x_g,point_y_g), 10)

    pygame.display.flip()
    t+=dt
    # Event handling
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running=False