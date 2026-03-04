import numpy as np
import pygame as pg


Force=0
cart_mass=5
ball_mass=1
g=9.80665
Length=3
#initial conditions
theta_double_dot=0
theta_dot=0
theta=0.1
F=0
a_x=0
dt=0.01
t=10
v_x=0
s_x=0
#EOM cart
# Force-sin(theta)*T=cart_mass*cart_ax

#EOM ball
#T*sin(theta)=ball_mass*ball_ax
#T*cos(theta)-ball_mass*ball_ay

#ball_ax=cart_ax-L*theta_dx2*cos(theta)+L*theta_dx**2*sin(theta)
#ball_ay=-L*theta_dx2*sin(theta)-L*theta_dx**2*cos(theta)

#governing differential equations
#-g*sin(theta)=a_x_cart*cos(theta)-L*theta_double_dot
#F+m_ball*L*theta_double_dot*cos(theta)-mass_ball*L*theta_dot**2*sin(theta)=(m_ball+m_cart)a_x




def equations_of_motion():
    global v_x, s_x, theta_dot, theta, a_x, theta_double_dot, F

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    total_mass = cart_mass + ball_mass
    temp = (F + ball_mass * Length * theta_dot**2 * sin_theta) / total_mass

    denominator = Length * (4/3 - (ball_mass * cos_theta**2) / total_mass)
    if abs(denominator) < 1e-6:
        denominator = 1e-6

    theta_double_dot = (g * sin_theta - cos_theta * temp) / denominator
    a_x = temp - (ball_mass * Length * theta_double_dot * cos_theta) / total_mass

    # Integrate
    v_x += a_x * dt
    s_x += v_x * dt

    theta_dot += theta_double_dot * dt
    theta += theta_dot * dt

    # Normalize theta
    theta = ((theta + np.pi) % (2 * np.pi)) - np.pi



Proportional=100
Integral=10
Derivative=10
sum=0

def controller(theta, P, I , D, initial_theta):
    Error=theta
    Force_required=0
    #proportional gain
    global F, sum
    Force_required+=Error*P
    #Integral controller
    sum+=Error*dt
    Force_required+=sum*I
    #derivative controller
    change=(theta-initial_theta)/dt
    Force_required+=change*D
    F=Force_required



pg.init()
WIDTH, HEIGHT = 1280, 800
screen = pg.display.set_mode((WIDTH, HEIGHT))
clock = pg.time.Clock()
start_ticks = pg.time.get_ticks()*dt  # milliseconds

running = True
while running:
    for event in pg.event.get():
        if event.type == pg.QUIT:
            running = False
    initial_theta=theta
    equations_of_motion()
    controller(theta, Proportional, Integral, Derivative, initial_theta)


    screen.fill((0, 0, 0))

    base_x = WIDTH // 2 + s_x

# Pendulum
    ball_x = base_x + np.sin(theta) * Length * 100
    ball_y = HEIGHT - 200 - np.cos(theta) * Length * 100

    pg.draw.circle(screen, (255, 255, 255), (ball_x, ball_y), 30)

# Cart
    square_rect = pg.Rect(base_x - 30, HEIGHT - 200 - 30, 60, 60)
    pg.draw.rect(screen, (255, 0, 0), square_rect)

    pg.display.flip()
    clock.tick(60)

pg.quit()
