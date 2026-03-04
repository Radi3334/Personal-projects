#when supersonic flow enters an engine, you want highest total pressures for subsonic flow
import numpy as np
import pygame as pg
from scipy.optimize import brentq
import matplotlib.pyplot as plt

#so there is a freestream mach number entering a duct with two walls that create shocks due to an angle theta to the flow
#essentially a diffuser, that will decrease to a smaller super sonic flow, and then a normal shock that will bring it to a final total pressure
#two shock locations, when they meet and at the wall

#alpha is free stream angle of attack, positive downwards deflection, so what would be pitch up, abs(alpha)< theta always

def running(free_stream_mach,alpha):
    gamma=1.4
    R=287
    free_stream_pressure=10**6
    free_stream_temperature=280
    def pressure_finder(mach_number, total_pressure):
        pressure=total_pressure/((1+(gamma-1)/2*mach_number**2)**(gamma/(gamma-1)))
        return pressure
    def temperature_finder(mach_number, chamber_temperature):
        Temperature= chamber_temperature*(1+(gamma-1)/2*mach_number**2)**(-1)
        return Temperature
    def total_pressure_finder(mach_number, pressure):
        tot_pressure=pressure*((1+((gamma-1)/2)*mach_number**2)**(gamma/(gamma-1)))
        return tot_pressure
    def tot_temperature_finder(mach_number, temperature):
        Temperature= temperature*(1+(gamma-1)/2*mach_number**2)
        return Temperature
    def mach_beta(theta, M):
        #first have to find theta max, so that you know beta max, then can check between mu and beta max
        theta_max, beta_max=theta_max_finder(M)
        if theta > theta_max+0.001:
            print(f"theta={theta:.2f} > theta_max={theta_max:.2f}, cannot solve beta")
            return 90 #
        if abs(theta) < 1e-7:
            return np.rad2deg(np.arcsin(1/M))
        theta=theta*np.pi/180
        mu=np.arcsin(1/M)
        lhs=np.tan(theta)
        def beta_f(beta):
            rhs=(2/np.tan(beta))*(((M**2)*(np.sin(beta)**2)-1)/(M**2*(gamma+np.cos(2*beta))+2))
            return rhs-lhs
        try:
            beta_weak=brentq(beta_f, mu-0.001, np.deg2rad(beta_max)+0.001)
            return beta_weak*180/np.pi
        except ValueError:
            return beta_max   
    def theta_max_finder(M):
        a=(-1)*2*gamma*M**4
        b=M**4*(gamma+1)-4*M**2
        c=M**2*(gamma+1)+2
        sin_beta_sq=(-b-np.sqrt(b**2-4*a*c))/(2*a)
        sin_beta_sq = np.clip(sin_beta_sq, 0, 1)  # Clamp to valid range for arcsin
        beta=np.arcsin(np.sqrt(sin_beta_sq))
        theta=np.arctan((2/np.tan(beta))*(((M**2)*(np.sin(beta)**2)-1)/(M**2*(gamma+np.cos(2*beta))+2)))
        theta=np.rad2deg(theta)
        beta=np.rad2deg(beta)
        #this is the max beta and max theta present for a given mach number
        return theta, beta
    def shock_changes(M, beta, theta, pressure): #M is machnumber, before normilastion, beta is in degrees
        M_normal=M*np.sin(np.deg2rad(beta))
        M_2_n=np.sqrt((1+(gamma-1)/2*M_normal**2)/(gamma*M_normal**2-(gamma-1)/2))

        m2=M_2_n/np.sin(np.deg2rad(beta-theta))

        p2=pressure*(1+(2*gamma)/(gamma+1)*(M_normal**2-1))
        return m2, p2
    def total_pressure_drop(M, p_tot):
        pzero2overp1=((((gamma+1)**2*M**2)/(4*gamma*M**2-2*gamma+2))**(gamma/(gamma-1)))*((1-gamma+2*gamma*M**2)/(gamma+1))
        p1=pressure_finder(M, p_tot)
        pzero2=pzero2overp1*p1
        return(pzero2)


    #region one is free stream, region 2 is after top shock, region 3 after bottom shock, region 4 is top of slip line, region 5 is below slip line
    #beta is always the beta of the region the shock goes into

    #flow is adiabatic so total temperature is constant

    def main(theta_top, theta_bottom):
        if np.abs(alpha)<theta_top and np.abs(alpha)<theta_bottom:
            #free stream properties
            total_temperature=tot_temperature_finder(free_stream_mach, free_stream_temperature)
            total_pressure1=total_pressure_finder(free_stream_mach, free_stream_pressure)

            #region two properties
            beta2=mach_beta(theta_top+alpha,free_stream_mach)
            mach2, pressure2=shock_changes(free_stream_mach, beta2, theta_top+alpha, free_stream_pressure)
            total_pressure2=total_pressure_finder(mach2, pressure2)

            #region three properties
            beta3=mach_beta(theta_bottom-alpha,free_stream_mach)
            mach3, pressure3=shock_changes(free_stream_mach, beta3, theta_bottom-alpha, free_stream_pressure)
            total_pressure3=total_pressure_finder(mach3, pressure3)

            #assuming that region 4,5 is behind the inlet, so flow has to be subsonic and is the flow that enters the engine
            #region 4,5 we now have to find the slip line, slip line is angle phi wrt the horizontal plane
            def slip_angle_finder(slip_angle):
                #region 4
                if (mach2<1)or mach3<1:
                    print("subsonic") 
                    return None
                beta4=mach_beta(theta_top+slip_angle,mach2)
                __, pressure4=shock_changes(mach2, beta4, theta_top+slip_angle, pressure2)
                if beta4 is None: return 1e10
                #region 5
                beta5=mach_beta(theta_bottom-slip_angle,mach3)
                __, pressure5=shock_changes(mach3, beta5, theta_bottom-slip_angle, pressure3)
                if beta5 is None: return -1e10
                p_diff = pressure4 - pressure5
                return p_diff
            

            theta_max1,__=theta_max_finder(mach2)
            theta_max2,__=theta_max_finder(mach3)
            slip_max = min(theta_max1 - theta_top, theta_bottom)  
            slip_min = -min(theta_max2 - theta_bottom, theta_top)
            
            if np.sign(slip_angle_finder(slip_min)*slip_angle_finder(slip_max))==-1:

                slip_angle=brentq(slip_angle_finder, slip_min, slip_max) #in degrees
                
                #region 4
                beta4=mach_beta(theta_top+slip_angle,mach2)
                mach4, pressure4=shock_changes(mach2, beta4, theta_top+slip_angle, pressure2)
                total_pressure4=total_pressure_finder(mach4, pressure4)
                temp4=temperature_finder(mach4, total_temperature)
                    #region 5
                beta5=mach_beta(theta_bottom-slip_angle,mach3)
                mach5, pressure5=shock_changes(mach3, beta5, theta_bottom-slip_angle, pressure3)
                total_pressure5=total_pressure_finder(mach5, pressure5)
                temp5=temperature_finder(mach5, total_temperature)
                if mach4>1 and mach5>1:
                    print(f"Slip-line angle               : {slip_angle:.4f} deg")
                    print(f"Region 4 total-pressure ratio : {total_pressure4/total_pressure1:.5f}")
                    print(f"Region 5 total-pressure ratio : {total_pressure5/total_pressure1:.5f}")
                    print(f"Region 4 Mach number          : {mach4:.4f} [-]")
                    print(f"Region 5 Mach number          : {mach5:.4f} [-]")
                #ok now last normal shock -> essentially an oblique shock with beta =90
                #assuming the flow will mix after the slip line,  af each region -> not entirely accurate but a close approximation
                #so we will average the properties which is fine for a preliminary assumption
                    pressuretot6=(total_pressure4+total_pressure5)/2
                    #region 7 is what the engine will see after the normal shock
                    pressuretot7=total_pressure_drop((mach4+mach5)/2, pressuretot6)
                    total_pressure_loss=pressuretot7/total_pressure1
                    mach7,__= shock_changes((mach4+mach5)/2, 90, 0, pressure5)

                    print("=== Terminal shock / region 7 results ===")
                    print(f"Region 7 Mach number             : {mach7:.4f} [-]")
                    print(f"Region 7 total-pressure ratio    : {total_pressure_loss:.5f} [-]")

                    return total_pressure_loss
                else: 
                    print("inlet shocks are too big and flow is now subsonic, so no normal shock will occur")
                    return 
            else:
                return 0 #since do not want the flow to be subsonic after the oblique shocks but only after the normal one



    #now we will optimize the inlet theta's 


    theta_max, beta_max= theta_max_finder(free_stream_mach)
    total_pressures={}
    lower_bound=alpha+1
    upper_bound=np.floor(theta_max)-1
    angle_range=np.arange(lower_bound, upper_bound,1)
    length=len(angle_range)
    grid=np.zeros((length,length))

    for i in np.arange(lower_bound, upper_bound, 1):
        for j in np.arange(alpha+1, np.floor(theta_max)-5, 1):
            print(f"theta_top:{i}, theta_bottom:{j}")
            try:
                pressuredrop=main(i,j)
                if pressuredrop is None:
                    pressuredrop=0
                total_pressures[i,j]=pressuredrop
                grid[int(i-lower_bound), int(j-lower_bound)] = pressuredrop
            except:
                grid[int(i-lower_bound), int(j-lower_bound)] = 0
    max_location = max(total_pressures, key=total_pressures.get)
    max_location=np.round(max_location)



    plt.imshow(grid, 
            extent=[angle_range[0], angle_range[-1], angle_range[0], angle_range[-1]], 
            origin='lower', 
            cmap='viridis')
    plt.colorbar(label='Total Pressure loss')
    plt.xlabel('Theta Bottom (degrees)')
    plt.ylabel('Theta Top (degrees)')
    plt.title(f'Total pressure drop (M={free_stream_mach}), (Aoa={alpha}), max at theta_top={max_location[0]}, theta_bottom={max_location[-1]}')
    plt.show()
    return max_location

Mach_input=5
positions=[]
aoas=[]
for i in range(0,10,1):
    positions.append(running(Mach_input,i))
    aoas.append(i)

positions = np.array([p for p in positions if p is not None], dtype=float)
plt.scatter(positions[:, 0], positions[:, 1], s=30, c='C0')
plt.plot(positions[:, 0], positions[:, 1], '-o', color='C0')
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('theta_top'); plt.ylabel('theta_bottom')
plt.title('2D positions')
plt.show()