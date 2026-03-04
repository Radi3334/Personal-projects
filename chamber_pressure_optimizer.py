import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
import math

#static, quasi one d, Ideal gas
#constants
gamma=1.4
R=287
Nozzle_length=1 # don't change this, since it is normalized, so 1 is the end of the nozzle

def area_function(x_position):
    #x=0 at throat
    A=0.1+0.1*x_position**2/Nozzle_length
    return A
def Mach_number_at_location(area, area_star):
    #only holds for fully isentropic flow
    ratio_sq=(area/area_star)**2
    
    def ratio(M):
        ratio=(1/M**2)*((2/(gamma+1))*(1+((gamma-1)/2)*M**2))**((gamma+1)/(gamma-1))-ratio_sq
        return ratio
    if ratio_sq != 1:
        supersonic=brentq(ratio, 1.001, 5)
        subsonic=brentq(ratio, 0.001, 0.9999)
    else:
        supersonic=1
        subsonic=1
    return subsonic, supersonic
def finding_station_values():
    #ratio is always p0/pb
    area_exit=area_function(Nozzle_length)
    area_throat=area_function(0)
    low, high=Mach_number_at_location(area_exit, area_throat)
    #p_e6
    #A_throat=A_star
    
    ratio_6=(1+(gamma-1)/2*high**2)**(gamma/(gamma-1))

    #p_e5
    #A_throat=A_star
    exit_pressure=(1+(2*gamma)/(gamma+1)*(high**2-1))**(-1)
    ratio_5=(1+(gamma-1)/2*high**2)**(gamma/(gamma-1))*exit_pressure
    
    #p_e3
    #A_throat=A_star
    ratio_3=(1+(gamma-1)/2*low**2)**(gamma/(gamma-1))

    return ratio_6, ratio_5, ratio_3
def region_finder(chamber_pressure, back_pressure,t):
    ratio=chamber_pressure/back_pressure
    
    ratio_6, ratio_5, ratio_3=finding_station_values()
    if ratio==1:
        number=0
    elif ratio>ratio_6:
        print(f"overexpanded at time ={t}")
        return 1
    elif ratio==ratio_6:
        print(f"ideally expanded at time={t}")
        return 2
    elif ratio< ratio_6 and ratio>ratio_5:
        print(f"underexpanded at time ={t}")
        return 3
    elif ratio<= ratio_5 and ratio>ratio_3:
        print(f"shockwaves present at time ={t}")
        return 4
    elif ratio==ratio_3:
        print(f"only supersonic in throat at time ={t}") 
        return 5
    else:
        print(f"fullysubsonic at time ={t}")
        return 6
def shock_location_finder( back_pressure, chamber_pressure):
    area_exit=area_function(Nozzle_length)
    area_throat=area_function(0)    
    ratio=chamber_pressure/back_pressure
    ratio_area=area_throat/area_exit
    M_e=np.sqrt((1/(gamma-1))*(np.sqrt(1+2*(gamma-1)*(2/(gamma+1))**((gamma+1)/(gamma-1))*(ratio_area*ratio)**2)-1))
    total_pressure_after= back_pressure*((1+(gamma-1)*0.5*M_e**2)**(gamma/(gamma-1)))
    ratio_total_pressure=total_pressure_after/chamber_pressure

    def mach_finder(M):
        term1 = ((gamma + 1) * M**2) / (2 + (gamma - 1) * M**2)
        term2 = (gamma + 1) / (2 * gamma * M**2 - (gamma - 1))
        ratio_from_M=term1**(gamma / (gamma - 1)) * term2**(1 / (gamma - 1))
        return ratio_from_M-ratio_total_pressure

    mach_solution=brentq(mach_finder, 1.0001, 5)

    ratio_area2=np.sqrt((1/mach_solution**2)*((2/(gamma+1))*(1+((gamma-1)/2)*mach_solution**2))**((gamma+1)/(gamma-1)))
    area_at_location=ratio_area2*area_throat
    def location_finder(x):
        return area_function(x)-area_at_location
    guess_location=brentq(location_finder, 0, Nozzle_length)
    return guess_location, total_pressure_after
def pressure_finder(mach_number, total_pressure):
    pressure=total_pressure/((1+(gamma-1)/2*mach_number**2)**(gamma/(gamma-1)))
    return pressure
def temperature_finder(mach_number, chamber_temperature):
    Temerature= chamber_temperature*(1+(gamma-1)/2*mach_number**2)**(-1)
    return Temerature
def mach_number_at_location_subsonic(back_pressure, total_pressure, total_temp, location):
    area_exit=area_function(Nozzle_length)
    m_exit= np.sqrt(2/(gamma-1)*((total_pressure/back_pressure)**((gamma-1)/gamma)-1))
    exit_temp=total_temp/(1+(gamma-1)/2*m_exit**2)
    mass_flow2=m_exit*np.sqrt(gamma/(R*exit_temp))*back_pressure*area_exit
    
    def mass_flow(M):
        pressure=pressure_finder(M, total_pressure)
        temperature=temperature_finder(M, total_temp)
        mass_flow=M*np.sqrt(gamma/(R*temperature))*pressure*area_function(location)
        return mass_flow-mass_flow2

    mach_number_solution=brentq(mass_flow, 0,0.99)
    return mach_number_solution
def main(chamber_pressure, back_pressure, chamber_temp,t):
    class layer:
        def __init__(self, x_location, pressure, Machnumber, Temperature, total_pressure, density):
            self.x_location=x_location
            self.pressure=pressure
            self.Machnumber=Machnumber
            self.Temperature=Temperature
            self.total_pressure=total_pressure,
            self.density=density
    layers=[]
    layer_count=100
    number=region_finder(chamber_pressure, back_pressure, t)
    if number==1 or number==2 or number==3:   
        for i in range(layer_count):
            x=i/layer_count
            area=area_function(x)
            mach1, mach2=Mach_number_at_location(area, area_function(0))
            p=pressure_finder(mach2, chamber_pressure)
            T=temperature_finder(mach2, chamber_temp)
            rho=p/(R*T)
            layers.append(layer(x_location=x, pressure=p, Machnumber=mach2, Temperature=T, total_pressure=chamber_pressure, density=rho))
    elif number==4:
        shock_location, t_p_after=shock_location_finder(back_pressure, chamber_pressure)
        for i in range(layer_count):
            x=i/layer_count
            area4=area_function(x)
            if x<shock_location:   
                mach1, mach2=Mach_number_at_location(area4, area_function(0))
                p=pressure_finder(mach2, chamber_pressure)
                T=temperature_finder(mach2, chamber_temp)
                rho=p/(R*T)
                layers.append(layer(x_location=x, pressure=p, Machnumber=mach2, Temperature=T, total_pressure=chamber_pressure, density=rho))
            else:
                mach_after=mach_number_at_location_subsonic(back_pressure,t_p_after,chamber_temp,x  )
                p=pressure_finder(mach_after, t_p_after)
                T=temperature_finder(mach_after, chamber_temp)
                rho=p/(R*T)
                layers.append(layer(x_location=x, pressure=p, Machnumber=mach_after, Temperature=T, total_pressure=t_p_after, density=rho))
    elif number==5:
        for i in range(layer_count):
            x=i/layer_count
            area=area_function(x)
            mach1, mach2=Mach_number_at_location(area, area_function(0))
            p=pressure_finder(mach1, chamber_pressure)
            T=temperature_finder(mach1, chamber_temp)
            rho=p/(R*T)
            layers.append(layer(x_location=x, pressure=p, Machnumber=mach1, Temperature=T, total_pressure=chamber_pressure, density=rho))
    elif number==6:
        for i in range(layer_count):
            x=i/layer_count
            area=area_function(x)
            mach1=mach_number_at_location_subsonic(back_pressure,chamber_pressure,chamber_temp,x  )
            p=pressure_finder(mach1, chamber_pressure)
            T=temperature_finder(mach1, chamber_temp)
            rho=p/(R*T)
            layers.append(layer(x_location=x, pressure=p, Machnumber=mach1, Temperature=T, total_pressure=chamber_pressure, density=rho))
    elif number==7:
        print("no flow")
    return layers
def plot(layers):
    x_locations = [layer.x_location for layer in layers]
    mach_numbers = [layer.Machnumber for layer in layers]
    pressures = [layer.pressure for layer in layers]
    temperatures = [layer.Temperature for layer in layers]
    densities = [layer.density for layer in layers]
    areas = [area_function(x) for x in x_locations]

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Converging-Diverging Nozzle Flow Analysis\nChamber Pressure: {chamber_pressure/1000:.1f} kPa, Back Pressure: {back_pressure/1000:.1f} kPa', fontsize=14, fontweight='bold')

    # Plot 1: Mach number
    axes[0, 0].plot(x_locations, mach_numbers, 'g-', linewidth=2)
    axes[0, 0].axhline(y=1, color='r', linestyle='--', alpha=0.5, label='Mach = 1')
    axes[0, 0].set_xlabel('Nozzle Position (normalized)')
    axes[0, 0].set_ylabel('Mach Number')
    axes[0, 0].set_title('Mach Number Distribution')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend()

    # Plot 2: Pressure distribution
    axes[0, 1].plot(x_locations, np.array(pressures)/1000, 'r-', linewidth=2)
    axes[0, 1].axhline(y=back_pressure/1000, color='k', linestyle='--', alpha=0.5, label='Back Pressure')
    axes[0, 1].set_xlabel('Nozzle Position (normalized)')
    axes[0, 1].set_ylabel('Pressure (kPa)')
    axes[0, 1].set_title('Pressure Distribution')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend()

    # Plot 3: Temperature distribution
    axes[1, 0].plot(x_locations, temperatures, 'orange', linewidth=2)
    axes[1, 0].set_xlabel('Nozzle Position (normalized)')
    axes[1, 0].set_ylabel('Temperature (K)')
    axes[1, 0].set_title('Temperature Distribution')
    axes[1, 0].grid(True, alpha=0.3)

    # Plot 4: Density distribution
    axes[1, 1].plot(x_locations, np.array(densities), 'purple', linewidth=2)
    axes[1, 1].set_xlabel('Nozzle Position (normalized)')
    axes[1, 1].set_ylabel('Density (kg/m³)')
    axes[1, 1].set_title('Density Distribution')
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


#set altitude an pressure profiles
def ISA(altitude):
    Temperature11000= 216.65
    Pressure11000= 22626.31034841721
    Pressure20000= 5469.405328263262
    Temperature20000= 216.65
    Temperature32000=228.65
    Pressure32000= 866.8540511142719
    Temperature47000=270.65
    Pressure47000= 110.71536817451754
    Pressure51000=66.80586596547876
    Temperature71000= 214.65
    Pressure71000= 3.946494350415755
    if 0 <= altitude <= 11000:
        Temperature= 288.15 - (0.0065 * altitude)
        Pressure= 101325*math.pow((Temperature)/288.15, 5.256767623)
  
    elif 11000 < altitude <= 20000:
        Temperature = Temperature11000
        Pressure = Pressure11000*math.exp(-(9.81/(287*Temperature))*(altitude-11000))
        
    elif 20000 < altitude <= 32000:
        Temperature= Temperature11000 + (0.001 * (altitude-20000))
        Pressure= Pressure20000*math.pow((Temperature)/Temperature20000, -34.1695122)

    elif 32000 < altitude <= 47000:
        Temperature= Temperature32000 + (0.0028 * (altitude-32000))
        Pressure= Pressure32000*math.pow(((Temperature)/Temperature32000), -12.20339721)
        
    elif 47000 < altitude <= 51000:
        Temperature= Temperature47000
        Pressure = Pressure47000*math.exp(-(9.81/(287*Temperature))*(altitude-47000))


    elif 51000 < altitude <= 71000:
        Temperature= Temperature47000 - (0.0028 * (altitude-51000))
        Pressure= Pressure51000*math.pow(((Temperature)/Temperature47000), 12.20339721)


    elif 71000 < altitude <= 86000:
        Temperature= Temperature71000 - (0.002 * (altitude-71000))
        Pressure= Pressure71000*math.pow(((Temperature)/Temperature71000), 17.0847561)
    else:
        Pressure=100
    return Pressure
def altitude_profile(t):
    if t>0:
        return max(0,1000*math.log(t)+1000*t)
    else: return 0



#Simulation
chamber_pressure=120000 #pa
chamber_temp=300 #k
t_end=500
dt=0.1

ideal_expanded_pressures=[]
altitudes=[]
time=[]
ideal_pressure_ratio,__,__=finding_station_values()
for i in np.arange(0,t_end, dt):
    
    altitude=altitude_profile(i)
    if altitude>86000:
        break
    back_pressure=ISA(altitude)
    time.append(i)
    ideal_expanded_pressures.append(ideal_pressure_ratio*back_pressure/1000)
    altitudes.append(altitude/1000)
    if chamber_pressure<back_pressure:
        print("wrong pressures")
        break
    list_of_layers=main(chamber_pressure, back_pressure, chamber_temp,np.round(i,2))
    print(f"pressure ratio at t={np.round(i,2)} at an altitude of h={np.round(altitude,2)}:",np.round(chamber_pressure/back_pressure,2))
    #you can change how many graphs you want in what time interval
    if i%10==0:
        plot(list_of_layers)

fig, axs = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)

# Time vs pressure
axs[0, 0].plot(time, ideal_expanded_pressures, linewidth=2)
axs[0, 0].set_xlabel('Time [s]')
axs[0, 0].set_ylabel('Ideal chamber pressure [kPa]')
axs[0, 0].set_title('Pressure vs Time')
axs[0, 0].grid(True, alpha=0.3)

# Altitude vs pressure
axs[0, 1].plot(altitudes, ideal_expanded_pressures, linewidth=2)
axs[0, 1].set_xlabel('Altitude [km]')
axs[0, 1].set_ylabel('Ideal chamber pressure [kPa]')
axs[0, 1].set_title('Pressure vs Altitude')
axs[0, 1].grid(True, alpha=0.3)

# Optional: turn off unused subplots
axs[1, 0].axis('off')
axs[1, 1].axis('off')

fig.suptitle('Ideal Chamber Pressure During Simulation', fontsize=14)

plt.show()

#next steps: import an actual area distribution of a rocket nozzle, and then program in ISA