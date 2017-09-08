#! /usr/bin/env python 

##########################################################################################
# mass_forceControl_odeint.py
#
# Script to a simulate a mass with a force input and disturbance 
#   using an ode solver
#
# NOTE: Plotting is set up for output, not viewing on screen.
#       So, it will likely be ugly on screen. The saved PDFs should look
#       better.
# 
# Created: 08/22/17
#   - Joshua Vaughan 
#   - joshua.vaughan@louisiana.edu
#   - http://www.ucs.louisiana.edu/~jev9637
#
# Modified:
#   * 
#
##########################################################################################

# Simple mass with a force input
#
#               +---> X
#               |
#            +-----+
#            |     |
#   F +--->  |  M  |<--- Fd
#            |     |
#            +-----+


import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint


def eq_of_motion(w, t, p):
    """
    Defines the differential equations for the coupled spring-mass system.

    Arguments:
        w :  vector of the state variables:
        t :  time
        p :  vector of the parameters:
    """
    x, x_dot, control_force, dist_force = w
    m, Distance, StartTime, Amax, Vmax, DistStart, F_amp = p

    # Create sysODE = (x', x_dot')
    sysODE = [x_dot,
             (U(w, t, p) - f_disturbance(w, t, p))/m,
             U(w, t, p),
             f_disturbance(w, t, p)]
    return sysODE


def f_disturbance(w, t, p):
    """
    defines the disturbance force input to the system
    """
    x, x_dot, control_force, dist_force = w
    m, Distance, StartTime, Amax, Vmax, DistStart, F_amp = p
    
    # Select one of the two inputs below
    # Be sure to comment out the one you're not using
    
    # No disturbance force 
    # f = 0
    
    # Input Option 1: 
    #    Just a step in force beginning at t=DistStart
    # f = F_amp * (t >= DistStart)
    
    # Disturbance Option 2:
    # Just a pulse - similar to running into something momentarily or a bump:
    #   A pulse in force beginning at t=DistStart and ending at t=(DistStart+0.1)
    f = F_amp * (t >= DistStart) * (t <= DistStart + 0.1)
    
    # Disturbance Option 3:
    #   viscous damper to ground
#     f = 0.05 # * (np.abs(x_dot) > 1e-2)
        
            
    return f


def U(w, t, p):
    """
    Defines the accel input to the system.
    
    Depending on the desired move distance, max accel, and max velocity, the input is either
    bang-bang or bang-coast-bang
    """
    x, x_dot, control_force, dist_force = w
    m, Distance, StartTime, Amax, Vmax, DistStart, F_amp = p
    
    
    # Bang-bang
    # These are the times for a bang-coast-bang input 
#     t1 = StartTime
#     t2 = (Vmax/Amax) + t1
#     t3 = (Distance/Vmax) + t1
#     t4 = (t2 + t3) - t1
#     
#     if t3 <= t2: # command should be bang-bang, not bang-coast-bang
#         t2 = np.sqrt(Distance/Amax)+t1
#         t3 = 2*np.sqrt(Distance/Amax)+t1
#         
#         accel = Amax*(t > t1) - 2*Amax*(t > t2) + Amax*(t > t3)
#     
#     else: # command is bang-coast-bang
#         accel = Amax*(t > t1) - Amax*(t > t2) - Amax*(t > t3) + Amax*(t > t4)
#     
#     return accel * m  * (0.01 * np.random.random() + 0.99)
# 
#     # On/off
#     if Distance - x > 1e-1:
#         f = Amax
#     elif Distance - x < -1e-1:
#         f = -Amax
#     else:
#         f = 0

    # This is a simple PD controller
    f = (kp * (Distance - x) + kd * (-x_dot)) * (t > StartTime)
    
    if f > Amax * m:
        f = Amax * m
    elif f < -Amax * m:
        f = -Amax * m

    return f
    



#---- Main script -----------------------------------------------------------------------

# Define the parameters for simluation
m = 1.0                      # mass (kg)


# ODE solver parameters
abserr = 1.0e-6
relerr = 1.0e-6
max_step = 0.01
stoptime = 5.0
numpoints = 5001

# Create the time samples for the output of the ODE solver
t = np.linspace(0.0, stoptime, numpoints)

# Initial conditions
x_init = 0.0                        # initial position
x_dot_init = 0.0                    # initial velocity


# Set up the parameters for the input function
Distance = 1.0                # Desired move distance (m)
Amax = 2.0                    # acceleration limit (m/s^2)
Vmax = 2.0                    # velocity limit (m/s)
StartTime = 0.5               # Time the y(t) input will begin
DistStart = 1.5               # Time the disturbance input will begin
F_amp = 0.5                   # Amplitude of Disturbance force (N)
kp = 10                       # Proportional Gain

wn = np.sqrt(kp/m)

zeta = np.sqrt(2)/2
kd = 2 * m * zeta * wn                        # Derivative Gain


# Pack the parameters and initial conditions into arrays 
p = [m, Distance, StartTime, Amax, Vmax, DistStart, F_amp]
x0 = [x_init, x_dot_init, 0, 0]

# Call the ODE solver.
resp = odeint(eq_of_motion, x0, t, args=(p,), atol=abserr, rtol=relerr,  hmax=max_step)
              

# Calculate the forces
control_force = np.zeros_like(t)
disturbance_force = np.zeros_like(t)
for index, time in enumerate(t):
    control_force[index] = U(resp[index,:], time, p)
    disturbance_force[index] = f_disturbance(resp[index,:], time, p)

#----- Plot the velocity response
# Set the plot size - 3x2 aspect ratio is best
fig = plt.figure(figsize=(6,4))
ax = plt.gca()
plt.subplots_adjust(bottom=0.17, left=0.17, top=0.96, right=0.96)

# Change the axis units font
plt.setp(ax.get_ymajorticklabels(),fontsize=18)
plt.setp(ax.get_xmajorticklabels(),fontsize=18)

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# Turn on the plot grid and set appropriate linestyle and color
ax.grid(True,linestyle=':', color='0.75')
ax.set_axisbelow(True)

# Define the X and Y axis labels
plt.xlabel('Time (s)',family='serif',fontsize=22,weight='bold',labelpad=5)
plt.ylabel('Velocity (m/s)',family='serif',fontsize=22,weight='bold',labelpad=10)
 
plt.plot(t, resp[:,1], linewidth=2, linestyle='-', label=r'$x$')

if F_amp > 0:
    plt.annotate(r'Force Disturbance Begins',
         xy=(DistStart,resp[int(DistStart/0.001),0]), xycoords='data',
         ha='center',
         xytext=(DistStart+1, 1.05*np.max(resp[:,0])), textcoords='data', fontsize=16,
         arrowprops=dict(arrowstyle="simple, head_width = 0.35, tail_width=0.05", connectionstyle="arc3", color="black"),color = "black")
    

# uncomment below and set limits if needed
# plt.xlim(0,5)
# plt.ylim(0,10)

# Create the legend, then fix the fontsize
# leg = plt.legend(loc='upper right', ncol = 1, fancybox=True)
# ltext  = leg.get_texts()
# plt.setp(ltext,fontsize=18)

# Adjust the page layout filling the page using the new tight_layout command
plt.tight_layout(pad=0.5)

# save the figure as a high-res pdf in the current folder
plt.savefig('mass_velocity_prop_Kp10_Z0p707_pulseDist_limited.pdf')



#----- Plot the response
# Set the plot size - 3x2 aspect ratio is best
fig = plt.figure(figsize=(6,4))
ax = plt.gca()
plt.subplots_adjust(bottom=0.17, left=0.17, top=0.96, right=0.96)

# Change the axis units font
plt.setp(ax.get_ymajorticklabels(),fontsize=18)
plt.setp(ax.get_xmajorticklabels(),fontsize=18)

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# Turn on the plot grid and set appropriate linestyle and color
ax.grid(True,linestyle=':', color='0.75')
ax.set_axisbelow(True)

# Define the X and Y axis labels
plt.xlabel('Time (s)',family='serif',fontsize=22,weight='bold',labelpad=5)
plt.ylabel('Position (m)',family='serif',fontsize=22,weight='bold',labelpad=10)
 
plt.plot(t, resp[:,0], linewidth=2, linestyle='-', label=r'$x$')

if F_amp > 0:
    plt.annotate('Force Disturbance\nBegins',
         xy=(DistStart,resp[int(DistStart/0.001),0]), xycoords='data',
         ha='center', va='top',
         xytext=(DistStart+2, .75*np.max(resp[:,0])), textcoords='data', fontsize=16,
         arrowprops=dict(arrowstyle="simple, head_width = 0.35, tail_width=0.05", connectionstyle="arc3", color="black"),color = "black")
    

# uncomment below and set limits if needed
# plt.xlim(0,5)
# plt.ylim(0,10)

# Create the legend, then fix the fontsize
# leg = plt.legend(loc='upper right', ncol = 1, fancybox=True)
# ltext  = leg.get_texts()
# plt.setp(ltext,fontsize=18)

# Adjust the page layout filling the page using the new tight_layout command
plt.tight_layout(pad=0.5)

# save the figure as a high-res pdf in the current folder
plt.savefig('mass_position_prop_Kp10_Z0p707_pulseDist_limited.pdf')


# Set the plot size - 3x2 aspect ratio is best
fig = plt.figure(figsize=(6,4))
ax = plt.gca()
plt.subplots_adjust(bottom=0.17, left=0.17, top=0.96, right=0.96)

# Change the axis units font
plt.setp(ax.get_ymajorticklabels(),fontsize=18)
plt.setp(ax.get_xmajorticklabels(),fontsize=18)

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# Turn on the plot grid and set appropriate linestyle and color
ax.grid(True,linestyle=':', color='0.75')
ax.set_axisbelow(True)

# Define the X and Y axis labels
plt.xlabel('Time (s)', fontsize=22, weight='bold', labelpad=5)
plt.ylabel('Force (N)', fontsize=22, weight='bold', labelpad=10)
 
plt.plot(t, control_force, linewidth=2, linestyle='-', label=r'Control')

if F_amp > 0:
    plt.plot(t, -disturbance_force, linewidth=2, linestyle='--', label=r'Disturbance')
    
    ## Create the legend, then fix the fontsize
    leg = plt.legend(loc='upper right', ncol = 1, fancybox=True)
    ltext  = leg.get_texts()
    plt.setp(ltext,fontsize=18)

# uncomment below and set limits if needed
# plt.xlim(0,5)
# plt.ylim(0,10)

# Adjust the page layout filling the page using the new tight_layout command
plt.tight_layout(pad=0.5)

# save the figure as a high-res pdf in the current folder
plt.savefig('mass_forces_prop_Kp10_Z0p707_pulseDist_limited.pdf')
# plt.savefig('bang_bang_command.pdf')

# show the figure
plt.show()





# 
# # If there is a non-zero force disturbance show where it began via an annotation
# if F_amp > 0:
#     annotate(r'Force Disturbance Begins',
#          xy=(DistStart,resp[-1,2]), xycoords='data',
#          ha='center',
#          xytext=(DistStart, 1.05*np.max(resp[:,0])), textcoords='data', fontsize=16,
#          arrowprops=dict(arrowstyle="simple, head_width = 0.35, tail_width=0.05", connectionstyle="arc3", color="black"),color = "black")
#     
# leg = legend(loc='upper right', fancybox=True)
# ltext  = leg.get_texts() 
# setp(ltext,family='Serif',fontsize=16)
# 
# # Adjust the page layout filling the page using the new tight_layout command
# tight_layout(pad=0.5)
# 
# # If you want to save the figure, uncomment the commands below. 
# # The figure will be saved in the same directory as your IPython notebook.
# # Save the figure as a high-res pdf in the current folder
# # savefig('MassSpringDamper_Disturbance_Resp.pdf')
# 
# 
# #----- Let's look at the forces acting on the mass
# Fsp = k * (resp[:,2] - resp[:,0])       # Spring Force (N)
# Fd =  c * (resp[:,3] - resp[:,1])       # Damping Force (N)
# F_pos = Fsp + Fd                        # Total force from position input
# 
# # Calculate the disturbance force over time by calling the disturbance force function
# F_dist = np.zeros_like(t)
# 
# for ii in range(len(t)):
#     F_dist[ii] = -f(t[ii],p)
#     
# # Now, let's plot the forces    
# # Make the figure pretty, then plot the results
# #   "pretty" parameters selected based on pdf output, not screen output
# #   Many of these setting could also be made default by the .matplotlibrc file
# fig = figure(figsize=(6,4))
# ax = gca()
# subplots_adjust(bottom=0.17,left=0.17,top=0.96,right=0.96)
# setp(ax.get_ymajorticklabels(),family='serif',fontsize=18)
# setp(ax.get_xmajorticklabels(),family='serif',fontsize=18)
# ax.spines['right'].set_color('none')
# ax.spines['top'].set_color('none')
# ax.xaxis.set_ticks_position('bottom')
# ax.yaxis.set_ticks_position('left')
# ax.grid(True,linestyle=':',color='0.75')
# ax.set_axisbelow(True)
# 
# xlabel('Time (s)',family='serif',fontsize=22,weight='bold',labelpad=5)
# ylabel('Force (N)',family='serif',fontsize=22,weight='bold',labelpad=10)
# 
# # You may need to reset these limits based on the forces in your simulation
# ymax = 1.1 * np.max([np.max(np.abs(F_pos)), np.max(np.abs(F_dist))])
# ylim(-ymax, ymax)
# 
# # plot the response
# plot(t,F_pos, linewidth=2, linestyle = '-', label=r'Spring-Damper')
# plot(t,F_dist, linewidth=2, linestyle = '--', label=r'Disturbance')
# 
# leg = legend(loc='best', fancybox=True)
# ltext  = leg.get_texts() 
# setp(ltext,family='Serif',fontsize=16)
# 
# # Adjust the page layout filling the page using the new tight_layout command
# tight_layout(pad=0.5)
# 
# # If you want to save the figure, uncomment the commands below. 
# # The figure will be saved in the same directory as your IPython notebook.
# # Save the figure as a high-res pdf in the current folder
# # savefig('MassSpringDamper_Disturbance_Forces.pdf')
# 
