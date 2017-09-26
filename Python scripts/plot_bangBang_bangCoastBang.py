#! /usr/bin/env python

###############################################################################
# plot_bangBang_bangCoastBang.py
#
# simple script to plot bang-bang or bang-coast-bang commands
#
# NOTE: Any plotting is set up for output, not viewing on screen.
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
# TODO:
#   * 
###############################################################################

import numpy as np
import matplotlib.pyplot as plt


def U(t, p):
    """
    Defines the accel input to the system.
    
    Depending on the desired move distance, max accel, and max velocity, the input is either
    bang-bang or bang-coast-bang
    """
    Distance, StartTime, Amax, Vmax = p
    
    # Bang-bang
    # These are the times for a bang-coast-bang input 
    t1 = StartTime
    t2 = (Vmax/Amax) + t1
    t3 = (Distance/Vmax) + t1
    t4 = (t2 + t3) - t1
    
    if t3 <= t2: # command should be bang-bang, not bang-coast-bang
        t2 = np.sqrt(Distance/Amax)+t1
        t3 = 2*np.sqrt(Distance/Amax)+t1
        
        accel = Amax*(t > t1) - 2*Amax*(t > t2) + Amax*(t > t3)
    
    else: # command is bang-coast-bang
        accel = Amax*(t > t1) - Amax*(t > t2) - Amax*(t > t3) + Amax*(t > t4)
    
    return accel

# Create the time samples for the output of the ODE solver
t = np.linspace(0.0, 5.0, 5001)

# Set up the parameters for the input function
Distance = 2.0                # Desired move distance (m)
Amax = 2.0                   # acceleration limit (m/s^2)
Vmax = 2.0                    # velocity limit (m/s)
StartTime = 0.5               # Time the y(t) input will begin

# Pack the parameters and initial conditions into arrays 
p = [Distance, StartTime, Amax, Vmax]

# Calculate the forces
control_force = U(t, p)


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
# plt.plot(t, disurbance_force, linewidth=2, linestyle='--', label=r'Disturbance')

plt.yticks([-Amax, -Amax/2, 0, Amax/2, Amax], ['$-u_{max}$', '', '$0$', '', '$u_{max}$'])

# uncomment below and set limits if needed
# plt.xlim(0,5)
# plt.ylim(0,10)

## Create the legend, then fix the fontsize
# leg = plt.legend(loc='upper right', ncol = 2, fancybox=True)
# ltext  = leg.get_texts()
# plt.setp(ltext,fontsize=18)

# Adjust the page layout filling the page using the new tight_layout command
plt.tight_layout(pad=0.5)

# save the figure as a high-res pdf in the current folder
plt.savefig('bang_bang_command.pdf')

# show the figure
plt.show()
