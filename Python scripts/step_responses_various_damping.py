#! /usr/bin/env python

###############################################################################
# step_responses_various_damping.py
#
# Script to plot the step response of a 2nd-order underdamped system for various damping ratios
#
# NOTE: Any plotting is set up for output, not viewing on screen.
#       So, it will likely be ugly on screen. The saved PDFs should look
#       better.
#
# Created: 09/26/17
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
import control

damping = [0, 0.1, 0.2, 0.4, 0.7, 1.0, 2.0]
wn = 2 * np.pi     # natural frequency (rad)

# The time range over which to simulate the response
t = np.linspace(0, 3, 301)

# Define a matrix to hold all the responses. Defining it as all zeros and filling it 
# later will be significantly faster for large problems and good practice otherwise
resp = np.zeros((len(t), len(damping)))

# loop through the values of damping ratio in the list damping, simluating the system
# and saving the response as column in resp
for index, zeta in enumerate(damping):

    # Define the system for this damping ratio
    num = wn**2
    den = [1, 2 * zeta * wn, wn**2]
    sys = control.tf(num, den)
    
    resp[:,index], t_out = control.step(sys, t)
  
# Now, plot the various responses  
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
plt.xlabel('Normalized Time $(\omega_n t)$', fontsize=22, weight='bold', labelpad=5)
plt.ylabel('$y(t)$', fontsize=22, weight='bold', labelpad=10)
 
plt.plot(wn * t_out, resp[:,0], linewidth=2, linestyle='-', label=r'$\zeta = 0$')
plt.plot(wn * t_out, resp[:,1], linewidth=2, linestyle='--', label=r'$\zeta = 0.1$')
plt.plot(wn * t_out, resp[:,2], linewidth=2, linestyle='-.', label=r'$\zeta = 0.2$')
plt.plot(wn * t_out, resp[:,3], linewidth=2, linestyle=':', label=r'$\zeta = 0.4$')

# uncomment below and set limits if needed
# plt.xlim(0,5)
plt.ylim(0, 3)

# Create the legend, then fix the fontsize
leg = plt.legend(loc='upper right', ncol = 2, fancybox=True)
ltext  = leg.get_texts()
plt.setp(ltext,fontsize=18)

# Adjust the page layout filling the page using the new tight_layout command
plt.tight_layout(pad=0.5)

# save the figure as a high-res pdf in the current folder
# plt.savefig('SecondOrder_StepResponses_LightDamping.pdf')


# Now, plot the various responses  
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
plt.xlabel('Normalized Time $(\omega_n t)$', fontsize=22, weight='bold', labelpad=5)
plt.ylabel('$y(t)$', fontsize=22, weight='bold', labelpad=10)

plt.plot(wn * t_out, resp[:,0], linewidth=2, linestyle='-', label=r'$\zeta = 0.0$')
plt.plot(wn * t_out, resp[:,4], linewidth=2, linestyle='--', label=r'$\zeta = 0.7$')
plt.plot(wn * t_out, resp[:,5], linewidth=2, linestyle='-.', label=r'$\zeta = 1.0$')
plt.plot(wn * t_out, resp[:,6], linewidth=2, linestyle=':', label=r'$\zeta = 2.0$')

# uncomment below and set limits if needed
# plt.xlim(0,5)
plt.ylim(0, 3)

# Create the legend, then fix the fontsize
leg = plt.legend(loc='upper right', ncol = 2, fancybox=True)
ltext  = leg.get_texts()
plt.setp(ltext,fontsize=18)

# Adjust the page layout filling the page using the new tight_layout command
plt.tight_layout(pad=0.5)

# save the figure as a high-res pdf in the current folder
# plt.savefig('SecondOrder_StepResponses_HighDamping.pdf')

# show the figure
plt.show()