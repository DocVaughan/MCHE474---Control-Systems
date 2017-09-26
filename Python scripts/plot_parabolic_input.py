#! /usr/bin/env python

###############################################################################
# plot_parabola_input.py
#
# simple script to plot a ramp input, written for MCHE474
#
# NOTE: Any plotting is set up for output, not viewing on screen.
#       So, it will likely be ugly on screen. The saved PDFs should look
#       better.
#
# Created: 09/25/17
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

STEP_TIME = 0.0     # time at which the step should occur

# Define the time vector
time = np.linspace(-0.25, 3, 3001)

# Define the ramp command, the boolean operator is what created the step
parabola = time**2 * (time > STEP_TIME)

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
plt.ylabel('$r(t)$', fontsize=22, weight='bold', labelpad=10)
 
plt.plot(time, parabola, linewidth=2, linestyle='-', label=r'Data 1')

# uncomment below and set limits if needed
plt.xlim(-0.25,2)
plt.ylim(0,4)

plt.yticks([0, 2, 4],[0, '', ''])

# Create the legend, then fix the fontsize
# leg = plt.legend(loc='upper right', ncol = 1, fancybox=True)
# ltext  = leg.get_texts()
# plt.setp(ltext,fontsize=18)

# Adjust the page layout filling the page using the new tight_layout command
plt.tight_layout(pad=0.5)

# save the figure as a high-res pdf in the current folder
plt.savefig('parabola_input.pdf')

# show the figure
plt.show()