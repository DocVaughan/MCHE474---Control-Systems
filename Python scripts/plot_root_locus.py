#! /usr/bin/env python

###############################################################################
# plot_root_locus.py
#
# script to plot the root locus of a system defined by a transfer function
#
# NOTE: Any plotting is set up for output, not viewing on screen.
#       So, it will likely be ugly on screen. The saved PDFs should look
#       better.
#
# Created: 11/08/17
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

# Define the system
num = [1]
den = [1, 9, 18, 0]
sys1 = control.tf(num, den)

num2 = [1, 2]
den2 = [1, 12]
sys2 = control.tf(num2,den2)

sys = sys1# * sys2

# Get the poles and zeros
poles, zeros = control.pzmap.pzmap(sys, Plot=False)

# Now, we can plot the root locus, using the control System Library's built in 
# `.root_locus()` function. We'll add a few options:
#  * `kvect=np.linspace(0,50,50001)` : Specify the range of gains to plot the locus over
#  * `Plot=False` : Don't plot the results (since we want to better control the styling)
roots, gains = control.root_locus(sys, kvect = np.linspace(0,500,100001), Plot=False)

# Set the plot size - 3x2 aspect ratio is best
fig = plt.figure(figsize=(6, 4))
ax = plt.gca()
plt.subplots_adjust(bottom=0.17, left=0.17, top=0.96, right=0.96)

# Change the axis units to serif
plt.setp(ax.get_ymajorticklabels(), family='serif', fontsize=18)
plt.setp(ax.get_xmajorticklabels(), family='serif', fontsize=18)

ax.spines['left'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_position('zero')

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('right')

# Turn on the plot grid and set appropriate linestyle and color
ax.grid(True,linestyle=':', color='0.75')
ax.set_axisbelow(True)

# Define the X and Y axis labels
# Define the X and Y axis labels
plt.xlabel('$\sigma$', family='serif', fontsize=22, weight='bold', labelpad=5)
ax.xaxis.set_label_coords(1.02, .525)


plt.ylabel('$j\omega$', family='serif', fontsize=22, weight='bold', rotation=0, labelpad=10)
ax.yaxis.set_label_coords(9/10, 1.05)

plt.plot(np.real(poles), np.imag(poles), linestyle='', 
         marker='x', markersize=10, markeredgewidth=5, 
         zorder = 10, label=r'Poles')

plt.plot(np.real(zeros), np.imag(zeros), linestyle='', 
         marker='o', markersize=10, markeredgewidth=3, markerfacecolor='none', 
         zorder = 10, label=r'Zeros')

plt.plot(roots.real, roots.imag, linewidth=2, linestyle='-', label=r'Data 1')


# plt.xticks([-5,-4,-3,-2,-1,0, 1],['','','','','','',''], bbox=dict(facecolor='black', edgecolor='None', alpha=0.65 ))
# plt.yticks([-1.0, -0.5, 0, 0.5, 1],['','', '0', '',''])
# plt.xticks(np.arange(-14, 2.1, 2))

# uncomment below and set limits if needed
plt.xlim(-9, 1)
plt.ylim(-6, 6)

# Adjust the page layout filling the page using the new tight_layout command
plt.tight_layout(pad=0.5)

# Uncomment to save the figure as a high-res pdf in the current folder
# It's saved at the original 6x4 size
plt.savefig('MidTerm2_Prob3_rootLocus.pdf')

plt.show()

# 
# # To find the root for a particular gain, we'll find the index of 
# # the gains vector that most closely matches the desired gain.
# # Then, we'll get the pol
# z0p6_roots = roots[np.argmin(np.abs(0.6-gains))]
# z2_roots = roots[np.argmin(np.abs(2-gains))]
# z4_roots = roots[np.argmin(np.abs(4-gains))]
# 
# # Set the plot size - 3x2 aspect ratio is best
# fig = plt.figure(figsize=(6, 4))
# ax = plt.gca()
# plt.subplots_adjust(bottom=0.17, left=0.17, top=0.96, right=0.96)
# 
# # Change the axis units to serif
# plt.setp(ax.get_ymajorticklabels(), family='serif', fontsize=18)
# plt.setp(ax.get_xmajorticklabels(), family='serif', fontsize=18)
# 
# ax.spines['left'].set_color('none')
# ax.spines['top'].set_color('none')
# ax.spines['bottom'].set_position('zero')
# ax.spines['right'].set_position('zero')
# 
# ax.xaxis.set_ticks_position('bottom')
# ax.yaxis.set_ticks_position('right')
# 
# # Turn on the plot grid and set appropriate linestyle and color
# ax.grid(True,linestyle=':', color='0.75')
# ax.set_axisbelow(True)
# 
# # Define the X and Y axis labels
# # Define the X and Y axis labels
# plt.xlabel('$\sigma$', family='serif', fontsize=22, weight='bold', labelpad=5)
# ax.xaxis.set_label_coords(1.02, .525)
# 
# 
# plt.ylabel('$j\omega$', family='serif', fontsize=22, weight='bold', rotation=0, labelpad=10)
# ax.yaxis.set_label_coords(4/5, 1.05)
# 
# plt.plot(np.real(poles), np.imag(poles), linestyle='', 
#          marker='x', markersize=10, markeredgewidth=5, 
#          zorder = 10, label=r'')
# 
# plt.plot(np.real(zeros), np.imag(zeros), linestyle='', 
#          marker='o', markersize=10, markeredgewidth=3, markerfacecolor='none', 
#          zorder = 10, label=r'')
# 
# plt.plot(roots.real, roots.imag, linewidth=2, linestyle='-', label=r'')
# 
# plt.plot(z0p6_roots.real, z0p6_roots.imag, linewidth=2, marker='o', markersize=10, linestyle='', label=r'$z = 0.1$')
# plt.plot(z2_roots.real, z2_roots.imag, linewidth=2, marker='*', markersize=10, linestyle='', label=r'$z = 1.0$')
# plt.plot(z4_roots.real, z4_roots.imag, linewidth=2, marker='d', markersize=10, linestyle='', label=r'$z = 5.0$')
# 
# # Create the legend, then fix the fontsize
# leg = plt.legend(loc='upper left', ncol = 1, fancybox=True)
# ltext  = leg.get_texts()
# plt.setp(ltext, family='serif', fontsize=20)
# 
# 
# # plt.xticks([-5,-4,-3,-2,-1,0, 1],['','','','','','',''], bbox=dict(facecolor='black', edgecolor='None', alpha=0.65 ))
# # plt.yticks([-1.0, -0.5, 0, 0.5, 1],['','', '0', '',''])
# 
# # uncomment below and set limits if needed
# plt.xlim(-4, 1)
# # plt.ylim(0, 10)
# 
# # Adjust the page layout filling the page using the new tight_layout command
# plt.tight_layout(pad=0.5)
# 
# # Uncomment to save the figure as a high-res pdf in the current folder
# # It's saved at the original 6x4 size
# plt.savefig('ProblemE7p13_rootLocus_withValues.pdf')
# 
# plt.show()