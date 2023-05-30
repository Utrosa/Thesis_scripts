# code for figures describing the KS (discrete) and PDF (continuous)
# distributions

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pdf_tonality import pdf_tonality as pdft

# set to None to show plot, or provide a path
dir_out = None
dir_out = os.path.join(os.getcwd(),'figures')

# define the ratings (1-7) collected by Krumhansl & Kessler (1982), as listed
# in the MIDI Toolbox code
# see: https://github.com/miditoolbox/1.1/blob/master/miditoolbox/refstat.m
weights_maj = [6.35,2.23,3.48,2.33,4.38,4.09,2.52,5.19,2.39,3.66,2.29,2.88]
intervals_maj = np.array([0,2,4,5,7,9,11])
weights_min = [6.33,2.68,3.52,5.38,2.6,3.53,2.54,4.75,3.98,2.69,3.34,3.17]
intervals_min = np.array([0,2,3,5,7,8,10])

# collect the evaluations of the probability density functions
PDF_maj, _, _ = pdft.build_diatonic_PDF('major')
PDF_min, _, _ = pdft.build_diatonic_PDF('minor')

# set up subplots
fig = plt.figure(figsize=(12, 6), dpi=100)
ax1 = fig.add_subplot(221) # K-S major
ax2 = fig.add_subplot(223) # PDF major
ax3 = fig.add_subplot(222) # K-S minor
ax4 = fig.add_subplot(224) # PDF minor
fig.tight_layout(pad=5.0)

# settings for color/fill
legend_elements = [
    Patch(facecolor='#404040', edgecolor='#404040',label='Diatonic'),
    Patch(facecolor='#BABABA', edgecolor='#BABABA',label='Non-diatonic'),
]

# subplot 1 - upper left - KK major
nondiatonic = ax1.bar(
    x=np.arange(0,1200,100),
    height=weights_maj,
    width=50,
    color='#BABABA',
)
for i in intervals_maj:
    ax1.bar(x=i*100, height=weights_maj[i], width=50, color='#404040')
ax1.set_xlim([-50, 1250])
ax1.set_ylim([1, 7])
ax1.set_yticks(range(1,8))
ax1.set_xticks(np.arange(0,1200,100))
ax1.set_xticklabels(range(0,12))
ax1.set_ylabel('Rating')
ax1.set_xlabel('Semitones from tonic')
ax1.set_title('Krumhansl-Kessler (1982) probe-tone ratings (Major)')
ax1.legend(handles=legend_elements, loc="upper right")

# subplot 2 - lower left - PDF major
ax2.plot(PDF_maj, color='#404040')
ax2.set_xlim([-50,1250])
ax2.axvline(x=0, c='k', linestyle='--')
ax2.axvline(x=1200, c='k', linestyle='--')
ax2.set_yticks([])
ax2.set_ylabel('Normalized probability')
ax2.set_xlabel('Cents from tonic')
ax2.set_title('Diatonic probability evaluated at each cent (Major)')

# subplot 3 - upper right - KK minor
nondiatonic = ax3.bar(
    x=np.arange(0,1200,100),
    height=weights_min,
    width=50,
    color='#BABABA',
)
for i in intervals_min:
    ax3.bar(x=i*100, height=weights_min[i], width=50, color='#404040')
ax3.set_xlim([-50, 1250])
ax3.set_ylim([1, 7])
ax3.set_yticks(range(1,8))
ax3.set_xticks(np.arange(0,1200,100))
ax3.set_xticklabels(range(0,12))
ax3.set_ylabel('Rating')
ax3.set_xlabel('Semitones from tonic')
ax3.set_title('Krumhansl-Kessler (1982) probe-tone ratings (Minor)')
ax3.legend(handles=legend_elements, loc="upper right")

# subplot 4 - lower right - PDF minor
ax4.plot(PDF_min, color='#404040')
ax4.set_xlim([-50,1250])
ax4.axvline(x=0, c='k', linestyle='--')
ax4.axvline(x=1200, c='k', linestyle='--')
ax4.set_yticks([])
ax4.set_ylabel('Normalized probability')
ax4.set_xlabel('Cents from tonic')
ax4.set_title('Diatonic probability evaluated at each cent (Minor)')

# export or show
if dir_out is not None:
    filepath = os.path.join(dir_out, 'KK_v_PDF.png')
    plt.savefig(filepath, dpi=100)
    plt.close()
else:
    plt.show()