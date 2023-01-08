import warnings
import os
from io import BytesIO
from zipfile import ZipFile

import seaborn as sns
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Ignore Seaborn future warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Figure layout
FIGURE_HEIGHT = 2  # inch
FIGURE_WIDTH = 8  # inch


def group_colors():
    return sns.color_palette("Dark2", 7)

def seaborn_style():
    """
    Set figures style
    """
    sns.set(style="ticks", context="paper",
            font="Arial",
            rc={"font.size": 9,
                "axes.titlesize": 9,
                "axes.labelsize": 9,
                "lines.linewidth": 1,
                "xtick.labelsize": 7,
                "ytick.labelsize": 7,
                "savefig.transparent": True,
                "xtick.major.size": 2.5,
                "ytick.major.size": 2.5,
                "xtick.minor.size": 2,
                "ytick.minor.size": 2,
                })
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42


def figpath():

    repo_dir = os.path.dirname(os.path.realpath(__file__))
    # Create figure directory
    fig_dir = os.path.join(repo_dir, 'exported_figs')
    
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)
    return fig_dir
















