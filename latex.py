

def get_figsize(columnwidth, wf=0.5, hf=(5.**0.5-1.0)/2.0, ):
    """Parameters:
      - wf [float]:  width fraction in columnwidth units
      - hf [float]:  height fraction in columnwidth units.
                     Set by default to golden ratio.
      - columnwidth [float]: width of the column in latex. Get this from LaTeX
                             using \showthe\columnwidth
    Returns:  [fig_width,fig_height]: that should be given to matplotlib
    """
    fig_width_pt = columnwidth*wf
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*hf      # height in inches
    return [fig_width, fig_height]


# obtained using \showthe\columnwidth:
COLWIDTH_BEAMER = 307.28987     # [pt]
COLWIDTH_SCRARTCL = 418.25555   # [pt]

print("beamer", get_figsize(COLWIDTH_BEAMER))
print("scrartcl", get_figsize(COLWIDTH_SCRARTCL))

