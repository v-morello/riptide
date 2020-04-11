import numpy as np

def generate_width_trials(nbins, ducy_max=0.20, wtsp=1.5):
    widths = []
    w = 1
    wmax = int(max(1, ducy_max * nbins))
    while w <= wmax:
        widths.append(w)
        w = int(max(w + 1, wtsp * w))
    return np.asarray(widths)
