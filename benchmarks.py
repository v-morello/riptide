import time

import numpy as np

from riptide.libffa import LibFFA




def get_width_trials(b, wmax=0.25, wtsp=1.5):
    """ Generate the list of width trials """
    widths = []
    w = 1
    wmax = max(1, wmax * b)
    while w <= wmax:
        widths.append(w)
        w = max(w + 1, int(w * wtsp))
    return np.asarray(widths)
    

def runtime_avg(m, b, snr=True, runs=20, threads=1):
    data = np.zeros(shape=(m, b), dtype=np.float32)
    out = np.zeros_like(data)
    widths = get_width_trials(b)
    snrs = np.zeros(shape=(m, widths.size), dtype=np.float32)
    varnoise = 1.0

    runtimes = []

    for irun in range(runs):
        start = time.time()
        LibFFA.py_transform(data, m, b, out)
        if snr:
            LibFFA.get_snr_2d(out, m, b, widths, widths.size, varnoise, threads, snrs)
        end = time.time()
        runtimes.append(end - start)

    return sum(runtimes) / runs


if __name__ == '__main__':

    b_trials = [256, 384, 512, 768, 1024]
    m_trials = 2 ** np.arange(7, 14)

    print('%8s %8s %12s %16s' % ('m', 'b', 'tobs (min)', 'runtime (ms)'))
    for b in b_trials:
        for m in m_trials:
            t = runtime_avg(m, b, runs=10, snr=True, threads=1)
            print('%8d %8d %12.1f %16.2f' % (m, b, m / 60.0, t * 1000.0))
