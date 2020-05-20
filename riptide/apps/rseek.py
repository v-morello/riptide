import logging
import argparse

import numpy as np
import pandas

from riptide import TimeSeries, ffa_search, find_peaks
from riptide.clustering import cluster1d


help_formatter = lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=16)


def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=help_formatter,
        description=(
        "FFA search a single time series and print a table of parameters of all significant peaks found."
        " Peaks found with nearly identical periods at different trial pulse widths are grouped,"
        " but no harmonic filtering is performed."
        )
    )
    parser.add_argument(
        "--Pmin", type=float, default=0.5, 
        help="Minimum trial period in seconds"
    )
    parser.add_argument(
        "--Pmax", type=float, default=10.0, 
        help="Maximum trial period in seconds"
    )
    parser.add_argument(
        "--bmin", type=int, default=240, 
        help="Minimum number of phase bins used in the search"
    )
    parser.add_argument(
        "--bmax", type=int, default=260, 
        help="Maximum number of phase bins used in the search"
    )
    parser.add_argument(
        "--smin", type=float, default=7.0, 
        help="Only report peaks above this minimum S/N"
    )
    parser.add_argument(
        "--clrad", type=float, default=0.2, 
        help=(
            "Clustering radius in units of 1/Tobs. Peaks with similar frequencies"
            " are grouped together, and only the brightest one of the group"
            " is printed"
        )
    )
    parser.add_argument(
        "-f", "--format", type=str, choices=('presto', 'sigproc'), required=True, 
        help="Input TimeSeries format"
    )
    parser.add_argument(
        "fname", type=str, 
        help="Input file name"
    )
    return parser


def run_program(args):
    """
    Run the rseek program and return a pandas DataFrame with the detected peak 
    parameters, or None if no significant peaks were found. This is used to 
    check the results in unit tests.

    Parameters
    ----------
    args : list
        List of command line arguments

    Returns
    -------
    peaks : pandas.DataFrame
        DataFrame with columns: 'period', 'freq', 'width', 'ducy', 'dm', 'snr'
    """
    logging.basicConfig(
        level='DEBUG',
        format='%(asctime)s %(filename)18s:%(lineno)-4s %(levelname)-8s %(message)s'
    )

    LOADERS = {
        'sigproc': TimeSeries.from_sigproc,
        'presto' : TimeSeries.from_presto_inf
    }

    loader = LOADERS.get(args.format, None)
    if not loader:
        raise ValueError(f"Unknown format {args.format!r}")

    # Search and find peaks
    ts = loader(args.fname)
    __, __, pgram = ffa_search(
        ts,
        period_min=args.Pmin,
        period_max=args.Pmax,
        bins_min=args.bmin,
        bins_max=args.bmax
    )
    peaks, __ = find_peaks(pgram, smin=args.smin, clrad=args.clrad)

    if not peaks:
        print(f"No peaks found above S/N = {args.smin:.2f}")
        return None

    # Cluster peaks, i.e. for each period keep only the trial width
    # that yield the highest S/N
    freqs = np.asarray([p.freq for p in peaks])
    cluster_indices = cluster1d(freqs, r=args.clrad/ts.length)
    peaks = [
        max([peaks[ii] for ii in indices], key=lambda p: p.snr)
        for indices in cluster_indices
    ]
    peaks = sorted(peaks, key=lambda p: p.snr, reverse=True)

    # DataFrame constructs from namedtuples nicely
    df = pandas.DataFrame(peaks)
    df = df.drop(columns=['iw', 'ip'])

    # Print this in a pleasing way
    # https://stackoverflow.com/questions/20937538/how-to-display-pandas-dataframe-of-floats-using-a-format-string-for-columns
    # NOTE: we have inserted a leading space to each format string on purpose
    # This makes the output table more readable
    formatters = {
        'period': '  {:.9f}'.format,
        'freq': '  {:.9f}'.format,
        'ducy': lambda x: '  {:#.2f}%'.format(100 * x),
        'dm': '  {:.2f}'.format,
        'snr': '  {:.1f}'.format,
    }
    output = df.to_string(
        columns=['period', 'freq', 'width', 'ducy', 'dm', 'snr'],
        formatters=formatters,
        index=False
    )
    print(output)
    return df


def main():
    """
    Console script entry point for 'rseek'
    """
    args = get_parser().parse_args()
    run_program(args)


if __name__ == "__main__":
    main()