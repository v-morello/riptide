import os
import numpy as np
from riptide import TimeSeries


INF_TEMPLATE = """
 Data file name without suffix          =  {basename:s}
 Telescope used                         =  Parkes
 Instrument used                        =  Multibeam
 Object being observed                  =  Pulsar
 J2000 Right Ascension (hh:mm:ss.ssss)  =  00:00:01.0000
 J2000 Declination     (dd:mm:ss.ssss)  =  -00:00:01.0000
 Data observed by                       =  Kenji Oba
 Epoch of observation (MJD)             =  59000.000000
 Barycentered?           (1=yes, 0=no)  =  1
 Number of bins in the time series      =  {nsamp:d}
 Width of each time series bin (sec)    =  {tsamp:.12e}
 Any breaks in the data? (1=yes, 0=no)  =  0
 Type of observation (EM band)          =  Radio
 Beam diameter (arcsec)                 =  981
 Dispersion measure (cm-3 pc)           =  {dm:.12f}
 Central freq of low channel (Mhz)      =  1182.1953125
 Total bandwidth (Mhz)                  =  400
 Number of channels                     =  1024
 Channel bandwidth (Mhz)                =  0.390625
 Data analyzed by                       =  Space Sheriff Gavan
 Any additional notes:
    Input filterbank samples have 2 bits.
"""

def generate_data_presto(outdir, basename, tobs=128.0, tsamp=256e-6, period=1.0, dm=0.0, amplitude=20.0, ducy=0.05):
    """
    Generate some time series data with a fake signal, and save it in PRESTO
    inf/dat format in the specified output directory.

    Parameters
    ----------
    outdir : str
        Path to the output directory
    basename : str
        Base file name (not path) under which the .inf and .dat files 
        will be saved.
    **kwargs: self-explanatory
    """
    ### IMPORTANT: seed the RNG to get reproducible results ###
    np.random.seed(0)

    ts = TimeSeries.generate(
        tobs, tsamp, period, 
        amplitude=amplitude, ducy=ducy, stdnoise=1.0
    )
    inf_text = INF_TEMPLATE.format(basename=basename, nsamp=ts.nsamp, tsamp=tsamp, dm=dm)

    inf_path = os.path.join(outdir, f"{basename}.inf")
    dat_path = os.path.join(outdir, f"{basename}.dat")
    with open(inf_path, 'w') as fobj:
        fobj.write(inf_text)
    ts.data.tofile(dat_path)
