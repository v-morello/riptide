Unless stated otherwise, all test time series:  
* Contain 16 samples: the integers from 0 to 15 inclusive
* Have a 64 us sampling time

The data are stored as `float32` for the PRESTO time series. For the `SIGPROC` data, the data type is stated in the file name.
The file `fake_sigproc_uint8_nosignedkey.tim` is a copy of `fake_sigproc_uint8.tim` where the `signed` header key has been deleted, and we test that riptide correctly refuses to process it since it cannot interpret the data with certainty.