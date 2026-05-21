read_win
========

A Python-based reader for WIN-format seismic data.

Overview
--------

``read_win`` provides a Python-native approach for reading WIN-format data, designed for users who cannot rely on system-level WIN libraries.

The standard ObsPy ``read`` function does not reliably support WIN data. Known limitations (reported by S. Nakagawa and A. Kato, see references) include:

1. Incorrect handling of 0.5- and 3-byte encoded data
2. Inability to properly handle missing data and variable sampling rates
3. Performance limitations

Nakagawa and Kato provide an alternative solution using `libwinsystem.so`, which is effective but requires installation of the WIN system. This repository offers a lightweight, Pythonic alternative that avoids that dependency.


Installation
------------

It is recommended to install ``read_win`` within a conda environment.

Create a new environment:
```
conda create -n my_env -c conda-forge python=3 obspy pandas ipython
```

Activate the environment and install:
```
conda activate my_env
git clone -b optimize https://github.com/darren-tpk/read_win
cd read_win
pip install -e .
```

The package is installed in editable mode, allowing updates via ``git_pull``


Usage
-----------

See ``~/read_win/scripts/read_win_example.py``. More detailed tutorials to follow.


Conversion to miniSEED format
-----------

The ``convert_win2mseed.py`` script is provided as convenience utility to convert WIN files into miniSEED format. 

Disclaimer: Please note that it has been written for a specific field campaign. While we took extra care to include various options and possibilities, it has not been extensively tested. 

To run the script, first ensure you have activated the correct environment:

``conda activate my_env``

or

``source my_env/bin/activate``

Then, you can run the script from the command line as:

``python convert_win2mseed.py /path/to/WIN_folder/ /path/to/mseed_folder/ /path/to/channel_table NW STA --fill-value 0 -v ``

where **NW** and **STA** are the network and station names, respectively. Locations codes (e.g. '01' will be automatically extracted from the channel table. SEED channel codes (e.g. 'HDF') will be automatically determined using the data sampling_rate.

Note that the script will try to convert EVERY file inside the input folder, regardless of the format. 

The output hourly miniSEED files will be named as: ``NW.STA.01.HDF.yyy.ddd.hh``
Note that the output miniSEED files will have physical units (e.g. Pa for infrasound), using the channel table for conversion.



References
----------

Nakagawa, S. and Kato, A. — Report on ObsPy WIN-format issues:<br>
https://www.eri.u-tokyo.ac.jp/GIHOU/archive/26_031-036.pdf<br>
https://www.eri.u-tokyo.ac.jp/people/nakagawa/win/

Maeda, Y. — WIN format documentation:<br>
https://www.seis.nagoya-u.ac.jp/~maeda/ymaeda_opentools_doc/include/win/data_format.html<br>
https://www.seis.nagoya-u.ac.jp/~maeda/ymaeda_opentools_doc/win_data/index.html<br>
https://www.seis.nagoya-u.ac.jp/~maeda/ymaeda_opentools_doc/win_data/wintosac.html

Ikeda, W. and Ozaki, T. — Earlier Python conversion efforts (shared via Ichihara, M.)

Authors and Contributors
------------------------

Darren Tan<br>
Gilles Seropian

