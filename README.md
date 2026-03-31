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
git clone https://github.com/darren-tpk/read_win
cd read_win
pip install -e .
```

The package is installed in editable mode, allowing updates via ``git_pull``


Usage
-----------

Tutorial and examples will be added.


References
----------

Nakagawa, S. and Kato, A. — Report on ObsPy WIN-format issues:
* https://www.eri.u-tokyo.ac.jp/GIHOU/archive/26_031-036.pdf
* https://www.eri.u-tokyo.ac.jp/people/nakagawa/win/

Ikeda, W. — Early Python conversion efforts (shared via Mie Ichihara

Maeda, Y. — WIN format documentation:
* https://www.seis.nagoya-u.ac.jp/~maeda/ymaeda_opentools_doc/include/win/data_format.html
* https://www.seis.nagoya-u.ac.jp/~maeda/ymaeda_opentools_doc/win_data/index.html
* https://www.seis.nagoya-u.ac.jp/~maeda/ymaeda_opentools_doc/win_data/wintosac.html

Ikeda, W. — Early Python conversion efforts (shared via Ichihara, M.)

Authors and Contributors
------------------------

Darren Tan
Gilles Seropian

