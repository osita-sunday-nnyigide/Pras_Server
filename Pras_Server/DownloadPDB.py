#!/usr/bin/env python

__doc__ = """

This program requires python 3.6 or higher.

This module, DownloadPDB.py is used to

automatically download a PDB file if it

does not exist in user's directory

where Pras_Server program is run.

"""

__author__     = "Osita Sunday Nnyigide"

__copyright__  = "Copyright 2022, Osita Sunday Nnyigide"

__credits__    = ["Tochukwu Olunna Nnyigide", "Lee Sun-Gu", "Hyun Kyu"]

__license__    = "MIT"

__version__    = "1.2.1"

__maintainer__ = "Osita Sunday Nnyigide"

__email__      = "osita@protein-science.com"

__status__     = "Production"

__date__       = "May 11, 2022"

import os
import sys
import gzip
from urllib.request import urlretrieve
from urllib.request import urlcleanup

def downloadFile(code):

    """
    This function downloads the PDB structure file
    in .pdb or .cif format using ftp protocol.
    The downloaded file is stored in  user's
    current working directory

    Arguments
    ----------
    code: this can be a 4 letter code or
          code+format, i.e., it accepts xxxx or
          xxxx.pdb or xxxx.cif. The default
          behaviour is to download .pdb if the
          format is not specified.

    Returns
    -------
    None:   it downloads the PDB structure file

    """

    code = code.lower()

    if len(code) ==4:
        f_format = 'pdb'
    else:
        f_format = code[-3:] if code.endswith('b') else 'mmCIF'

    code = code[:4]
    f_in =   'pdb'+code+'.ent.gz' if f_format == 'pdb' else code+'.cif.gz'
    url =  "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/%s/%s/%s" % (f_format,code[1:3],f_in)

    path = os.getcwd()
    f_in = os.path.join(path, f_in)
    f_format = 'cif' if f_format == 'mmCIF' else f_format
    f_out = os.path.join(path, 'pdb'+'{}.{}'.format(code, f_format))

    if os.path.exists(f_out):
        return

    try:
        urlcleanup()
        urlretrieve(url, f_in)
    except:
        print('The specified file does not exist.'+'\n' 'Program has terminated abnormally')
        sys.exit(1)

    else:
        with gzip.open(f_in, "rb") as gz:
            with open(f_out, "wb") as out:
                out.writelines(gz)
        os.remove(f_in)

    return
