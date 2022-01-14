# extract.py: script to extract 1D spectra from 2D STIS data (from the flt fits files)
# based on the python notebook created by Sean Lockwood:
# 2d_combination_echelle.ipynb

import glob
import os
import stistools


def extract(filename):
    outname = filename.replace(".fits", "_x1d.fits")
    # remove the file if it already exists
    if os.access(outname, os.F_OK):
        os.remove(outname)

    # extract the 1D spectrum
    res = stistools.x1d.x1d(filename, output=outname, ctecorr="omit")


if __name__ == "__main__":
    datapath = "/Users/mdecleir/Documents/Depletions/HST_data/"
    stars = [
        "HD014434",
        "HD037903",
        "HD038087",
        "HD073882",
        "HD093250",
        "HD149404",
        "HD152249",
        "HD152723",
        "HD168076",
        "HD190603",
        "HD192639",
        "HD197770",
        "HD200775",
        "HD203938",
        "HD206267",
        "HD216898",
        "HD239729",
    ]

    for star in stars:
        file_list = glob.glob(datapath + star + "/*comb.fits")
        if len(file_list) > 0:
            for filename in file_list:
                extract(filename)
        else:
            print(star + " has no combined exposures.")
