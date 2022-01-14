# combine_2d.py: script to combine STIS exposures in 2D (from the flt fits files)
# based on the python notebook created by Sean Lockwood:
# 2d_combination_echelle.ipynb

import glob
import os
import shutil
import datetime
import warnings
import numpy as np

from astropy.io import fits


def get_shifts(filename):
    # collect the wavelength shifts
    with fits.open(filename) as f:
        a1shifts = []
        a2shifts = []
        for i in range(1, len(f) // 3 + 1):
            a1shifts.append(f[("SCI", i)].header["SHIFTA1"])  # X-axis / dispersion
            a2shifts.append(
                f[("SCI", i)].header["SHIFTA2"]
            )  # Y-axis / cross-dispersion

    # calculate the range in the shifts, and the average shifts
    delta_a1 = np.ptp(a1shifts)
    delta_a2 = np.ptp(a2shifts)
    ave_a1 = np.mean(a1shifts)
    ave_a2 = np.mean(a2shifts)
    print(delta_a1, delta_a2)

    return delta_a1, delta_a2, ave_a1, ave_a2


def coadd(filename, ave_a1, ave_a2):
    # create arrays with zeros
    summed = np.zeros((1024, 1024), dtype=float)
    summed_err_sq = np.zeros((1024, 1024), dtype=float)
    summed_dq = np.zeros((1024, 1024), dtype=int)
    total_exptime = 0.0

    # sum the data, sum the squared errors, and bitwise-or combine the quality flags
    with fits.open(filename) as f:
        for i in range(1, len(f) // 3 + 1):
            summed += f[("SCI", i)].data
            summed_err_sq += f[("ERR", i)].data ** 2
            summed_dq |= f[("DQ", i)].data
            total_exptime += f[("SCI", i)].header["EXPTIME"]

    # calculate the error on the summed data
    summed_err = np.sqrt(summed_err_sq)

    # copy the input file to another file
    outname = filename.replace("flt", "comb")
    shutil.copy(filename, outname)

    # update the file with the new summed data
    with fits.open(outname, "update") as f:
        # update the data
        f[("SCI", 1)].data = summed
        f[("ERR", 1)].data = summed_err
        f[("DQ", 1)].data = summed_dq

        # update the header
        f[0].header["FILENAME"] = os.path.basename(outname)
        f[0].header["ROOTNAME"] = "combined"
        f[0].header["TEXPTIME"] = total_exptime
        f[("SCI", 1)].header["EXPTIME"] = total_exptime

        # add the average wavelength shifts
        f[("SCI", 1)].header["SHIFTA1"] = ave_a1
        f[("SCI", 1)].header["SHIFTA2"] = ave_a2
        f[0].header["WAVECORR"] = "OMIT"  # use the hard-coded SHIFTAn value above

        # add some HISTORY description of what we've done
        f[0].header.add_history("")
        f[0].header.add_history(
            f"Written {datetime.datetime.now().isoformat(' ').rsplit('.',1)[0]}:"
        )
        f[0].header.add_history(f"Modified to contain the sum of {i} exposures of")
        f[0].header.add_history(filename)
        f[0].header.add_history("Not all header keywords here will be correct!")

        # Remove other (SCI, ERR, DQ) extension sets
        while len(f) > 4:
            del f[-1]


def combine_2d(filename):
    # get the range in wavelength calibration shifts between different exposures
    delta_a1, delta_a2, ave_a1, ave_a2 = get_shifts(filename)

    # if the shifts are not too different, co-add the different exposures
    if delta_a1 < 0.5 and delta_a2 < 0.5:
        coadd(filename, ave_a1, ave_a2)
    else:
        warnings.warn(
            "The difference in wavelength shift is too big for these exposures to be co-added."
        )


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
        file_list = glob.glob(datapath + star + "/*flt.fits")
        for filename in file_list:
            with fits.open(filename) as f:
                # if the file has multiple exposures
                if len(f) > 4:
                    combine_2d(filename)
                else:
                    print(filename + " only contains 1 exposure")
