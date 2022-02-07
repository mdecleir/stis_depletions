# merge.py: This script merges the different orders in a STIS spectrum.

import glob
import numpy as np

from astropy.io import fits
from astropy.table import Table
from matplotlib import pyplot as plt
from scipy.interpolate import *


def merge_orders(filename):
    """
    Function to merge the orders of the spectrum

    Parameters
    ----------
    filename : string
        full path of the file

    Returns
    -------
    Fits file with the merged spectrum
    """

    # obtain the data
    data = fits.getdata(filename, ext=1)
    sdqflags = fits.getval(filename, ext=1, keyword="SDQFLAGS")

    # create a new wavelength grid between the minimum and maximum wavelength of the spectrum, with dispersion (step size) equal to the largest dispersion of all the orders
    # find the largest dispersion
    new_delta = 0
    for j in range(len(data)):
        delta_wave = np.diff(data["WAVELENGTH"][j])
        if np.max(delta_wave) > new_delta:
            new_delta = np.max(delta_wave)

    # find the minimum (last order, first wavelength) and maximum (first order, last wavelength) wavelength
    good_last = (data["DQ"][-1] & sdqflags) == 0
    min_wave = data["WAVELENGTH"][-1][good_last][0]
    good_first = (data["DQ"][0] & sdqflags) == 0
    max_wave = data["WAVELENGTH"][0][good_first][-1]

    # create the new wavelength grid
    new_waves = np.arange(min_wave, max_wave, new_delta, dtype="float64")

    # place all orders (len(data)) on this wavelength grid, using nearest neighbors
    # do this for the fluxes, uncertainties, data quality flags and net count rates (i.e. 4 arrays)
    # only use the good data
    new_data = np.full((len(data), 4, len(new_waves)), np.nan)
    for j in range(len(data)):
        good = (data["DQ"][j] & sdqflags) == 0
        func = interpolate.interp1d(
            data["WAVELENGTH"][j][good],
            [
                data["FLUX"][j][good],
                data["ERROR"][j][good],
                data["DQ"][j][good],
                data["NET"][j][good],
            ],
            kind="nearest",
            bounds_error=False,
        )
        new_data[j] = func(new_waves)
    # save the new data as an Astropy Table
    new_data_table = Table(new_data, names=["FLUX", "ERROR", "DQ", "NET"])

    # calculate the sensivity
    new_data_table["SENS"] = new_data_table["FLUX"] / new_data_table["NET"]

    # co-add orders with weights
    # use 1/sensivity as the weight
    weights = 1 / new_data_table["SENS"]
    summed_weights = np.nansum(weights, axis=0)
    # if all orders are nan, the summed weight should be nan (0 by default from nansum)
    summed_weights[summed_weights == 0] = np.nan
    # weighted mean flux = sum(fi*wi)/sum(wi)
    coadded_fluxes = (
        np.nansum(new_data_table["FLUX"] * weights, axis=0) / summed_weights
    )
    # standard deviation of weighted mean = sqrt[sum(si**2*wi**2)]/sum(wi)
    coadded_uncs = (
        np.sqrt(np.nansum(new_data_table["ERROR"] ** 2 * weights ** 2, axis=0))
        / summed_weights
    )

    # save the merged spectrum
    coadded_table = Table(
        [new_waves, coadded_fluxes, coadded_uncs], names=["WAVELENGTH", "FLUX", "ERROR"]
    )
    outname = filename.replace("x1f", "merged")
    coadded_table.write(outname, format="fits", overwrite=True)


if __name__ == "__main__":
    datapath = "/Users/mdecleir/Documents/Depletions/HST_data/"
    stars = ["HD203938"]

    # merge the orders of the blaze function corrected spectra
    for star in stars:
        file_list = glob.glob(datapath + star + "/*_x1f.fits")
        for file in file_list:
            merge_orders(file)
