# cont_fit.py: Script to fit the continuum flux around a line

import glob
import os
import numpy as np

from astropy.io import ascii, fits
from astropy.table import Table
from matplotlib import pyplot as plt

from plotting.plot_spectra import velo


def norm(datapath, star, linewaves, degree):
    """
    Function to normalize the spectrum around every line

    Parameters
    ----------
    datapath : string
        Path to the data files

    star : string
        Name of the star

    linewaves : np
        Wavelengths of the lines

    degree:
        Degree of the polynomial

    Returns
    -------
    Plot with the continuum
    """
    # check if the file with windows to fit the continuum exists
    window_file = datapath + star + "/" + star + "_cont_win.dat"
    if not os.path.isfile(window_file):
        print("Please, create windows file for star " + star, window_file)
    else:
        # read the windows file
        windows = ascii.read(window_file)

        # create the Table to store the continuum coefficients
        ncols = degree + 1
        colnames = ["C" + str(i + 1) for i in range(ncols)]
        colnames.insert(0, "line_wave")
        colnames.insert(1, "deg")
        coeff_table = Table(names=colnames)

        # obtain the data files with the merged spectra for this star
        file_list = glob.glob(datapath + star + "/*merged.fits")
        for filename in file_list:
            # obtain the data
            data = fits.getdata(filename, ext=1)
            waves = data["WAVELENGTH"]
            fluxes = data["FLUX"]

            # for all absorption lines
            for lwave in linewaves:
                # check if the line is in the wavelength range of the spectrum
                if np.logical_and(lwave > waves[0], lwave < waves[-1]):
                    # convert the wavelengths to velocities
                    velos = velo(waves, lwave)

                    # check if the line is included in the windows file
                    if not lwave in windows["line_wave"]:
                        print(
                            "Please, add continuum windows for line with wavelength "
                            + str(lwave)
                            + " \u212B to the windows file."
                        )
                        # open a plot to determine appropriate continuum windows
                        fig, ax = plt.subplots()
                        plt.plot(velos, fluxes, "k")
                        plt.xlim([-550, 550])
                        plt.ylim(0, 1.5 * np.nanmedian(fluxes))
                        plt.text(0.84, 0.9, lwave, transform=ax.transAxes)
                        plt.show()
                    else:
                        # obtain the velocities in the fit windows
                        rows = windows["line_wave"] == lwave
                        vmins = windows["vmin"][rows]
                        vmaxs = windows["vmax"][rows]
                        mask = np.full(len(velos), False)
                        for vmin, vmax in zip(vmins, vmaxs):
                            mask += (velos > vmin) & (velos < vmax)
                        # exclude NaNs from the fitting
                        mask *= ~np.isnan(fluxes)

                        # fit the continuum with a Legendre polynomial
                        coefs = np.polynomial.legendre.legfit(
                            velos[mask], fluxes[mask], degree
                        )
                        cont = np.polynomial.legendre.legval(velos[mask], coefs)

                        # add polynomial coefficients to the table
                        coeff_table.add_row(np.insert(coefs, 0, [lwave, degree]))

                        # normalize the spectrum around the line
                        norm = fluxes[mask] / cont

                        # make a plot
                        fix, ax = plt.subplots()
                        # plot the data
                        plt.plot(velos, fluxes, c="k", alpha=0.7, label="data")
                        # plot the continuum fit
                        plt.plot(
                            velos[mask],
                            cont,
                            c="tab:orange",
                            label="cont. fit (" + str(degree) + ")",
                        )
                        # plot the data points that were used in the continuum fit
                        plt.plot(
                            velos[mask],
                            fluxes[mask],
                            "r.",
                            markersize=2.5,
                            alpha=0.6,
                            label="fit points",
                        )
                        # finalize and save the plot
                        plt.xlim(-510, 510)
                        plt.ylim([0, 1.5 * np.median(fluxes[mask])])
                        plt.xlabel("velocity (km/s)", fontsize=fs)
                        plt.ylabel(
                            r"flux ($ergs\ cm^{-2}\ s^{-1}\ \AA^{-1}$)", fontsize=fs
                        )
                        plt.legend(loc=3)
                        plt.savefig(
                            datapath
                            + star
                            + "/"
                            + star
                            + "_cont_fit_"
                            + str(int(lwave))
                            + "_"
                            + str(degree)
                            + ".pdf",
                            bbox_inches="tight",
                        )
        # finalize and save the table
        coeff_table["deg"] = coeff_table["deg"].astype(int)
        coeff_table.write(
            datapath + star + "/" + star + "_cont_coeff.dat",
            format="ascii",
            overwrite=True,
        )


if __name__ == "__main__":
    datapath = "/Users/mdecleir/Documents/Depletions/HST_data/"
    stars = ["HD203938"]

    # plotting properties for uniform plots
    fs = 18

    # define the wavelenghts of the lines of interest
    line_waves = np.array([1239.9253, 1240.3947, 1355.5977, 2249.8768, 2260.7805])

    # normalize the spectrum around every line
    degree = 9

    for star in stars:
        norm(datapath, star, line_waves, degree)
