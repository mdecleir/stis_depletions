# norm.py: Script to fit the continuum flux around a line and normalize the spectrum around that line

import glob
import os
import numpy as np

from astropy.io import ascii, fits
from astropy.table import Table
from matplotlib import pyplot as plt

from plotting.plot_spectra import velo


def plot_cont(velos, fluxes, cont, degree, fitmask, rangemask, outname):
    """
    Function to plot the continuum flux

    Parameters
    ----------
    velos : np.ndarray
        Velocities (in km/s)

    fluxes : np.ndarray
        Fluxes (in erg/s/cm2/AA)

    cont : np.ndarray
        Continuum flux at velocities in rangemask

    degree : int
        Degree of the polynomial used for the continuum fit

    fitmask : np.ndarray
        Mask for data points used in the fitting

    rangemask : np.ndarray
        Mask for data points between minimum and maximum velocity

    outname : string
        Path and name of the plot

    Returns
    -------
    Plot with data fluxes and continuum fit
    """
    # create a plot
    fix, ax = plt.subplots()

    # plot the data
    plt.plot(velos, fluxes, c="k", alpha=0.7, label="data")

    # plot the continuum fit
    plt.plot(
        velos[rangemask],
        cont,
        c="tab:orange",
        label="cont. fit (" + str(degree) + ")",
    )

    # plot the data points that were used in the fitting
    plt.plot(
        velos[fitmask],
        fluxes[fitmask],
        "r.",
        markersize=2.5,
        alpha=0.6,
        label="fit points",
    )

    # finalize and save the plot
    plt.xlim(-510, 510)
    plt.ylim([0, 1.5 * np.median(fluxes[fitmask])])
    plt.xlabel(r"velocity (km s$^{-1}$)", fontsize=fs)
    plt.ylabel(r"flux ($erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$)", fontsize=fs)
    plt.legend(loc=3)
    plt.savefig(
        outname,
        bbox_inches="tight",
    )


def norm(velos, fluxes, lwave, windows, outname):
    """
    Function to normalize the spectrum around the given line

    Parameters
    ----------
    velos : np.ndarray
        Velocities (in km/s)

    fluxes : np.ndarray
        Fluxes at velos (in erg/s/cm2/Angstrom)

    lwave : float
        Wavelength (in Angstrom) of absorption line

    windows : Astropy Table
        Table with the continuum windows

    outname : string
        Path and name of the plot

    Returns
    -------
    degree : int
        Degree of the polynomial used for the continuum fit

    coefs : np.ndarray
        Array of polynomial coefficients of the continuum fit

    cont : np.ndarray
        Fitted continuum flux (in erg/s/cm2/Angstrom) at velos in a certain range

    norm_fluxes : np.ndarray
        Normalized flux (in erg/s/cm2/Angstrom) at velos in a certain range
    """
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
        fitmask = np.full(len(velos), False)
        for vmin, vmax in zip(vmins, vmaxs):
            fitmask += (velos > vmin) & (velos < vmax)
        # exclude NaNs from the fitting
        fitmask *= ~np.isnan(fluxes)

        # fit the continuum with a Legendre polynomial
        # obtain the degree of the polynomial to fit
        degree = windows["degree"][rows][0]
        if degree == 0:
            print(
                "Please, add polynomial degree for line with wavelength "
                + str(lwave)
                + " \u212B to the windows file."
            )
        else:
            coefs, diags = np.polynomial.legendre.legfit(
                velos[fitmask], fluxes[fitmask], degree, full=True
            )
            # obtain the sum of squared residuals of the fit
            sqres = diags[0][0]
            print("Sum of squared residuals: ", sqres)

            # evaluate the continuum at all velocities in the fitting range, i.e. from the minimum to the maximum velocity used in the fitting (including velocities in between, that were not used in the continuum fit)
            rangemask = (velos > np.min(vmins)) & (velos < np.max(vmaxs))
            cont = np.polynomial.legendre.legval(velos[rangemask], coefs)

            # normalize the spectrum to the continuum
            norm_fluxes = fluxes[rangemask] / cont

            # plot the continuum
            plot_cont(
                velos,
                fluxes,
                cont,
                degree,
                fitmask,
                rangemask,
                outname.replace(".pdf", str(degree) + ".pdf"),
            )

    return degree, coefs, cont, norm_fluxes


def norm_all(datapath, star, linewaves):
    """
    Function to normalize the spectrum around every line

    Parameters
    ----------
    datapath : string
        Path to the data files

    star : string
        Name of the star

    linewaves : np.ndarray
        Wavelengths of the lines

    Returns
    -------
    Plot with the continuum fit
    """
    # check if the file with windows to fit the continuum exists
    window_file = datapath + star + "/" + star + "_cont_win.dat"
    if not os.path.isfile(window_file):
        print("Please, create windows file for star " + star, window_file)
    else:
        # read the windows file
        windows = ascii.read(window_file)

        # create the Table to store the continuum coefficients
        ncols = max(windows["degree"]) + 1
        colnames = ["c" + str(i) for i in range(ncols)]
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

                    # normalize the spectrum around the line
                    degree, coefs, cont, norm_fluxes = norm(
                        velos,
                        fluxes,
                        lwave,
                        windows,
                        datapath
                        + star
                        + "/"
                        + star
                        + "_cont_fit_"
                        + str(int(lwave))
                        + "_.pdf",
                    )

                    # add polynomial coefficients to the table
                    if len(coefs) < ncols:
                        coefs = np.pad(coefs, (0, ncols - len(coefs)))
                    coeff_table.add_row(np.insert(coefs, 0, [lwave, degree]))

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
    for star in stars:
        norm_all(datapath, star, line_waves)
