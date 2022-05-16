import glob
import os
import numpy as np

from astropy.io import ascii, fits
from astropy.table import Table
from matplotlib import pyplot as plt

import norm
from plotting.plot_spectra import velo


def AOD(velos, norm_fluxes, log_flam, vmin, vmax, plot=False):
    """
    Function to calculate the apparent column density for the given line, using the apparent optical depth method

    Parameters
    ----------
    velos : np.ndarray
        Velocities (in km/s)

    norm_fluxes : np.ndarray
        Normalized fluxes at velos

    log_flam : float
        Logarithm of the oscillator strength of the absorption line log(f*lambda)

    vmin : float
        Minimum velocity to be included in the integration of the line profile

    vmax : float
        Maximum velocity to be included in the integration of the line profile

    plot : boolean [default=False]
        Whether or not to plot the column density profile

    Returns
    -------
    N : float
        Apparent column density for the given line
    """
    # calculate the apparent optical depth, Savage & Sembach 1991, eq. 5:
    # tau = ln(I0/Iobs), I0=continuum flux, Iobs=observed flux
    # tau = ln(1/Inorm), Inorm = Iobs/I0
    tau = np.log(1 / norm_fluxes)

    # calculate the apparent column density, Savage & Sembach 1991, eq. 10:
    # log(N) = log(tau) - log(f*lamda) + 14.576
    log_Nv = np.log10(tau) - log_flam + 14.576
    Nv = 10 ** log_Nv

    # plot the apparent column density profile if requested
    if plot:
        fig, ax = plt.subplots()
        fs = 18
        plt.plot(velos, log_Nv, c="k")
        plt.axvline(
            vmin,
            ls="--",
        )
        plt.axvline(
            vmax,
            ls="--",
        )
        plt.xlim(-100, 100)
        plt.xlabel(r"velocity (km s$^{-1}$)", fontsize=fs)
        plt.ylabel(r"log(N)", fontsize=fs)
        plt.show()

    # integrate (sum) the apparent column density profile over the relevant velocity range
    intmask = (velos > vmin) & (velos < vmax)
    # check if the velocities are equally distant
    dv_array = np.around(np.diff(velos[intmask]), 12)
    if np.all(dv_array == dv_array[0]):
        dv = np.median(np.diff(velos[intmask]))
        log_N = np.log10(np.nansum(Nv[intmask] * dv))
        print(velos[intmask], Nv[intmask])
        return log_N
    else:
        print("The velocities are not all equally distant.")


def AOD_all(datapath, star, linewaves, log_flams, vmin, vmax):
    """
    Function to calculate and save the column density for every line

    Parameters
    ----------
    datapath : string
        Path to the data files

    star : string
        Name of the star

    linewaves : np.ndarray
        Wavelengths (in Angstrom) of the lines

    log_flams : np.ndarray
        Array with the logarithm of the oscillator strengths of the absorption lines log(f*lambda) in linewaves

    vmin : float
        Minimum velocity to be included in the integration of the line profiles

    vmax : float
        Maximum velocity to be included in the integration of the line profiles

    Returns
    -------
    File with the column densities for all the lines
    """
    # check if the file with windows to fit the continuum exists
    window_file = datapath + star + "/" + star + "_cont_win.dat"
    if not os.path.isfile(window_file):
        print("Please, create windows file for star " + star, window_file)
    else:
        # read the continuum windows file
        windows = ascii.read(window_file)

        # create the Table to store the column densities
        coldens_table = Table(names=["line_wave", "log(Na)"])

        # obtain the data files with the merged spectra for this star
        file_list = glob.glob(datapath + star + "/*merged.fits")
        for filename in file_list:
            # obtain the data
            data = fits.getdata(filename, ext=1)
            waves = data["WAVELENGTH"]
            fluxes = data["FLUX"]

            # for all absorption lines
            for lwave, log_flam in zip(linewaves, log_flams):
                # check if the line is in the wavelength range of the spectrum
                if np.logical_and(lwave > waves[0], lwave < waves[-1]):
                    # convert the wavelengths to velocities
                    velos = velo(waves, lwave)

                    # normalize the spectrum around the line
                    degree, coefs, norm_fluxes, rangemask = norm.norm(
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

                    # calculate the apparent column density for this line
                    log_N = AOD(velos[rangemask], norm_fluxes, log_flam, vmin, vmax)

                    # add the column density to the table
                    coldens_table.add_row([lwave, log_N])
                    print(coldens_table)

        # save the table
        coldens_table.write(
            datapath + star + "/" + star + "_N.dat",
            format="ascii",
            overwrite=True,
        )


if __name__ == "__main__":
    # define the data path
    datapath = "/Users/mdecleir/Documents/Depletions/HST_data/"

    # define the star names
    stars = ["HD203938"]

    # plotting properties for uniform plots
    fs = 18

    # define the wavelenghts and oscillator strengths of the lines of interest
    line_waves = np.array([1239.9253, 1240.3947, 1355.5977, 2249.8768, 2260.7805])
    line_waves = np.array([1239.9253, 1240.3947, 1355.5977, 2249.8768, 2260.7805])
    log_flams = np.array([0.106, -0.355, -2.805, 0.612, 0.742])

    # normalize the spectrum around every line
    for star in stars:
        AOD_all(datapath, star, line_waves, log_flams, -26.2, 1.2)
