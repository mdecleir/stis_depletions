# plot_spectra: different functions to create plots of the spectra
import glob
import os
import numpy as np

from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.constants import c


# convert from wavelength to velocity space (in km/s)
def velo(waves, lwave):
    return c / 1000 * (waves - lwave) / lwave


# plot the absorption lines in velocity space
def plot_lines(datapath, star):
    """
    Function to plot the absorption lines in velocity space
    Parameters
    ----------
    datapath : string
        the path of the data
    star : string
        the name of the star for which the lines are plotted
    Returns
    -------
    Plot with all the absorption lines underneath each other
    """
    file_list = glob.glob(datapath + star + "/*merged.fits")

    # define the lines of interest and their wavelengths
    line_names = np.array(
        [
            r"Mg\textsc{ii}",
            r"Mg\textsc{ii}",
            r"S\textsc{ii}",
            r"S\textsc{ii}",
            r"Ni\textsc{ii}",
            r"O\textsc{i}",
            r"Cu\textsc{ii}",
            r"Ni\textsc{ii}",
            r"Fe\textsc{ii}",
            r"Fe\textsc{ii}",
            r"C\textsc{ii}",
            r"Si\textsc{ii}",
        ]
    )
    line_waves = np.array(
        [
            1239.9253,
            1240.3947,
            1250.578,
            1253.805,
            1317.217,
            1355.5977,
            1358.7730,
            1370.132,
            2249.8768,
            2260.7805,
            2325.4029,
            2335.123,
        ]
    )

    # check which data are available
    line_masks = []
    for filename in file_list:
        data = fits.getdata(filename, ext=1)
        waves = data["WAVELENGTH"]
        # check which absorption lines are in the wavelength range of this spectrum
        line_masks.append(np.logical_and(line_waves > waves[0], line_waves < waves[-1]))

    # create the figure
    nlines = np.sum(line_masks)
    if nlines > 0:
        fig, ax = plt.subplots(nlines, figsize=(5, nlines), sharex=True)

    # plot the absorption lines
    # initialize the plot number
    plotn = 0
    for filename, line_mask in zip(file_list, line_masks):
        # obtain the data
        data = fits.getdata(filename, ext=1)
        waves = data["WAVELENGTH"]
        fluxes = data["FLUX"]
        # for all the lines within the spectrum
        for lwave, lname in zip(line_waves[line_mask], line_names[line_mask]):
            # convert the wavelengths to velocities
            velos = velo(waves, lwave)

            # plot the spectrum
            x_mask = (velos > -100) & (velos < 100)
            ax[plotn].plot(velos[x_mask], fluxes[x_mask])
            ax[plotn].text(
                0.02,
                0.5,
                lname + " (" + str(int(lwave)) + "Å)",
                transform=ax[plotn].transAxes,
                fontsize=fs * 0.8,
                usetex=True,
            )
            ax[plotn].axhline(c="k", lw=1, ls="--")

            # increase the plot number
            plotn += 1

    # finalize and save the plot
    plt.xlim(-100, 100)
    plt.subplots_adjust(hspace=0)
    plt.xlabel("velocity (km/s)")
    fig.supylabel("normalized flux", fontsize=fs)
    plt.savefig(datapath + star + "/" + star + "_all_lines.pdf", bbox_inches="tight")


if __name__ == "__main__":
    fs = 16
    plt.rc("font", size=fs)
    plt.rc("xtick", top=True, direction="in", labelsize=fs * 0.8)
    plt.rc("ytick", right=True, direction="in", labelsize=fs * 0.8)
    plt.rc("legend", fontsize=fs * 0.6)
    datapath = "/Users/mdecleir/Documents/Depletions/HST_data/"

    stars = ["HD203938"]

    # plot all the absorption lines of the star, using the merged spectra
    for star in stars:
        plot_lines(datapath, star)
