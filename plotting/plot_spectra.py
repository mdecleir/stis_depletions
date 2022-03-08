# plot_spectra: different functions to create plots of the spectra
import glob
import os
import numpy as np

from astropy.io import fits
from matplotlib import pyplot as plt
from scipy import stats
from scipy.constants import c


def velo(waves, lwave):
    """
    Function to convert from wavelength to velocity space (in km/s)

    Parameters
    ----------
    waves : np.ndarray
        Array with wavelengths

    lwave : float
        Central rest wavelength of the line in vacuum

    Note: waves and lwave should have the same units

    Returns
    -------
    Array with velocities (in km/s)
    """

    return c / 1000 * (waves - lwave) / lwave


def rebin(waves, fluxes, rebin_fac):
    """
    Function to rebin a spectrum with a certain factor

    Parameters
    ----------
    waves : np.ndarray
        Numpy array with the wavelengths of the spectrum

    fluxes : np.ndarray
        Numpy array with the fluxes of the spectrum

    rebin_fac : int
        The factor by which to rebin the spectrum

    Returns
    -------
    The rebinned wavelengths and fluxes
    """
    # calculate the number of bins
    nbins = int(len(waves) / rebin_fac)

    # mask NaNs from the fluxes before rebinning, otherwise every bin with at least one NaN value would have NaN as the mean flux
    # also mask the corresponding wavelengths, otherwise all the wavelengths in a bin would be used to calculate the mean wavelength, while only the non-NaN fluxes would be used to calculate the mean flux
    mask = ~np.isnan(fluxes)
    waves = waves[mask]
    fluxes = fluxes[mask]

    # calculate the mean wavelength and mean flux in every bin
    # note that using np.nanmean as the statistic instead of masking the NaNs beforehand would result in using all the wavelengths in a bin to calculate the mean wavelength, while only using the non-NaN fluxes to calculate the mean flux
    new_waves, new_fluxes = stats.binned_statistic(
        waves,
        (waves, fluxes),
        statistic="mean",
        bins=nbins,
    )[0]

    return new_waves, new_fluxes


def plot_lines(datapath, star, rebin_fac=None):
    """
    Function to plot the absorption lines in velocity space

    Parameters
    ----------
    datapath : string
        The path of the data

    star : string
        The name of the star for which the lines will be plotted

    rebin_fac : int [default=None]
        The factor by which to rebin the spectrum. The rebinned spectrum will be overplotted on the full spectrum

    Returns
    -------
    Plot with all the absorption lines underneath each other
    """
    # obtain the data files with the merged spectra for this star
    file_list = glob.glob(datapath + star + "/*merged.fits")

    # define the lines of interest and their wavelengths
    line_names = np.array(
        [
            r"Mg\textsc{ii}",
            r"Mg\textsc{ii}",
            r"S\textsc{ii}",
            r"S\textsc{ii}",
            r"S\textsc{ii}",
            r"Si\textsc{ii}",
            r"Ni\textsc{ii}",
            r"C\textsc{ii}",
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
            1259.518,
            1260.4221,
            1317.217,
            1334.5323,
            1355.5977,
            1358.7730,
            1370.132,
            2249.8768,
            2260.7805,
            2325.4029,
            2335.123,
        ]
    )

    # check which absorption lines are in the wavelength range of each spectrum
    line_masks = []
    for filename in file_list:
        data = fits.getdata(filename, ext=1)
        waves = data["WAVELENGTH"]
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

        # rebin the spectrum if requested
        if rebin_fac is not None:
            binned_waves, binned_fluxes = rebin(waves, fluxes, rebin_fac)

        # for all the lines within the spectrum
        for lwave, lname in zip(line_waves[line_mask], line_names[line_mask]):
            # convert the wavelengths to velocities
            velos = velo(waves, lwave)

            # plot the spectrum
            x_mask = (velos > -100) & (velos < 100)
            ax[plotn].plot(velos[x_mask], fluxes[x_mask])

            # convert the rebinned wavelengths to velocities and overplot the rebinned spectrum if requested
            if rebin_fac is not None:
                binned_velos = velo(binned_waves, lwave)
                binned_x_mask = np.logical_and(
                    binned_velos > -100,
                    binned_velos < 100,
                    where=~np.isnan(binned_velos),
                )  # keep the NaNs in the velocities between -100 and 100, to create a gap in the plot where there is no data
                ax[plotn].plot(
                    binned_velos[binned_x_mask], binned_fluxes[binned_x_mask]
                )

            # add the name and wavelength of the line
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
    ax[-2].set_ylim(1.15e-12, 1.4e-12)
    ax[-1].set_ylim(1.3e-12, 1.55e-12)
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
        plot_lines(datapath, star, 5)
