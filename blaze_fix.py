# blaze_fix.py: script to correct STIS echelle data for misalignment in the blaze function
# more information can be found here:
# - https://stisblazefix.readthedocs.io/en/latest/
# - https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/stis/documentation/instrument-science-reports/_documents/2018_01.pdf

import stisblazefix
import glob


def fix_blaze(path, star, datalist):
    stisblazefix.fluxfix(datalist, path + star + "_blaze_fix_diagnostic.pdf", nxplot=5)


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
        datalist = glob.glob(datapath + star + "/*x1d.fits")
        if len(datalist) > 0:
            fix_blaze(datapath + star + "/", star, datalist)
        else:
            print(star + " has no x1d data available to perform a blaze fix.")
