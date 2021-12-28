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
    star = "HD203938"

    datalist = glob.glob(datapath + star + "/*x1d.fits")
    fix_blaze(datapath + star + "/", star, datalist)
