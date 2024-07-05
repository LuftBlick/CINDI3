from os import path as ospath
import os
from numpy import array, loadtxt, hstack, mean, zeros, ones, where, savetxt, nan, isnan
from glob import glob
from copy import copy
from matplotlib.dates import datestr2num, num2date
from datetime import date
from datetime import datetime
from copy import deepcopy

from sys import path as syspath
syspath.append("C:/Blick/scripts/auxilary")
import os
# Define the desired path
path = os.path.normpath("C:/Blick/src")
os.chdir(path)
from GetDataOfRoutineBlick import GetDataOfRoutine

from blick_params import blick_params
gp = blick_params()

from pylab import *

from matplotlib import rcParams

rcParams.update({'font.size': 6})
rcParams.update({'lines.markersize': 2})

ylimNO2 = [-10, 100]
ylimO3 = [-0.1, 1.0]
ylimO2O2 = [-100, 4000]
ylimH2O = [-1, 1]
ylimHCHO = [-1, 10]
ylimRMS = [1e-4, 1e-1]
ylimRI = [1e1, 5e5]
ylimCI = [1.0, 1.25]
ylimWVL0 = [-0.001, 0.001]
ylimOffs = [-0.01, 0.01]
ylimOffsLin = [-0.01, 0.01]
ylimRing = [-0.1, 0.1]

proctype2ref = {
    'ONLYL1': 0,
    'SUN': 2,
    'MOON': 3,
    'SKY': 4,
    'TARGET': 5,
    'PROFILE': 6,
    'ALMUCANTAR': 7,
    'LAMP': 8,
    'SPECIAL': 9
}

reftypeconf = {
    'Ref': 'DAILYREF',
    'SyntOPEN': 'DSREF',
    'SyntU340': 'DSREF',
    'ExtOPEN': 'DSREF',
    'ExtU340': 'DSREF',
    'MeasLow': 'SEQREF',
    'MeasHigh': 'SEQREF',
}

colAssignL2Fit ={
    'DOY': 'UT date and time for center-time of measurement, yyyymmddThhmmssZ (ISO 8601)',
    'UTC': 'Fractional days since 1-Jan-2000 midnight for center-time of measurement',
    'ACQT': 'Total duration of measurement set in seconds',
    'SZA': 'Solar zenith angle for center-time of measurement in degree',
    'SAA': 'Solar azimuth for center-time of measurement in degree, 0=north, increases clockwise',
    'VEA': 'Pointing zenith angle in degree, absolute or relative (see next column), 999=tracker not used',
    'VEAMODE': 'Zenith pointing mode: zenith angle is... 0=absolute, 1=relative to sun, 2=relative to moon',
    'VAA': 'Pointing azimuth in degree, increases clockwise, absolute (0=north) or relative (see next column), 999=tracker not used',
    'VAAMODE': 'Azimuth pointing mode: like zenith angle mode but also fixed scattering angles relative to sun (3) or moon (4)',
    'NO2_DSCD_294': 'Nitrogen dioxide slant column amount [moles per square meter], -9e99=fitting not successful',
    'NO2_DSCD_220': 'Nitrogen dioxide slant column amount [moles per square meter], -9e99=fitting not successful',
    'NO2_DSCD_294_Error': 'Independent uncertainty of nitrogen dioxide slant column amount [moles per square meter], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'NO2_DSCD_220_Error': 'Independent uncertainty of nitrogen dioxide slant column amount [moles per square meter], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'NO2_T': 'Nitrogen dioxide effective temperature [K]',
    'H2O_T': 'Water vapor effective temperature [K]',
    'O4_DSCD_293': 'Oxygen dimer slant column amount [moles squared per meter to the 5th], -9e99=fitting not successful',
    'O4_DSCD_223': 'Oxygen dimer slant column amount [moles squared per meter to the 5th], -9e99=fitting not successful',
    'O4_DSCD_293_Error': 'Independent uncertainty of oxygen dimer slant column amount [moles squared per meter to the 5th], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'O4_DSCD_223_Error': 'Independent uncertainty of oxygen dimer slant column amount [moles squared per meter to the 5th], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'O4_T': 'Oxygen dimer effective temperature [K]',
    'O3_DSCD_243': 'Ozone slant column amount [moles per square meter], -9e99=fitting not successful',
    'O3_DSCD_223': 'Ozone slant column amount [moles per square meter], -9e99=fitting not successful',
    'O3_DSCD_243_Error': 'Independent uncertainty of ozone slant column amount [moles per square meter], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'O3_DSCD_223_Error': 'Independent uncertainty of ozone slant column amount [moles per square meter], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'O3_T': 'Ozone effective temperature [K]',
    'H2O_DSCD_273': 'Water vapor slant column amount [moles per square meter], -9e99=fitting not successful',
    'H2O_DSCD_273_Error': 'Independent uncertainty of water vapor slant column amount [moles per square meter], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'HCHO_DSCD_298': 'Formaldehyde slant column amount [moles per square meter], -9e99=fitting not successful',
    'HCHO_DSCD_298_Error': 'Independent uncertainty of formaldehyde slant column amount [moles per square meter], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'HCHO_T': 'Formaldehyde effective temperature [K]',
    'SO2_DSCD_298': 'Sulfur dioxide [moles per square meter], -9e99=fitting not successful',
    'SO2_DSCD_298_Error': 'Independent uncertainty of sulfur dioxide slant column amount [moles per square meter], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'CHOCHO_DSCD_296': 'Glyoxal slant column amount [moles per square meter], -9e99=fitting not successful',
    'CHOCHO_DSCD_296_Error': 'Independent uncertainty of Glyoxal slant column amount [moles per square meter], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'CHOCHO_T': 'Glyoxal effective temperature [K]',
    'BrO_DSCD_223': 'Bromine oxide slant column amount [moles per square meter], -9e99=fitting not successful',
    'BrO_DSCD_223_Error': 'Independent uncertainty of bromine oxide slant column amount [moles per square meter], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'BrO_T': 'Bromine oxide effective temperature [K]',
    'HONO_DSCD_296': 'Nitrous acid slant column amount [moles per square meter], -9e99=fitting not successful',
    'HONO_DSCD_296_Error': 'Independent uncertainty of nitrous acid slant column amount [moles per square meter], -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting not successful',
    'HONO_T': 'Nitrous acid effective temperature [K]',
    'Ring': 'Fitted Ring spectrum',
    'Ring_Error': 'Independent uncertainty of fitted Ring spectrum',
    'RMS': 'rms of unweighted fitting residuals',
    'SHIFT': 'Retrieved total wavelength shift [nm], -9=no wavelength change fitting',
    'OFFSET': 'Offset polynomial coefficient, order 0',
}

colAssignCindi3 ={
    'DOY': 'Col {:02d}: DOY: Day of year 2024, start with 1.0 for January 1st, 00:00 UTC (days since 2023-12-31 00:00 UTC)',
    'UTC': 'Col {:02d}: UTC: UTC Time of day (hours)',
    'ACQT': 'Col {:02d}: ACQT: Total Acquisition Time (s)',
    'SZA': 'Col {:02d}: SZA: Solar Zenith Angle (degree)',
    'SAA': 'Col {:02d}: SAA: Solar Azimuth Angle (degree) North=0, East=90',
    'VEA': 'Col {:02d}: VEA: Viewing Elevation Angle (degree)',
    'VAA': 'Col {:02d}: VAA: Viewing Azimuth Angle (degree) North=0, East=90',
    'NO2_DSCD': 'Col {:02d}: NO2_DSCD_{}: (1E15 molec/cm2)',
    'NO2_DSCD_Error': 'Col {:02d}: NO2_DSCD_{}_Error: (1E15 molec/cm2)',
    'O4_DSCD': 'Col {:02d}: O4_DSCD_{}: (1E40 molec2/cm5)',
    'O4_DSCD_Error': 'Col {:02d}: O4_DSCD_{}_Error: (1E40 molec2/cm5)',
    'O3_DSCD': 'Col {:02d}: O3_DSCD_{}: (1E20 molec/cm2)',
    'O3_DSCD_Error': 'Col {:02d}: O3_DSCD_{}_Error: (1E20 molec/cm2)',
    'HCHO_DSCD': 'Col {:02d}: HCHO_DSCD_{}: (1E15 molec/cm2)',
    'HCHO_DSCD_Error':'Col {:02d}: HCHO_DSCD_{}_Error: (1E15 molec/cm2)',
    'BrO_DSCD': 'Col {:02d}: BrO_DSCD_{}: (1E15 molec/cm2)',
    'BrO_DSCD_Error': 'Col {:02d}: BrO_DSCD_{}_Error: (1E15 molec/cm2)',
    'SO2_DSCD': 'Col {:02d}: SO2_DSCD_{}: (1E15 molec/cm2)',
    'SO2_DSCD_Error': 'Col {:02d}: SO2_DSCD_{}_Error: (1E15 molec/cm2)',
    'CHOCHO_DSCD': 'Col {:02d}: CHOCHO_DSCD_{}: (1E15 molec/cm2)',
    'CHOCHO_DSCD_Error': 'Col {:02d}: CHOCHO_DSCD_{}_Error: (1E15 molec/cm2)',
    'HONO_DSCD': 'Col {:02d}: HONO_DSCD_{}: (1E15 molec/cm2)',
    'HONO_DSCD_Error': 'Col {:02d}: HONO_DSCD_{}_Error: (1E15 molec/cm2)',
    'H2O_DSCD': 'Col {:02d}: H2O_DSCD_{}: (1E20 molec/cm2)',
    'H2O_DSCD_Error': 'Col {:02d}: H2O_DSCD_{}_Error: (1E20 molec/cm2)',
    'Ring': 'Col {:02d}: Ring: (1)',
    'Ring_Error': 'Col {:02d}: Ring_Error: (1)',
    'RMS': 'Col {:02d}: RMS: Fit RMS in OD (1)',
    'SHIFT': 'Col {:02d}: SHIFT: Spectrum shift (nm), against FRS reference',
    'OFFSET': 'Col {:02d}: OFFSET: Intensity offset normalised by the mean intensity (1)',
    'INORM': 'Col {:02d}: INORM_{}:  Normalised Intensity counts/integration time at wavelength lambda {} nm (1/s)',
}

colAssignL1 = {
    'RTN': 'Two letter code of measurement routine',
    'RTNC': 'Routine count',
    'REPC': 'Repetition count',
    'INT': 'Integration time',
    'UTC': 'UT date and time for beginning of measurement',
    'NCYB': 'Number of bright count cycles',
    'SAT': 'Saturation index',
    'ACQT': 'Total duration of measurement set in seconds',
    'WVL0': 'Retrieved wavelength change, order 0',
    'WVL1': 'Retrieved wavelength change, order 1',
    'WVLm': 'Mean wavelength correction applied',
    'WVLfi': 'Wavelength change fitting result index',
    'FW1ep': 'Effective position of filterwheel #1',
    'FW2ep': 'Effective position of filterwheel #2',
    'INORM': 'L1 data for each pixel',
    'INORNu': 'Independent instrumental uncertainty of L1 data for each pixel',
}

gasRefTemps = {
    'NO2': 254.5,
    'O3': 225.0,
    'O4': 262.0,
}

prodMainProd = {
    'NO2VIS':'NO2_DSCD_294',
    'O4VIS':'O4_DSCD_293',
    'NO2VIS-SMALL':'NO2_DSCD_294',
    'CHOCHO':'CHOCHO_DSCD_296',
    'O3VIS':'O3_DSCD_223',
    'NO2UV':'NO2_DSCD_294',
    'O4UV':'O4_DSCD_293',
    'HCHO-WIDE':'HCHO_DSCD_298',
    'HCHO':'HCHO_DSCD_298',
    'HONO':'HONO_DSCD_296',
    'O3UV':'O3_DSCD_223',
    'BrO':'BrO_DSCD_223',
    'NO2EXT':'NO2_DSCD_294',
    'H2OEXT':'H2O_DSCD_273',
    'O3EXT':'O3_DSCD_223',
    'SO2UV':'SO2_DSCD_298',
}

refTypeSyn = {
    'Ref': 'DAILYREF',
    'RefFix': 'FIXREF',
    'MeasLow': 'SEQREF',
    'MeasHigh': 'SEQREF',
    'SyntRef': 'DSREF',
    'SyntOPEN': 'DSREF',
    'SyntU340': 'DSREF',
    'ExtRef': 'DSREF',
    'ExtOPEN': 'DSREF',
    'ExtU340': 'DSREF',
}

saveFmt ={
    'DOY': '%.7e',
    'UTC': '%.7e',
    'ACQT': '%.7e',
    'SZA': '%.7e',
    'SAA': '%.7e',
    'VEA': '%.7e',
    'VEAMODE': '%.7e',
    'VAA': '%.7e',
    'VAAMODE': '%.7e',
    'NO2_DSCD_294': '%.7e',
    'NO2_DSCD_220': '%.7e',
    'NO2_DSCD_294_Error': '%.7e',
    'NO2_DSCD_220_Error': '%.7e',
    'NO2_T': '%.7e',
    'H2O_T': '%.7e',
    'O4_DSCD_293': '%.7e',
    'O4_DSCD_223': '%.7e',
    'O4_DSCD_293_Error': '%.7e',
    'O4_DSCD_223_Error': '%.7e',
    'O4_T': '%.7e',
    'O3_DSCD_243': '%.7e',
    'O3_DSCD_223': '%.7e',
    'O3_DSCD_243_Error': '%.7e',
    'O3_DSCD_223_Error': '%.7e',
    'O3_T': '%.7e',
    'H2O_DSCD_273': '%.7e',
    'H2O_DSCD_273_Error': '%.7e',
    'HCHO_DSCD_298': '%.7e',
    'HCHO_DSCD_298_Error': '%.7e',
    'SO2_DSCD_298': '%.7e',
    'SO2_DSCD_298_Error': '%.7e',
    'CHOCHO_DSCD_296': '%.7e',
    'CHOCHO_DSCD_296_Error': '%.7e',
    'BrO_DSCD_223': '%.7e',
    'BrO_DSCD_223_Error': '%.7e',
    'HONO_DSCD_296': '%.7e',
    'HONO_DSCD_296_Error': '%.7e',
    'RMS': '%.7e',
    'SHIFT': '%.7e',
    'OFFSET': '%.7e',
}

def LatestFile(lPth, bIsL2Tot=False):
    lVersC = []
    lVersP1 = []
    lVersP2 = []
    for sPthC in lPth:
        sPrt = sPthC.partition('/')[-1].split('_')[-1].partition('.txt')[0]
        # > Get calibration file version
        if not bIsL2Tot:
            # lVersC.append(int(sPrt[:-4].split('c')[-1]))  # remove dependency on 'p'!
            sPrt2 = sPrt[:-4].split('c')[-1]
            sPrt2 = sPrt2.split('d')[0]
            lVersC.append(int(sPrt2))  # remove dependency on 'p'!
        # > Get processing version
        lVersP1.append(int(sPrt[-3:].replace('-', '.').split('.')[0]) * 100)  # * 100 to enhance importance of main version
        lVersP2.append(int(sPrt[-3:].replace('-', '.').split('.')[1]))
    # > Find maximum
    if not bIsL2Tot:
        lPth = [lPth[argmax(sum(asarray([lVersC, lVersP1, lVersP2]), axis=0))]]
    else:
        lPth = [lPth[argmax(sum(asarray([lVersP1, lVersP2]), axis=0))]]
    return lPth

def GetPths(sDate, sPan, sLoc, sInstituteC, sInstNum, sProcGas, sRef, sVers, sSCode, sFCode):
    # > File name CINDI2 format
    sPthCindi = '{}_PANDORA_{}_{}_{}_CINDI3_{}_v{}.asc'.format(sInstituteC, sInstNum, sProcGas, sRef, sDate, sVers)
    # > File name Blick format
    sPthBlick = '{}_{}_{}_L2Fit_f{}*.txt'.format(sPan, sLoc, sDate, sFCode)
    # > L1 file in Blick format
    sPthL1Blick = '{}_{}_{}_L1_s{}*.txt'.format(sPan, sLoc, sDate, sSCode)

    return sPthCindi, sPthBlick, sPthL1Blick

class CINDI3SemiBlind:

    def __init__(self, cols, iwvls, iwvlref, l1, l2fit, inp, pth, fmtColdDesc, ovrwVza, ovrwVaa, missValue):
        self.cols = cols
        self.iwvls = iwvls
        self.iwvlref = iwvlref
        self.l1 = l1
        self.l2fit = l2fit
        self.inp = inp
        self.fitPars = inp['fitPars']
        self.pth = pth
        self.fmtColdDesc = fmtColdDesc
        self.ovrwVza = ovrwVza
        self.ovrwVaa = ovrwVaa
        self.missValue = missValue

    def buildupColumns(self, prodMainProd):

        # loop columns
        dataCols = []
        descrCols = ['Data format:', '------------']

        ## fitted columns
        for icol, datavar in enumerate(self.cols):
            if datavar == 'DOY':
                # Extract components
                day_of_year = self.l2fit[datavar].time.dt.dayofyear
                hour = self.l2fit[datavar].time.dt.hour
                minute = self.l2fit[datavar].time.dt.minute
                second = self.l2fit[datavar].time.dt.second

                # Calculate the decimal day of the year
                decimal_day_of_year = day_of_year + (hour + minute / 60.0 + second / 3600.0) / 24.0
                dataCols.append(list(decimal_day_of_year.data))
                # write column description:
                descrCols.append(self.fmtColdDesc[datavar].format(icol + 1))
            elif datavar == 'UTC':
                # Extract hour and minute
                hours = self.l2fit[datavar].time.dt.hour
                minutes = self.l2fit[datavar].time.dt.minute
                seconds = self.l2fit[datavar].time.dt.second
                # Convert to decimal hour
                decimal_hour = hours + minutes / 60.0 + seconds / 3600.0

                dataCols.append(list(decimal_hour.data))
                # write column description:
                descrCols.append(self.fmtColdDesc[datavar].format(icol + 1))
            elif datavar == 'VEA':
                vza = self.l2fit['VEA'].data
                vza[self.l2fit['ZPM'].data == 1.] = self.l2fit['SZA'].data[self.l2fit['ZPM'].data == 1.]
                if self.ovrwVza:
                    vza = self.nominalAngles(self.ovrwVza)
                vea = 90. - vza
                dataCols.append(list(vea))
                # write column description:
                descrCols.append(self.fmtColdDesc[datavar].format(icol + 1))
            elif datavar == 'VAA':
                vaa = self.l2fit['VAA'].data
                vaa[self.l2fit['APM'].data == 1.] = self.l2fit['SAA'].data[self.l2fit['APM'].data == 1.]
                if self.ovrwVaa:
                    vaa = self.nominalAngles(self.ovrwVaa)
                dataCols.append(list(vaa))
                # write column description:
                descrCols.append(self.fmtColdDesc[datavar].format(icol + 1))
            elif '_Error' in datavar and len(datavar.split('_')) == 4:
                gas, Tx = datavar.split('_DSCD_')
                if gas == 'O4':
                    unitcf = 2.75*10**-38  # from mol**2/m**5 to molec**2/cm**5
                else:
                    unitcf = gp.unitcf[1][1]
                Tx = Tx[:-6]
                colDescKey = ''.join(datavar.split('_{}'.format(int(Tx))))
                # get the scale factor
                l, m, r = self.fmtColdDesc[colDescKey].rpartition('E')
                scl = float(l[-1] + m + r[:2])
                # extract main temperature value
                mainProd = [p for p in unique(unique(prodMainProd.values())).tolist() if p.startswith(gas)][0]
                _, _, Tmain = mainProd.split('_')
                Tmain = float(Tmain)
                # only include the "main temperature"
                if Tmain == float(Tx):
                    dataNew = self.l2fit[datavar].data / unitcf / scl
                else:
                    dataNew = array(ones(self.l2fit[datavar].data.shape[0])) * self.missValue
                dataCols.append(list(dataNew))
                # write column description:
                descrCols.append(self.fmtColdDesc[colDescKey].format(icol + 1, int(Tx)))
            elif '_DSCD_' in datavar:  # retrieved column amounts are (re)written
                gas, Tx = datavar.split('_DSCD_')
                if gas == 'O4':
                    unitcf = 2.75*10**-38  # from mol**2/m**5 to molec**2/cm**5
                else:
                    unitcf = gp.unitcf[1][1]
                Tx = float(Tx)
                colDescKey = datavar.split('_{}'.format(int(Tx)))[0]
                # check gas fitting method
                tfitmeth = self.fitPars['Gas temps'][self.fitPars['Fitted gases'] == gas]
                # get the scale factor
                l, m, r = self.fmtColdDesc[colDescKey].partition('E')
                scl = float(l[-1] + m + r[:2])
                if tfitmeth == 'FIT':  # split DSCDs to two temperatures
                    Tref = gasRefTemps[gas]
                    # get second temperature
                    for datavar2 in self.cols:
                        parts = datavar2.split('_')
                        if (len(parts) == 3) and (parts[0] == gas) and (datavar2 != datavar):
                            Ty = float(parts[2])
                            Tlowhigh = sorted(array([Tx, Ty]))
                    # convert DSCD and temp fit to two DSCD
                    DSCD_T = self.DSCDaTfit2DSCDs(self.l2fit[datavar].data, self.l2fit['{}_T'.format(gas)].data,
                                                  Tlowhigh[0], Tlowhigh[1], Tref, scl*unitcf)
                    dataCols.append(list(DSCD_T[Tx]))
                else:  # use DSCD
                    # extract main temperature value
                    mainProd = [p for p in unique(unique(prodMainProd.values())).tolist() if p.startswith(gas)][0]
                    _, _, Tmain = mainProd.split('_')
                    Tmain = float(Tmain)
                    # only include the "main temperature"
                    if Tmain == Tx:
                        dataNew = self.l2fit[datavar].data  / unitcf / scl
                    else:
                        dataNew = array(ones(self.l2fit[datavar].data.shape[0])) * self.missValue
                    dataCols.append(list(dataNew))
                # write column description:
                descrCols.append(self.fmtColdDesc[datavar.split('_{}'.format(int(Tx)))[0]].format(icol + 1, int(Tx)))
            elif '_T' in datavar:  # retrieved temperature is not used directly
                pass
            else:
                dataCols.append(list(self.l2fit[datavar].data))
                # write column description:
                descrCols.append(self.fmtColdDesc[datavar].format(icol + 1))
        dataCols = asarray(dataCols).T

        ## spectral columns
        ### reduce L1 to same values as L2Fit
        #### reduce to equal routines
        l1strip = self.l1.isel(time=isin(self.l1.RTNC, self.l2fit.RTNC))
        #### reduce to equal functional filter
        dtEqual = array([], dtype='datetime64[ns]')
        for rtnc in unique(self.l2fit.RTNC):
            repcL2fit = self.l2fit.REPC.isel(time=self.l2fit.RTNC == rtnc)
            l1stripRTNC = l1strip.isel(time=l1strip.RTNC == rtnc)
            dt = l1stripRTNC.isel(time=isin(l1stripRTNC.REPC.data, repcL2fit.data)).time.data
            dtEqual = append(dtEqual, dt)
        ### average
        l1stripavg = self.getAveragedSpectra(l1strip.sel(time=dtEqual))
        ### extract columns
        specCols = l1stripavg.sel(lam=self.iwvls, method='nearest').data
        # add reference wavelength to the top
        specCols = hstack((l1stripavg.sel(lam=self.iwvlref, method='nearest').data[:, None], specCols))
        # write column description
        for j, wvl in enumerate(append(self.iwvlref, self.iwvls)):  # first wavelength is reference wavelength for horizon scan. It should have no extension in the column name
            icol += 1
            if j == 0:
                descrCols.append(self.fmtColdDesc['INORM'].format(icol + 1, int(wvl), int(wvl)).replace("_{}".format(int(wvl)), ''))
            else:
                descrCols.append(self.fmtColdDesc['INORM'].format(icol + 1, int(wvl), int(wvl)))
        descrCols.append(' ')

        # stack columns
        allCols = hstack((dataCols, specCols))

        return allCols, descrCols

    def nominalAngles(self, nomAng):

        nza = []
        for rtns in unique(self.l2fit['RTNC']):
            idx = self.l2fit['RTNC'] == rtns
            nza = append(
                nza, nomAng[self.l2fit['RTN'].isel(time=idx).data[0]][:self.l2fit['REPC'].isel(time=idx).shape[0]]
            )
        return nza

    @staticmethod
    def DSCDaTfit2DSCDs(SC, T, Tlow, Thigh, Tref, scl):

        xT = (T - Tlow) / (Thigh - Tlow)  # is it Tlow instead of Tref?

        SChigh = SC * xT
        SChigh /= scl

        SClow = SC * (1 - xT)
        SClow /= scl

        return {Tlow: SClow, Thigh: SChigh}

    @staticmethod
    def getAveragedSpectra(l1strip, avg=1.0):
        lam_diff = diff(l1strip.INORM['lam'].values)
        mean_lam_diff = mean(lam_diff)
        window_size = int(round(avg / mean_lam_diff))

        return l1strip.INORM.rolling(lam=window_size, center=True).mean()

    def headerRetrSet(self):

        retrSet = [
            "Retrieval settings:",
            "-------------------",
            "Fitting Window: {}-{} nm".format(self.fitPars['WL-starts'], self.fitPars['WL-ends']),
            "Polynomial degree: {} ({} coefficients)".format(self.fitPars['npol'], int(self.fitPars['npol']) + 1),
            "Offset order: {} ({} coefficients)".format(self.fitPars['noffs'], int(self.fitPars['noffs']) + 1),
            "Wavelength correction order: {} ({} coefficients)".format(self.fitPars['nwlc'], int(self.fitPars['nwlc']) + 1),
            "Resolution correction order: {} ({} coefficients)".format(self.fitPars['nresc'], int(self.fitPars['nresc']) + 1),
            "Reference spectrum: average of zenith spectra between 11:30 and 11:40",
        ]

        for gas, src in zip(self.fitPars['Fitted gases'], self.fitPars['Gas sources']):
            retrSet.append('   {}: {}'.format(gas, src))

        if self.fitPars['Ring'] != 'NO':
            retrSet.append("   Ring: High Resolution calculation with QDOAS according to Chance and Spurr (1997) and normalized as in Wagner et al. (2009)",)

        retrSet.append(' ')

        return retrSet
    
    def headerGeneral(self, prodMainProd, refTypeSyn, comment):

        generalInfo = [
            "CAMPAIGNNAME: CINDI3",
            "SITE: Cabauw, The Netherlands (CESAR site)",
            "LATITUDE: 51.968 degree N",
            "LONGITUDE: 4.927 degree E",
            "ALTITUDE: -1 m asl",
            "INSTITUTE: {}".format(self.inp['institute']),
            "INSTRUMENTTYPE: PANDORA",
            "INSTRUMENTNUMBER: {}".format(self.inp['instNumber']),
            "DATAPRODUCT: {}".format(self.inp['prod']),
            "PRODUCTDSCD: {}".format(prodMainProd[self.inp['prod']]),
            "REFTYPE: {}".format(refTypeSyn[self.inp['refType']]),
            "Missing value: {}".format(self.missValue),
            "Retrieval code: BlickP (v1.8.62, 1 April 2024)",
            "Created by: Martin Tiefengraber (LuftBlick)",
            "Version: v{}".format(self.inp['prodVers']),
            "comment: {}".format(comment),
            " ",
        ]

        return generalInfo