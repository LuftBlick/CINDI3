from os import path as ospath
from numpy import array, loadtxt, hstack, mean, zeros, ones, where, savetxt, nan, isnan
from glob import glob
from copy import copy
from matplotlib.dates import datestr2num, num2date
from datetime import date
from datetime import datetime
from copy import deepcopy

from sys import path as syspath
syspath.append("C:/Blick/scripts/auxilary")

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

def LatestFile(lPth, bIsL2Tot=False):
    lVersC = []
    lVersP1 = []
    lVersP2 = []
    for sPthC in lPth:
        sPrt = sPthC.partition('\\')[-1].split('_')[-1].partition('.txt')[0]
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

def ConvertNominalVZA(par, a1Rtn, a1VZA, iSpec, doElev=False):

    a1VZACon = ones(a1VZA.shape) * nan

    dRtnVZA = deepcopy(par['dOvrwVZA'])

    #> double occurance of each measurement for s2 spectrometers (OPEN + U340 measure the same)
    if iSpec > 1:
        for sRtn in dRtnVZA.keys():
            lTmp = [2 * [c] for c in dRtnVZA[sRtn]]
            lTmp2 = []
            for c in lTmp:
                lTmp2 += c
            dRtnVZA[sRtn] = array(lTmp2)

    dRtn = {}
    for sRtn in dRtnVZA.keys():
        iRtn = array([ord(c)**2 for c in sRtn]).sum()
        dRtn[iRtn] = sRtn

    iRtnPrev = a1Rtn[0]
    iCnt = 0
    for iRtnCnt, iRtn in enumerate(a1Rtn):
        if iRtn in dRtn.keys():  # convert only mentioned routines
            iRtnLen = len(dRtnVZA[dRtn[iRtn]])
            if iRtnLen > 1:
                if iRtn != iRtnPrev:
                    a1VZACon[iRtnCnt] = dRtnVZA[dRtn[iRtn]][0]
                    iCnt = 1
                else:
                    if iCnt < iRtnLen:
                        a1VZACon[iRtnCnt] = dRtnVZA[dRtn[iRtn]][iCnt]
                        iCnt += 1
                    else:
                        a1VZACon[iRtnCnt] = dRtnVZA[dRtn[iRtn]][0]
                        iCnt = 1
            else:
                a1VZACon[iRtnCnt] = dRtnVZA[dRtn[iRtn]][0]
        else:  # otherwise use measured VZA
            a1VZACon[iRtnCnt] = a1VZA[iRtnCnt]
        iRtnPrev = iRtn

    # Convert to elevation angle
    if doElev:
        a1VZACon = 90. - a1VZACon

    return a1VZACon


def GetPths(sDate, sPan, sLoc, sInstituteC, sInstNum, sProcGas, sVers, sSCode, sFCode):
    # > File name CINDI2 format
    sPthCindi = '{}_MAXDOAS_{}_{}_CINDI2_{}_v{}.asc'.format(sInstituteC, sInstNum, sProcGas, sDate, sVers)
    # > File name Blick format
    sPthBlick = '{}_{}_{}_L2Fit_f{}*.txt'.format(sPan, sLoc, sDate, sFCode)
    # > L1 file in Blick format
    sPthL1Blick = '{}_{}_{}_L1_s{}*.txt'.format(sPan, sLoc, sDate, sSCode)

    return sPthCindi, sPthBlick, sPthL1Blick


def CINDI2SemiBlind_NO2vis(par, sInstituteC, sDate, sLoc, sInstNum, sPan, sProcGas, sSCode, sFCode, sVers, nHead, nHeadL1, iCcCol):

    #> Get paths
    sPthCindi, sPthBlick, sPthL1Blick = GetPths(sDate, sPan, sLoc, sInstituteC, sInstNum, sProcGas, sVers, sSCode, sFCode)

    # > Read Blick data
    sPthDat = glob(ospath.join(par['sL2FitPth'], sPthBlick))
    if len(sPthDat) > 0:
        sPthDat = LatestFile(sPthDat)
        for sPthBlickGlob in sPthDat:
            a2Dat = loadtxt(sPthBlickGlob,
                            converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                            skiprows=nHead)

        #> Get wavelength vector
        GD = GetDataOfRoutine(int(sPan.split('Pandora')[-1][:-2]), int(sPan.split('Pandora')[-1][-1]),
                              sLoc, par['sBlickRootPth'], par['sL0Pth'],
                              par['sOFPth'], par['sCFPth'], [sSCode, -1, -1], [-1, -1, -1])
        _, a1Wvl, _, panPars = GD.GetCalData(num2date(datestr2num(str(sDate))).replace(tzinfo=None))

        #> Read Blick L1 data
        a2Cc = None
        lPthBlickL1Glob = LatestFile(glob(ospath.join(par['sL1Pth'], sPthL1Blick)))
        for sPthBlickL1Glob in lPthBlickL1Glob:
            a2DatL1 = loadtxt(sPthBlickL1Glob,
                              converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                              skiprows=nHeadL1)
            bIdG = zeros((a2DatL1.shape[0]), dtype=bool)
            #> Reduce to L2Fit data
            for iRtnI in xrange(a2Dat.shape[0]):
                idx = where((a2DatL1[:, [3, 4]] == a2Dat[iRtnI, [3, 4]]).sum(axis=1) == 2)[0]
                if len(idx) > 0:
                    if len(idx) == 1:
                        bIdG[idx] = True
                    else:
                        assert True, 'routine doubled!'
                else:
                    assert True, 'corresponding routine not in L1 file!'
            assert (bIdG.shape[0] == a2DatL1.shape[0]) & (bIdG.sum() == a2Dat.shape[0]), 'L1 file filtering issue!'
            a2Cc = a2DatL1[bIdG, iCcCol-1:iCcCol-1+2048]
            a1Scl = a2DatL1[bIdG, iCcCol-1-3]
        fWvlRelInt = 440.
        a1RelInt = ones((a2Cc.shape[0]))
        fWvlColIdx0 = 425.
        fWvlColIdx1 = 440.
        a1ColIdx = ones((a2Cc.shape[0]))
        for iRtnI in xrange(a2Cc.shape[0]):
            #> Relative intensity
            bWvlIdx = (a1Wvl >=  fWvlRelInt - par['fWvlIntAvg']) & (a1Wvl <=  fWvlRelInt + par['fWvlIntAvg'])
            a1RelInt[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx] / a1Scl[iRtnI])
            #> Color index
            bWvlIdx0 = (a1Wvl >= fWvlColIdx0 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx0 + par['fWvlIntAvg'])
            bWvlIdx1 = (a1Wvl >= fWvlColIdx1 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx1 + par['fWvlIntAvg'])
            a1ColIdx[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx0] / a1Scl[iRtnI]) / mean(a2Cc[iRtnI, bWvlIdx1] / a1Scl[iRtnI])

        #> Write data
        ##> Meta info
        a2DatAct = array([]).reshape((a2Dat.shape[0], 0))
        for addCol in [3, 6, 8, 9, 12, 14]:
            if addCol == 3:
                ###> Decimal doy
                delta = (date(int(str(sDate)[:4]), 1, 1) - date(2000, 1, 1)).days
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] - delta))
                ###> Decimal hour of day
                a2DatAct = hstack((a2DatAct, a2Dat[:, [2]] % 1 * 24.))
            elif addCol == 12:
                ###> zenith viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VZA = a2Dat[:, [addCol - 1]]
                a1VZAMode = a2Dat[:, [addCol]]
                a1VZA[a1VZAMode == 0] = a1VZA[a1VZAMode == 0]
                a1VZA[a1VZAMode == 1] = a2Dat[(a1VZAMode == 1).squeeze(), 8 - 1] + a1VZA[a1VZAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVZA = -0.428
                    a1VZA[a1VZA > 5.] += fCorVZA
                ###> Convert to nomial elevation angles
                if par['doOvrwVZA']:
                    a1Elev = ConvertNominalVZA(par, a2Dat[:, 0], a1VZA, int(sPan.split('Pandora')[-1][-1]), doElev=True)
                else:
                    a1Elev = 90. - a1VZA
                a2DatAct = hstack((a2DatAct, a1Elev))
            elif addCol == 14:
                ###> azimuth viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VAA = a2Dat[:, [addCol - 1]]
                a1VAAMode = a2Dat[:, [addCol]]
                a1VAA[a1VZAMode == 0] = a1VAA[a1VAAMode == 0]
                a1VAA[a1VZAMode == 1] = a2Dat[(a1VAAMode == 1).squeeze(), 9 - 1] + a1VAA[a1VAAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVAA = 0.246
                    a1VAA[a1VZA > 5.] += fCorVAA
                a2DatAct = hstack((a2DatAct, a1VAA))
            else:
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))

        ##> Gases
        bOvrw = a2Dat[:, 15] > 2
        ###> NO2 DSC with Temp fit to two NO2 DSCs and rescale
        xs1 = 27
        uxs1 = 28
        xst2 = 29
        uxst2 = 30
        TX = 254.5
        T1 = 220.
        T2 = 298.
        ####> T = T2
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1))
                           * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + (((a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/((a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-15))
        ###> O2O2 DSC and rescale
        for addCol in [32, 33]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * 1.e-40))
            a2DatAct[bOvrw, -1] = -999
        ####> T = T1
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))
                           * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + ((1 - (a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/(1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-15))
        ###> O3 DSC and rescale
        for addCol in [24, 25]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * gp.unitcf[1][1] * 1.e-20))
            a2DatAct[bOvrw, -1] = -999
        ###> H2O DSC and rescale
        for addCol in [35, 36]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * 3.3428e22 * 1.e-23))
            a2DatAct[bOvrw, -1] = -999
        ###> Ring DSC and rescale
        for addCol in [40, 41]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999

        ##> Aux fit info
        for addCol in [18, 56]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999
        ##> Add relative intensity
        a1RelInt[a1RelInt == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1RelInt.reshape((a1RelInt.shape[0], 1))))
        ##> Add color index
        a1ColIdx[isinf(a1ColIdx)] = -999
        a1ColIdx[a1ColIdx == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1ColIdx.reshape((a1RelInt.shape[0], 1))))
        ##> Add intensity offset
        addCol = 54
        a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
        a2DatAct[bOvrw, -1] = -999
        ##> NO2 from retrieval output
        xs1 = 27
        uxs1 = 28
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [xs1 - 1]] * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [uxs1 - 1]] * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999

        # > Write header
        c = 0
        h = ''
        #> meta info
        h += '* NofHeaderlines: {}'.format(33) + '\n'
        h += '* NofColumns: {}'.format(a2DatAct.shape[1]) + '\n'
        h += '* Instrument identifier: {}_MAXDOAS'.format(sInstituteC) + '\n'
        h += '* Retrieval code: BlickP (v1.3.1, Feb 2017)' + '\n'
        h += '* Created by: Martin Tiefengraber' + '\n'
        h += '* Version: {}_v{}'.format(sProcGas, sVers) + '\n'; c += 1
        h += '* X-Axis (Col {}) = Day of year (DOY) 2016'.format(c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Time of day in hours (UTC)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Total Integration Time(s)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Zenith Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Azimuth Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Elevation Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Viewing Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        #> gases
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O4_DSCD_293 (1*10^40 molec2/cm5)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O4_DSCD_293_Error (1*10^40 molec2/cm5)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_220 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_220_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223 (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = H2O_DSCD_296 (1*10^23 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = H2O_DSCD_296_Error (1*10^23 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        #> aux fit info
        h += '* Y{}-Axis (Col {}) = Ring'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Ring_Error'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Fit RMS (in OD)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Spectrum shift (nm, against FRS reference)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Relative Intensity (counts/integration time @ {}nm)'.format(c-1, c, fWvlRelInt) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Colour index: ({} nm / {} nm)'.format(c-1, c, fWvlColIdx0, fWvlColIdx1) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = intensity offset with normalisation by I, I is the mean intensity in the spectral analysis windows, constant term'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '*DOY UTC Tint SZA SAA Elev Viewing_angle NO2_DSCD_298 NO2_DSCD_298_Error O4_DSCD_293 O4_DSCD_293_Error NO2_DSCD_220 NO2_DSCD_220_Error O3_DSCD_223 O3_DSCD_223_Error H2O_DSCD_296 H2O_DSCD_296_Error Ring Ring_Error RMS Spectrum_shift Intens(440) CI(425/440) offset_cst NO2_DSCD NO2_DSCD_Error'

        if par['dProdAna'][sProcGas]:
            f, ax = subplots(10, 2, figsize=(3, 8))
            idxNoon = argmin(a2DatAct[:, 3])
            fltAm = a2DatAct[:idxNoon, 5] > 85.
            fltPm = a2DatAct[idxNoon:, 5] > 85.
            f.suptitle(sPan + ', ' + sProcGas + ', ' + str(sDate))
            ax[0, 0].set_title('AM')
            ax[0, 1].set_title('PM')
            ylab = 'NO2 DSCD'
            ylim = ylimNO2
            scT1 = 7
            scT2 = 11
            sc = 24
            usc = 25
            row = 0
            ax[row, 0].plot(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, scT1][fltAm], 'o', label='298')
            ax[row, 0].plot(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, scT2][fltAm], 'o', label='220')
            ax[row, 0].errorbar(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], yerr=a2DatAct[:idxNoon, usc][fltAm], fmt='o', label='fit', capsize=0)
            ax[row, 1].plot(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, scT1][fltPm], 'o', label='298')
            ax[row, 1].plot(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, scT2][fltPm], 'o', label='220')
            ax[row, 1].errorbar(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], yerr=a2DatAct[idxNoon:, usc][fltPm], fmt='o', label='fit', capsize=0)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ax[row, 0].set_ylabel(ylab)
            ylab = 'O3 DSCD'
            ylim = ylimO3
            sc = 13
            usc = 14
            row = 1
            ax[row, 0].errorbar(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], yerr=a2DatAct[:idxNoon, usc][fltAm], fmt='o', capsize=0)
            ax[row, 1].errorbar(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], yerr=a2DatAct[idxNoon:, usc][fltPm], fmt='o', capsize=0)
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'O2O2 DSCD'
            ylim = ylimO2O2
            sc = 9
            usc = 10
            row = 2
            ax[row, 0].errorbar(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], yerr=a2DatAct[:idxNoon, usc][fltAm], fmt='o', capsize=0)
            ax[row, 1].errorbar(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], yerr=a2DatAct[idxNoon:, usc][fltPm], fmt='o', capsize=0)
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'H2O DSCD'
            ylim = ylimH2O
            sc = 15
            usc = 16
            row = 3
            ax[row, 0].errorbar(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], yerr=a2DatAct[:idxNoon, usc][fltAm], fmt='o', capsize=0)
            ax[row, 1].errorbar(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], yerr=a2DatAct[idxNoon:, usc][fltPm], fmt='o', capsize=0)
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'Ring DSCD'
            ylim = ylimRing
            sc = 17
            usc = 18
            row = 4
            ax[row, 0].errorbar(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], yerr=a2DatAct[:idxNoon, usc][fltAm], fmt='o', capsize=0)
            ax[row, 1].errorbar(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], yerr=a2DatAct[idxNoon:, usc][fltPm], fmt='o', capsize=0)
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'RMS'
            ylim = ylimRMS
            sc = 19
            row = 5
            ax[row, 0].semilogy(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], 'o')
            ax[row, 1].semilogy(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], 'o')
            # [ax[row, i].yaxis.set_scale('log') for i in xrange(ax.shape[1])]
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'WVL shift'
            ylim = ylimWVL0
            sc = 20
            row = 6
            ax[row, 0].plot(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], 'o')
            ax[row, 1].plot(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], 'o')
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'Offset'
            ylim = ylimOffs
            sc = 23
            row = 7
            ax[row, 0].plot(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], 'o')
            ax[row, 1].plot(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], 'o')
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'RI'
            ylim = ylimRI
            sc = 21
            row = 8
            ax[row, 0].semilogy(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], 'o')
            ax[row, 1].semilogy(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], 'o')
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'CI'
            ylim = ylimCI
            sc = 22
            row = 9
            ax[row, 0].plot(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], 'o')
            ax[row, 1].plot(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], 'o')
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]

            [ax[i, 0].set_xlim([100, 40]) for i in xrange(ax.shape[0])]
            [ax[i, 1].set_xlim([40, 100]) for i in xrange(ax.shape[0])]

            [ax[i, 1].get_yaxis().set_visible(False) for i in xrange(ax.shape[0])]
            [ax[i, j].get_xaxis().set_visible(False) for i in xrange(ax.shape[0]-1) for j in xrange(ax.shape[1])]
            # f.subplots_adjust(vspace=0)

            f.savefig(sPan + sProcGas + '_' + str(sDate) + '.png', dpi=300., bbox_inches='tight')
            close()

        #> save file
        sFmtMI = ['%.5f', '%.5f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'] # meta info
        sFmtDat = ['%.4e' for _ in xrange(a2DatAct.shape[1] - len(sFmtMI))]
        savetxt(ospath.join(par['sCampPth'], sPthCindi), a2DatAct, fmt=sFmtMI+sFmtDat, comments='', header=h)
    else:
        print '        ... file "{}" not available!'.format(sPthBlick)


def CINDI2SemiBlind_NO2visSmall(par, sInstituteC, sDate, sLoc, sInstNum, sPan, sProcGas, sSCode, sFCode, sVers, nHead, nHeadL1, iCcCol):

    #> Get paths
    sPthCindi, sPthBlick, sPthL1Blick = GetPths(sDate, sPan, sLoc, sInstituteC, sInstNum, sProcGas, sVers, sSCode, sFCode)

    # > Read Blick data
    sPthDat = glob(ospath.join(par['sL2FitPth'], sPthBlick))
    if len(sPthDat) > 0:
        sPthDat = LatestFile(sPthDat)
        for sPthBlickGlob in sPthDat:
            a2Dat = loadtxt(sPthBlickGlob,
                            converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                            skiprows=nHead)

        #> Get wavelength vector
        GD = GetDataOfRoutine(int(sPan.split('Pandora')[-1][:-2]), int(sPan.split('Pandora')[-1][-1]),
                              sLoc, par['sBlickRootPth'], par['sL0Pth'],
                              par['sOFPth'], par['sCFPth'], [sSCode, -1, -1], [-1, -1, -1])
        _, a1Wvl, _, panPars = GD.GetCalData(num2date(datestr2num(str(sDate))).replace(tzinfo=None))

        #> Read Blick L1 data
        a2Cc = None
        lPthBlickL1Glob = LatestFile(glob(ospath.join(par['sL1Pth'], sPthL1Blick)))
        for sPthBlickL1Glob in lPthBlickL1Glob:
            a2DatL1 = loadtxt(sPthBlickL1Glob,
                              converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                              skiprows=nHeadL1)
            bIdG = zeros((a2DatL1.shape[0]), dtype=bool)
            #> Reduce to L2Fit data
            for iRtnI in xrange(a2Dat.shape[0]):
                idx = where((a2DatL1[:, [3, 4]] == a2Dat[iRtnI, [3, 4]]).sum(axis=1) == 2)[0]
                if len(idx) > 0:
                    if len(idx) == 1:
                        bIdG[idx] = True
                    else:
                        assert True, 'routine doubled!'
                else:
                    assert True, 'corresponding routine not in L1 file!'
            assert (bIdG.shape[0] == a2DatL1.shape[0]) & (bIdG.sum() == a2Dat.shape[0]), 'L1 file filtering issue!'
            a2Cc = a2DatL1[bIdG, iCcCol-1:iCcCol-1+2048]
            a1Scl = a2DatL1[bIdG, iCcCol-1-3]
        fWvlRelInt = 440.
        a1RelInt = ones((a2Cc.shape[0]))
        fWvlColIdx0 = 412.
        fWvlColIdx1 = 440.
        a1ColIdx = ones((a2Cc.shape[0]))
        for iRtnI in xrange(a2Cc.shape[0]):
            #> Relative intensity
            bWvlIdx = (a1Wvl >=  fWvlRelInt - par['fWvlIntAvg']) & (a1Wvl <=  fWvlRelInt + par['fWvlIntAvg'])
            a1RelInt[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx] / a1Scl[iRtnI])
            #> Color index
            bWvlIdx0 = (a1Wvl >= fWvlColIdx0 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx0 + par['fWvlIntAvg'])
            bWvlIdx1 = (a1Wvl >= fWvlColIdx1 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx1 + par['fWvlIntAvg'])
            a1ColIdx[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx0] / a1Scl[iRtnI]) / mean(a2Cc[iRtnI, bWvlIdx1] / a1Scl[iRtnI])

        #> Write data
        ##> Meta info
        a2DatAct = array([]).reshape((a2Dat.shape[0], 0))
        for addCol in [3, 6, 8, 9, 12, 14]:
            if addCol == 3:
                ###> Decimal doy
                delta = (date(int(str(sDate)[:4]), 1, 1) - date(2000, 1, 1)).days
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] - delta))
                ###> Decimal hour of day
                a2DatAct = hstack((a2DatAct, a2Dat[:, [2]] % 1 * 24.))
            elif addCol == 12:
                ###> zenith viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VZA = a2Dat[:, [addCol - 1]]
                a1VZAMode = a2Dat[:, [addCol]]
                a1VZA[a1VZAMode == 0] = a1VZA[a1VZAMode == 0]
                a1VZA[a1VZAMode == 1] = a2Dat[(a1VZAMode == 1).squeeze(), 8 - 1] + a1VZA[a1VZAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVZA = -0.428
                    a1VZA[a1VZA > 5.] += fCorVZA
                ###> Convert to nomial elevation angles
                if par['doOvrwVZA']:
                    a1Elev = ConvertNominalVZA(par, a2Dat[:, 0], a1VZA, int(sPan.split('Pandora')[-1][-1]), doElev=True)
                else:
                    a1Elev = 90. - a1VZA
                a2DatAct = hstack((a2DatAct, a1Elev))
            elif addCol == 14:
                ###> azimuth viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VAA = a2Dat[:, [addCol - 1]]
                a1VAAMode = a2Dat[:, [addCol]]
                a1VAA[a1VZAMode == 0] = a1VAA[a1VAAMode == 0]
                a1VAA[a1VZAMode == 1] = a2Dat[(a1VAAMode == 1).squeeze(), 9 - 1] + a1VAA[a1VAAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVAA = 0.246
                    a1VAA[a1VZA > 5.] += fCorVAA
                a2DatAct = hstack((a2DatAct, a1VAA))
            else:
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))

        ##> Gases
        bOvrw = a2Dat[:, 15] > 2
        ###> NO2 DSC with Temp fit to two NO2 DSCs and rescale
        ####> True
        xs1 = 27
        uxs1 = 28
        xst2 = 29
        uxst2 = 30
        TX = 254.5
        T1 = 220.
        T2 = 298.
        ####> T = T2
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1))
                           * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + (((a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/((a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-15))
        ###> O2O2 DSC and rescale
        for addCol in [32, 33]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * 1.e-40))
            a2DatAct[bOvrw, -1] = -999
        ####> T = T1
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))
                           * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + ((1 - (a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/(1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-15))
        ###> O3 DSC and rescale
        for addCol in [24, 25]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * gp.unitcf[1][1] * 1.e-20))
            a2DatAct[bOvrw, -1] = -999
        ###> H2O DSC and rescale
        for addCol in [35, 36]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * 3.3428e22 * 1.e-23))
            a2DatAct[bOvrw, -1] = -999
        ###> Ring DSC and rescale
        for addCol in [40, 41]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999

        ##> Aux fit info
        for addCol in [18, 54]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999
        ##> Add relative intensity
        a1RelInt[a1RelInt == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1RelInt.reshape((a1RelInt.shape[0], 1))))
        ##> Add color index
        a1ColIdx[isinf(a1ColIdx)] = -999
        a1ColIdx[a1ColIdx == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1ColIdx.reshape((a1RelInt.shape[0], 1))))
        ##> Add intensity offset
        addCol = 52
        a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
        a2DatAct[bOvrw, -1] = -999
        ##> NO2 from retrieval output
        xs1 = 27
        uxs1 = 28
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [xs1 - 1]] * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [uxs1 - 1]] * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999

        # > Write header
        c = 0
        h = ''
        #> meta info
        h += '* NofHeaderlines: {}'.format(33) + '\n'
        h += '* NofColumns: {}'.format(a2DatAct.shape[1]) + '\n'
        h += '* Instrument identifier: {}_MAXDOAS'.format(sInstituteC) + '\n'
        h += '* Retrieval code: BlickP (v1.3.1, Feb 2017)' + '\n'
        h += '* Created by: Martin Tiefengraber' + '\n'
        h += '* Version: {}_v{}'.format(sProcGas, sVers) + '\n'; c += 1
        h += '* X-Axis (Col {}) = Day of year (DOY) 2016'.format(c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Time of day in hours (UTC)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Total Integration Time(s)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Zenith Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Azimuth Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Elevation Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Viewing Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        #> gases
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O4_DSCD_293 (1*10^40 molec2/cm5)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O4_DSCD_293_Error (1*10^40 molec2/cm5)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_220 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_220_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223 (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = H2O_DSCD_296 (1*10^23 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = H2O_DSCD_296_Error (1*10^23 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        #> aux fit info
        h += '* Y{}-Axis (Col {}) = Ring'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Ring_Error'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Fit RMS (in OD)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Spectrum shift (nm, against FRS reference)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Relative Intensity (counts/integration time @ {}nm)'.format(c-1, c, fWvlRelInt) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Colour index: ({} nm / {} nm)'.format(c-1, c, fWvlColIdx0, fWvlColIdx1) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = intensity offset with normalisation by I, I is the mean intensity in the spectral analysis windows, constant term'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '*DOY UTC Tint SZA SAA Elev Viewing_angle NO2_DSCD_298 NO2_DSCD_298_Error O4_DSCD_293 O4_DSCD_293_Error NO2_DSCD_220 NO2_DSCD_220_Error O3_DSCD_223 O3_DSCD_223_Error H2O_DSCD_296 H2O_DSCD_296_Error Ring Ring_Error RMS Spectrum_shift Intens(440) CI(412/440) offset_cst NO2_DSCD NO2_DSCD_Error'

        #> save file
        sFmtMI = ['%.5f', '%.5f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'] # meta info
        sFmtDat = ['%.4e' for _ in xrange(a2DatAct.shape[1] - len(sFmtMI))]
        savetxt(ospath.join(par['sCampPth'], sPthCindi), a2DatAct, fmt=sFmtMI+sFmtDat, comments='', header=h)
    else:
        print '        ... file "{}" not available!'.format(sPthBlick)


def CINDI2SemiBlind_NO2uv(par, sInstituteC, sDate, sLoc, sInstNum, sPan, sProcGas, sSCode, sFCode, sVers, nHead, nHeadL1, iCcCol):

    #> Get paths
    sPthCindi, sPthBlick, sPthL1Blick = GetPths(sDate, sPan, sLoc, sInstituteC, sInstNum, sProcGas, sVers, sSCode, sFCode)

    # > Read Blick data
    sPthDat = glob(ospath.join(par['sL2FitPth'], sPthBlick))
    if len(sPthDat) > 0:
        sPthDat = LatestFile(sPthDat)
        for sPthBlickGlob in sPthDat:
            a2Dat = loadtxt(sPthBlickGlob,
                            converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                            skiprows=nHead)
            # > Read VIS data for EH routine
            sPthDat = glob(sPthBlickGlob.partition('fC')[0] + 'fC07Oc*')
            sPthBlickGlob = LatestFile(sPthDat)[0]
            a2DatEH = loadtxt(sPthBlickGlob,
                              converters={0: lambda s: array([ord(c) ** 2 for c in s]).sum(), 1: datestr2num},
                              skiprows=103)
            bIsEH = a2DatEH[:, 0] == array([ord(c) ** 2 for c in 'EH']).sum()
            a2DatEH = a2DatEH[bIsEH, :a2Dat.shape[1]]  # extract only EH routines and cut by length of original data file
            if len(a2DatEH) > 0:
                # a2DatEH[:, 15:] = -999
                # > Merge UV with VIS data
                iMergPoint = where(a2Dat[:, 1] >= a2DatEH[-1, 1])[0][0]
                a2Dat = vstack((a2Dat[:iMergPoint, :], a2DatEH, a2Dat[iMergPoint:, :]))
                bIsEH = a2Dat[:, 0] == array([ord(c) ** 2 for c in 'EH']).sum()

        #> Get wavelength vector
        GD = GetDataOfRoutine(int(sPan.split('Pandora')[-1][:-2]), int(sPan.split('Pandora')[-1][-1]),
                              sLoc, par['sBlickRootPth'], par['sL0Pth'],
                              par['sOFPth'], par['sCFPth'], [sSCode, -1, -1], [-1, -1, -1])
        _, a1Wvl, _, panPars = GD.GetCalData(num2date(datestr2num(str(sDate))).replace(tzinfo=None))

        #> Read Blick L1 data
        a2Cc = None
        lPthBlickL1Glob = LatestFile(glob(ospath.join(par['sL1Pth'], sPthL1Blick)))
        for sPthBlickL1Glob in lPthBlickL1Glob:
            a2DatL1 = loadtxt(sPthBlickL1Glob,
                              converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                              skiprows=nHeadL1)
            bIdG = zeros((a2DatL1.shape[0]), dtype=bool)
            for iRtnI in xrange(a2Dat.shape[0]):
                idx = where((a2DatL1[:, [3, 4]] == a2Dat[iRtnI, [3, 4]]).sum(axis=1) == 2)[0]
                if len(idx) > 0:
                    if len(idx) == 1:
                        bIdG[idx] = True
                    else:
                        assert True, 'routine doubled!'
                else:
                    assert True, 'corresponding routine not in L1 file!'
            assert (bIdG.shape[0] == a2DatL1.shape[0]) & (bIdG.sum() == a2Dat.shape[0]), 'L1 file filtering issue!'
            a2Cc = a2DatL1[bIdG, iCcCol-1:iCcCol-1+2048]
            a1Scl = a2DatL1[bIdG, iCcCol-1-3]
        fWvlRelInt = 340.
        a1RelInt = ones((a2Cc.shape[0]))
        fWvlColIdx0 = 340.
        fWvlColIdx1 = 370.
        a1ColIdx = ones((a2Cc.shape[0]))
        for iRtnI in xrange(a2Cc.shape[0]):
            #> Relative intensity
            bWvlIdx = (a1Wvl >=  fWvlRelInt - par['fWvlIntAvg']) & (a1Wvl <=  fWvlRelInt + par['fWvlIntAvg'])
            a1RelInt[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx] / a1Scl[iRtnI])
            #> Color index
            bWvlIdx0 = (a1Wvl >= fWvlColIdx0 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx0 + par['fWvlIntAvg'])
            bWvlIdx1 = (a1Wvl >= fWvlColIdx1 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx1 + par['fWvlIntAvg'])
            a1ColIdx[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx0] / a1Scl[iRtnI]) / mean(a2Cc[iRtnI, bWvlIdx1] / a1Scl[iRtnI])

        #> Write data
        ##> Meta info
        a2DatAct = array([]).reshape((a2Dat.shape[0], 0))
        for addCol in [3, 6, 8, 9, 12, 14]:
            if addCol == 3:
                ###> Decimal doy
                delta = (date(int(str(sDate)[:4]), 1, 1) - date(2000, 1, 1)).days
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] - delta))
                ###> Decimal hour of day
                a2DatAct = hstack((a2DatAct, a2Dat[:, [2]] % 1 * 24.))
            elif addCol == 12:
                ###> zenith viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VZA = a2Dat[:, [addCol - 1]]
                a1VZAMode = a2Dat[:, [addCol]]
                a1VZA[a1VZAMode == 0] = a1VZA[a1VZAMode == 0]
                a1VZA[a1VZAMode == 1] = a2Dat[(a1VZAMode == 1).squeeze(), 8 - 1] + a1VZA[a1VZAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVZA = -0.428
                    a1VZA[a1VZA > 5.] += fCorVZA
                ###> Convert to nomial elevation angles
                if par['doOvrwVZA']:
                    a1Elev = ConvertNominalVZA(par, a2Dat[:, 0], a1VZA, int(sPan.split('Pandora')[-1][-1]), doElev=True)
                else:
                    a1Elev = 90. - a1VZA
                a2DatAct = hstack((a2DatAct, a1Elev))
            elif addCol == 14:
                ###> azimuth viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VAA = a2Dat[:, [addCol - 1]]
                a1VAAMode = a2Dat[:, [addCol]]
                a1VAA[a1VZAMode == 0] = a1VAA[a1VAAMode == 0]
                a1VAA[a1VZAMode == 1] = a2Dat[(a1VAAMode == 1).squeeze(), 9 - 1] + a1VAA[a1VAAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVAA = 0.246
                    a1VAA[a1VZA > 5.] += fCorVAA
                a2DatAct = hstack((a2DatAct, a1VAA))
            else:
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))

        ##> Gases
        bOvrw = a2Dat[:, 15] > 2
        if bIsEH.sum() > 0:
            bOvrw |= bIsEH
        ###> NO2 DSC with Temp fit to two NO2 DSCs and rescale
        ####> True
        xs1 = 29
        uxs1 = 30
        xst2 = 31
        uxst2 = 32
        TX = 254.5
        T1 = 220.
        T2 = 298.
        ####> T = T2
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1))
                           * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + (((a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/((a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-15))
        ###> O2O2 DSC and rescale
        for addCol in [34, 35]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * 1.e-40))
            a2DatAct[bOvrw, -1] = -999
        ####> T = T1
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))
                           * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + ((1 - (a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/(1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-15))
        ###> O3 DSC with Temp fit to two NO2 DSCs and rescale
        ####> True
        xs1 = 24
        uxs1 = 25
        xst2 = 26
        uxst2 = 27
        TX = 225.
        T1 = 223.
        T2 = 243.
        ####> T = T1
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))
                           * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + ((1 - (a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/(1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-20))
        ####> T = T2
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1))
                           * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + (((a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/((a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-20))
        ###> BrO DSC and rescale
        for addCol in [40, 41]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * 1.e-15))
            a2DatAct[bOvrw, -1] = -999
        ###> HCHO DSC and rescale
        for addCol in [37, 38]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * gp.unitcf[1][1] * 1.e-15))
            a2DatAct[bOvrw, -1] = -999
        ###> Ring DSC and rescale
        for addCol in [45, 46]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999

        ##> Aux fit info
        for addCol in [18, 61]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999
        ##> Add relative intensity
        a1RelInt[a1RelInt == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1RelInt.reshape((a1RelInt.shape[0], 1))))
        ##> Add color index
        a1ColIdx[isinf(a1ColIdx)] = -999
        a1ColIdx[a1ColIdx == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1ColIdx.reshape((a1RelInt.shape[0], 1))))
        ##> Add intensity offset
        addCol = 59  # constant term
        a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
        a2DatAct[bOvrw, -1] = -999
        ##> NO2 from retrieval output
        xs1 = 29
        uxs1 = 30
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [xs1 - 1]] * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [uxs1 - 1]] * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        ##> O3 from retrieval output
        xs1 = 24
        uxs1 = 25
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [xs1 - 1]] * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [uxs1 - 1]] * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999

        # > Write header
        c = 0
        h = ''
        #> meta info
        h += '* NofHeaderlines: {}'.format(39) + '\n'
        h += '* NofColumns: {}'.format(a2DatAct.shape[1]) + '\n'
        h += '* Instrument identifier: {}_MAXDOAS'.format(sInstituteC) + '\n'
        h += '* Retrieval code: BlickP (v1.3.1, Feb 2017)' + '\n'
        h += '* Created by: Martin Tiefengraber' + '\n'
        h += '* Version: {}_v{}'.format(sProcGas, sVers) + '\n'; c += 1
        h += '* X-Axis (Col {}) = Day of year (DOY) 2016'.format(c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Time of day in hours (UTC)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Total Integration Time(s)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Zenith Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Azimuth Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Elevation Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Viewing Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        #> gases
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O4_DSCD_293 (1*10^40 molec2/cm5)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O4_DSCD_293_Error (1*10^40 molec2/cm5)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_220 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_220_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223 (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_243 (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_243_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = BrO_DSCD_223 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = BrO_DSCD_223_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = HCHO_DSCD_297 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = HCHO_DSCD_297_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        #> aux fit info
        h += '* Y{}-Axis (Col {}) = Ring'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Ring_Error'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Fit RMS (in OD)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Spectrum shift (nm, against FRS reference)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Relative Intensity (counts/integration time @ {}nm)'.format(c-1, c, fWvlRelInt) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Colour index: ({} nm / {} nm)'.format(c-1, c, fWvlColIdx0, fWvlColIdx1) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = intensity offset with normalisation by I, I is the mean intensity in the spectral analysis windows, constant term'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '*DOY UTC Tint SZA SAA Elev Viewing_angle NO2_DSCD_298 NO2_DSCD_298_Error O4_DSCD_293 O4_DSCD_293_Error NO2_DSCD_220 NO2_DSCD_220_Error O3_DSCD_223 O3_DSCD_223_Error O3_DSCD_243 O3_DSCD_243_Error BrO_DSCD_223 BrO_DSCD_223_Error HCHO_DSCD_297 HCHO_DSCD_297_Error Ring Ring_Error RMS Spectrum_shift Intens(340) CI(340/370) offset_cst NO2_DSCD NO2_DSCD_Error O3_DSCD O3_DSCD_Error'

        #> save file
        sFmtMI = ['%.5f', '%.5f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'] # meta info
        sFmtDat = ['%.4e' for _ in xrange(a2DatAct.shape[1] - len(sFmtMI))]
        savetxt(ospath.join(par['sCampPth'], sPthCindi), a2DatAct, fmt=sFmtMI+sFmtDat, comments='', header=h)
    else:
        print '        ... file "{}" not available!'.format(sPthBlick)


def CINDI2SemiBlind_HCHO(par, sInstituteC, sDate, sLoc, sInstNum, sPan, sProcGas, sSCode, sFCode, sVers, nHead, nHeadL1, iCcCol):

    #> Get paths
    sPthCindi, sPthBlick, sPthL1Blick = GetPths(sDate, sPan, sLoc, sInstituteC, sInstNum, sProcGas, sVers, sSCode, sFCode)

    # > Read Blick data
    sPthDat = glob(ospath.join(par['sL2FitPth'], sPthBlick))
    if len(sPthDat) > 0:
        sPthDat = LatestFile(sPthDat)
        for sPthBlickGlob in sPthDat:
            a2Dat = loadtxt(sPthBlickGlob,
                            converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                            skiprows=nHead)
            # > Read VIS data for EH routine
            sPthDat = glob(sPthBlickGlob.partition('fC')[0] + 'fC07Oc*')
            sPthBlickGlob = LatestFile(sPthDat)[0]
            a2DatEH = loadtxt(sPthBlickGlob,
                              converters={0: lambda s: array([ord(c) ** 2 for c in s]).sum(), 1: datestr2num},
                              skiprows=103)
            bIsEH = a2DatEH[:, 0] == array([ord(c) ** 2 for c in 'EH']).sum()
            a2DatEH = a2DatEH[bIsEH, :a2Dat.shape[1]]  # extract only EH routines and cut by length of original data file
            if len(a2DatEH) > 0:
                # a2DatEH[:, 15:] = -999
                # > Merge UV with VIS data
                iMergPoint = where(a2Dat[:, 1] >= a2DatEH[-1, 1])[0][0]
                a2Dat = vstack((a2Dat[:iMergPoint, :], a2DatEH, a2Dat[iMergPoint:, :]))
                bIsEH = a2Dat[:, 0] == array([ord(c) ** 2 for c in 'EH']).sum()

        #> Get wavelength vector
        GD = GetDataOfRoutine(int(sPan.split('Pandora')[-1][:-2]), int(sPan.split('Pandora')[-1][-1]),
                              sLoc, par['sBlickRootPth'], par['sL0Pth'],
                              par['sOFPth'], par['sCFPth'], [sSCode, -1, -1], [-1, -1, -1])
        _, a1Wvl, _, panPars = GD.GetCalData(num2date(datestr2num(str(sDate))).replace(tzinfo=None))

        #> Read Blick L1 data
        a2Cc = None
        lPthBlickL1Glob = LatestFile(glob(ospath.join(par['sL1Pth'], sPthL1Blick)))
        for sPthBlickL1Glob in lPthBlickL1Glob:
            a2DatL1 = loadtxt(sPthBlickL1Glob,
                              converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                              skiprows=nHeadL1)
            bIdG = zeros((a2DatL1.shape[0]), dtype=bool)
            #> Reduce to L2Fit data
            for iRtnI in xrange(a2Dat.shape[0]):
                idx = where((a2DatL1[:, [3, 4]] == a2Dat[iRtnI, [3, 4]]).sum(axis=1) == 2)[0]
                if len(idx) > 0:
                    if len(idx) == 1:
                        bIdG[idx] = True
                    else:
                        assert False, 'routine doubled!'
                else:
                    assert False, 'corresponding routine not in L1 file!'
            assert (bIdG.shape[0] == a2DatL1.shape[0]) & (bIdG.sum() == a2Dat.shape[0]), 'L1 file filtering issue!'
            a2Cc = a2DatL1[bIdG, iCcCol-1:iCcCol-1+2048]
            a1Scl = a2DatL1[bIdG, iCcCol-1-3]
        fWvlRelInt = 340.
        a1RelInt = ones((a2Cc.shape[0]))
        fWvlColIdx0 = 340.
        fWvlColIdx1 = 359.
        a1ColIdx = ones((a2Cc.shape[0]))
        for iRtnI in xrange(a2Cc.shape[0]):
            #> Relative intensity
            bWvlIdx = (a1Wvl >=  fWvlRelInt - par['fWvlIntAvg']) & (a1Wvl <=  fWvlRelInt + par['fWvlIntAvg'])
            a1RelInt[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx] / a1Scl[iRtnI])
            #> Color index
            bWvlIdx0 = (a1Wvl >= fWvlColIdx0 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx0 + par['fWvlIntAvg'])
            bWvlIdx1 = (a1Wvl >= fWvlColIdx1 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx1 + par['fWvlIntAvg'])
            a1ColIdx[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx0] / a1Scl[iRtnI]) / mean(a2Cc[iRtnI, bWvlIdx1] / a1Scl[iRtnI])

        #> Write data
        ##> Meta info
        a2DatAct = array([]).reshape((a2Dat.shape[0], 0))
        for addCol in [3, 6, 8, 9, 12, 14]:
            if addCol == 3:
                ###> Decimal doy
                delta = (date(int(str(sDate)[:4]), 1, 1) - date(2000, 1, 1)).days
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] - delta))
                ###> Decimal hour of day
                a2DatAct = hstack((a2DatAct, a2Dat[:, [2]] % 1 * 24.))
            elif addCol == 12:
                ###> zenith viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VZA = a2Dat[:, [addCol - 1]]
                a1VZAMode = a2Dat[:, [addCol]]
                a1VZA[a1VZAMode == 0] = a1VZA[a1VZAMode == 0]
                a1VZA[a1VZAMode == 1] = a2Dat[(a1VZAMode == 1).squeeze(), 8 - 1] + a1VZA[a1VZAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVZA = -0.428
                    a1VZA[a1VZA > 5.] += fCorVZA
                ###> Convert to nomial elevation angles
                if par['doOvrwVZA']:
                    a1Elev = ConvertNominalVZA(par, a2Dat[:, 0], a1VZA, int(sPan.split('Pandora')[-1][-1]), doElev=True)
                else:
                    a1Elev = 90. - a1VZA
                a2DatAct = hstack((a2DatAct, a1Elev))
            elif addCol == 14:
                ###> azimuth viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VAA = a2Dat[:, [addCol - 1]]
                a1VAAMode = a2Dat[:, [addCol]]
                a1VAA[a1VZAMode == 0] = a1VAA[a1VAAMode == 0]
                a1VAA[a1VZAMode == 1] = a2Dat[(a1VAAMode == 1).squeeze(), 9 - 1] + a1VAA[a1VAAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVAA = 0.246
                    a1VAA[a1VZA > 5.] += fCorVAA
                a2DatAct = hstack((a2DatAct, a1VAA))
            else:
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))

        ##> Gases
        bOvrw = a2Dat[:, 15] > 2
        if bIsEH.sum() > 0:
            bOvrw |= bIsEH
        ###> HCHO DSC and rescale
        for addCol in [35, 36]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * gp.unitcf[1][1] * 1.e-15))
            a2DatAct[bOvrw, -1] = -999
        ###> O2O2 DSC and rescale
        for addCol in [32, 33]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * 1.e-40))
            a2DatAct[bOvrw, -1] = -999
        ###> NO2 DSC and rescale
        for addCol in [29, 30]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * gp.unitcf[1][1] * 1.e-15))
            a2DatAct[bOvrw, -1] = -999
        ###> O3 DSC with Temp fit to two NO2 DSCs and rescale
        xs1 = 24
        uxs1 = 25
        xst2 = 26
        uxst2 = 27
        TX = 225.
        T1 = 223.
        T2 = 243.
        ####> T = T1
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))
                           * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + ((1 - (a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/(1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-20))
        ####> T = T2
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1))
                           * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + (((a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/((a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-20))
        ###> BrO DSC and rescale
        for addCol in [38, 39]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * 1.e-15))
            a2DatAct[bOvrw, -1] = -999
        ###> Ring DSC and rescale
        for addCol in [43, 44]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999

        ##> Aux fit info
        for addCol in [18, 61]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999
        ##> Add relative intensity
        a1RelInt[a1RelInt == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1RelInt.reshape((a1RelInt.shape[0], 1))))
        ##> Add color index
        a1ColIdx[isinf(a1ColIdx)] = -999
        a1ColIdx[a1ColIdx == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1ColIdx.reshape((a1RelInt.shape[0], 1))))
        ##> Add intensity offset
        addCol = 57  # constant term
        a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
        a2DatAct[bOvrw, -1] = -999
        addCol = 59  # linear term
        a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
        a2DatAct[bOvrw, -1] = -999
        ##> Add O3 from retrieval output
        xs1 = 24
        uxs1 = 25
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [xs1 - 1]] * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [uxs1 - 1]] * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999

        # > Write header
        c = 0
        h = ''
        #> meta info
        h += '* NofHeaderlines: {}'.format(36) + '\n'
        h += '* NofColumns: {}'.format(a2DatAct.shape[1]) + '\n'
        h += '* Instrument identifier: {}_MAXDOAS'.format(sInstituteC) + '\n'
        h += '* Retrieval code: BlickP (v1.3.1, Feb 2017)' + '\n'
        h += '* Created by: Martin Tiefengraber' + '\n'
        h += '* Version: {}_v{}'.format(sProcGas, sVers) + '\n'; c += 1
        h += '* X-Axis (Col {}) = Day of year (DOY) 2016'.format(c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Time of day in hours (UTC)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Total Integration Time(s)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Zenith Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Azimuth Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Elevation Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Viewing Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        #> gases
        h += '* Y{}-Axis (Col {}) = HCHO_DSCD_297 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = HCHO_DSCD_297_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O4_DSCD_293 (1*10^40 molec2/cm5)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O4_DSCD_293_Error (1*10^40 molec2/cm5)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223 (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_243 (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_243_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = BrO_DSCD_223 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = BrO_DSCD_223_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        #> aux fit info
        h += '* Y{}-Axis (Col {}) = Ring'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Ring_Error'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Fit RMS (in OD)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Spectrum shift (nm, against FRS reference)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Relative Intensity (counts/integration time @ {}nm)'.format(c-1, c, fWvlRelInt) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Colour index: ({} nm / {} nm)'.format(c-1, c, fWvlColIdx0, fWvlColIdx1) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = intensity offset with normalisation by I, I is the mean intensity in the spectral analysis windows, constant term'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = intensity offset, linear term'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '*DOY UTC Tint SZA SAA Elev Viewing_angle HCHO_DSCD_297 HCHO_DSCD_297_Error O4_DSCD_293 O4_DSCD_293_Error NO2_DSCD_298 NO2_DSCD_298_Error O3_DSCD_223 O3_DSCD_223_Error O3_DSCD_243 O3_DSCD_243_Error BrO_DSCD_223 BrO_DSCD_223_Error Ring Ring_Error RMS Spectrum_shift Intens(340) CI(340/359) offset_cst offset_lin O3_DSCD O3_DSCD_Error '

        #> save file
        sFmtMI = ['%.5f', '%.5f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'] # meta info
        sFmtDat = ['%.4e' for _ in xrange(a2DatAct.shape[1] - len(sFmtMI))]
        savetxt(ospath.join(par['sCampPth'], sPthCindi), a2DatAct, fmt=sFmtMI+sFmtDat, comments='', header=h)
    else:
        print '        ... file "{}" not available!'.format(sPthBlick)


def CINDI2SemiBlind_O3vis(par, sInstituteC, sDate, sLoc, sInstNum, sPan, sProcGas, sSCode, sFCode, sVers, nHead, nHeadL1, iCcCol):

    #> Get paths
    sPthCindi, sPthBlick, sPthL1Blick = GetPths(sDate, sPan, sLoc, sInstituteC, sInstNum, sProcGas, sVers, sSCode, sFCode)

    # > Read Blick data
    sPthDat = glob(ospath.join(par['sL2FitPth'], sPthBlick))
    if len(sPthDat) > 0:
        sPthDat = LatestFile(sPthDat)
        for sPthBlickGlob in sPthDat:
            a2Dat = loadtxt(sPthBlickGlob,
                            converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                            skiprows=nHead)

        #> Get wavelength vector
        GD = GetDataOfRoutine(int(sPan.split('Pandora')[-1][:-2]), int(sPan.split('Pandora')[-1][-1]),
                              sLoc, par['sBlickRootPth'], par['sL0Pth'],
                              par['sOFPth'], par['sCFPth'], [sSCode, -1, -1], [-1, -1, -1])
        _, a1Wvl, _, panPars = GD.GetCalData(num2date(datestr2num(str(sDate))).replace(tzinfo=None))

        #> Read Blick L1 data
        a2Cc = None
        lPthBlickL1Glob = LatestFile(glob(ospath.join(par['sL1Pth'], sPthL1Blick)))
        for sPthBlickL1Glob in lPthBlickL1Glob:
            a2DatL1 = loadtxt(sPthBlickL1Glob,
                              converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                              skiprows=nHeadL1)
            bIdG = zeros((a2DatL1.shape[0]), dtype=bool)
            #> Reduce to L2Fit data
            for iRtnI in xrange(a2Dat.shape[0]):
                idx = where((a2DatL1[:, [3, 4]] == a2Dat[iRtnI, [3, 4]]).sum(axis=1) == 2)[0]
                if len(idx) > 0:
                    if len(idx) == 1:
                        bIdG[idx] = True
                    else:
                        assert True, 'routine doubled!'
                else:
                    assert True, 'corresponding routine not in L1 file!'
            assert (bIdG.shape[0] == a2DatL1.shape[0]) & (bIdG.sum() == a2Dat.shape[0]), 'L1 file filtering issue!'
            a2Cc = a2DatL1[bIdG, iCcCol-1:iCcCol-1+2048]
            a1Scl = a2DatL1[bIdG, iCcCol-1-3]
        fWvlRelInt = 500.
        a1RelInt = ones((a2Cc.shape[0]))
        fWvlColIdx0 = 440.
        fWvlColIdx1 = 500.
        a1ColIdx = ones((a2Cc.shape[0]))
        for iRtnI in xrange(a2Cc.shape[0]):
            #> Relative intensity
            bWvlIdx = (a1Wvl >=  fWvlRelInt - par['fWvlIntAvg']) & (a1Wvl <=  fWvlRelInt + par['fWvlIntAvg'])
            a1RelInt[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx] / a1Scl[iRtnI])
            #> Color index
            bWvlIdx0 = (a1Wvl >= fWvlColIdx0 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx0 + par['fWvlIntAvg'])
            bWvlIdx1 = (a1Wvl >= fWvlColIdx1 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx1 + par['fWvlIntAvg'])
            a1ColIdx[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx0] / a1Scl[iRtnI]) / mean(a2Cc[iRtnI, bWvlIdx1] / a1Scl[iRtnI])

        #> Write data
        ##> Meta info
        a2DatAct = array([]).reshape((a2Dat.shape[0], 0))
        for addCol in [3, 6, 8, 9, 12, 14]:
            if addCol == 3:
                ###> Decimal doy
                delta = (date(int(str(sDate)[:4]), 1, 1) - date(2000, 1, 1)).days
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] - delta))
                ###> Decimal hour of day
                a2DatAct = hstack((a2DatAct, a2Dat[:, [2]] % 1 * 24.))
            elif addCol == 12:
                ###> zenith viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VZA = a2Dat[:, [addCol - 1]]
                a1VZAMode = a2Dat[:, [addCol]]
                a1VZA[a1VZAMode == 0] = a1VZA[a1VZAMode == 0]
                a1VZA[a1VZAMode == 1] = a2Dat[(a1VZAMode == 1).squeeze(), 8 - 1] + a1VZA[a1VZAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVZA = -0.428
                    a1VZA[a1VZA > 5.] += fCorVZA
                ###> Convert to nomial elevation angles
                if par['doOvrwVZA']:
                    a1Elev = ConvertNominalVZA(par, a2Dat[:, 0], a1VZA, int(sPan.split('Pandora')[-1][-1]), doElev=True)
                else:
                    a1Elev = 90. - a1VZA
                a2DatAct = hstack((a2DatAct, a1Elev))
            elif addCol == 14:
                ###> azimuth viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VAA = a2Dat[:, [addCol - 1]]
                a1VAAMode = a2Dat[:, [addCol]]
                a1VAA[a1VZAMode == 0] = a1VAA[a1VAAMode == 0]
                a1VAA[a1VZAMode == 1] = a2Dat[(a1VAAMode == 1).squeeze(), 9 - 1] + a1VAA[a1VAAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVAA = 0.246
                    a1VAA[a1VZA > 5.] += fCorVAA
                a2DatAct = hstack((a2DatAct, a1VAA))
            else:
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))

        ##> Gases
        bOvrw = a2Dat[:, 15] > 2
        ###> O3 DSC with Temp fit to two NO2 DSCs and rescale
        ####> True
        xs1 = 24
        uxs1 = 25
        xst2 = 26
        uxst2 = 27
        TX = 225.
        T1 = 223.
        T2 = 293.
        ####> T = T1
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))
                           * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + ((1 - (a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/(1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-20))
        ####> T = T2
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1))
                           * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + (((a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/((a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-20))
        ###> O2O2 DSC and rescale
        for addCol in [34, 35]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * 1.e-40))
            a2DatAct[bOvrw, -1] = -999
        ###> NO2 DSC with Temp fit to two NO2 DSCs and rescale
        ####> True
        xs1 = 29
        uxs1 = 30
        xst2 = 31
        uxst2 = 32
        TX = 254.5
        T1 = 220.
        T2 = 298.
        ####> T = T2
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1))
                           * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + (((a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/((a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-15))
        ####> T = T1
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))
                           * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + ((1 - (a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/(1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-15))
        ###> H2O DSC and rescale
        for addCol in [37, 38]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * 3.3428e22 * 1.e-23))
            a2DatAct[bOvrw, -1] = -999
        ###> Ring DSC and rescale
        for addCol in [42, 43]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999

        ##> Aux fit info
        for addCol in [18, 60]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999
        ##> Add relative intensity
        a1RelInt[a1RelInt == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1RelInt.reshape((a1RelInt.shape[0], 1))))
        ##> Add color index
        a1ColIdx[isinf(a1ColIdx)] = -999
        a1ColIdx[a1ColIdx == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1ColIdx.reshape((a1RelInt.shape[0], 1))))
        ##> Add intensity offset
        addCol = 56  # constant term
        a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
        a2DatAct[bOvrw, -1] = -999
        addCol = 58  # linear term
        a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
        a2DatAct[bOvrw, -1] = -999
        ##> O3 from retrieval output
        xs1 = 24
        uxs1 = 25
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [xs1 - 1]] * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [uxs1 - 1]] * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        ##> NO2 from retrieval output
        xs1 = 29
        uxs1 = 30
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [xs1 - 1]] * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [uxs1 - 1]] * gp.unitcf[1][1] * 1.e-15))
        a2DatAct[bOvrw, -1] = -999

        # > Write header
        c = 0
        h = ''
        #> meta info
        h += '* NofHeaderlines: {}'.format(38) + '\n'
        h += '* NofColumns: {}'.format(a2DatAct.shape[1]) + '\n'
        h += '* Instrument identifier: {}_MAXDOAS'.format(sInstituteC) + '\n'
        h += '* Retrieval code: BlickP (v1.3.1, Feb 2017)' + '\n'
        h += '* Created by: Martin Tiefengraber' + '\n'
        h += '* Version: {}_v{}'.format(sProcGas, sVers) + '\n'; c += 1
        h += '* X-Axis (Col {}) = Day of year (DOY) 2016'.format(c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Time of day in hours (UTC)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Total Integration Time(s)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Zenith Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Azimuth Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Elevation Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Viewing Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        #> gases
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223 (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_293 (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_293_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O4_DSCD_293 (1*10^40 molec2/cm5)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O4_DSCD_293_Error (1*10^40 molec2/cm5)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_220 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_220_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = H2O_DSCD_296 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = H2O_DSCD_296_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        #> aux fit info
        h += '* Y{}-Axis (Col {}) = Ring'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Ring_Error'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Fit RMS (in OD)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Spectrum shift (nm, against FRS reference)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Relative Intensity (counts/integration time @ {}nm)'.format(c-1, c, fWvlRelInt) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Colour index: ({} nm / {} nm)'.format(c-1, c, fWvlColIdx0, fWvlColIdx1) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = intensity offset with normalisation by I, I is the mean intensity in the spectral analysis windows, constant term'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = intensity offset, linear term'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '*DOY UTC Tint SZA SAA Elev Viewing_angle O3_DSCD_223 O3_DSCD_223_Error O3_DSCD_293 O3_DSCD_293_Error O4_DSCD_293 O4_DSCD_293_Error NO2_DSCD_298 NO2_DSCD_298_Error NO2_DSCD_220 NO2_DSCD_220_Error H2O_DSCD_296 H2O_DSCD_296_Error Ring Ring_Error RMS Spectrum_shift Intens(500) CI(440/500) offset_cst offset_lin O3_DSCD O3_DSCD_Error NO2_DSCD NO2_DSCD_Error'

        #> save file
        sFmtMI = ['%.5f', '%.5f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'] # meta info
        sFmtDat = ['%.4e' for _ in xrange(a2DatAct.shape[1] - len(sFmtMI))]
        savetxt(ospath.join(par['sCampPth'], sPthCindi), a2DatAct, fmt=sFmtMI+sFmtDat, comments='', header=h)
    else:
        print '        ... file "{}" not available!'.format(sPthBlick)


def CINDI2SemiBlind_O3uv(par, sInstituteC, sDate, sLoc, sInstNum, sPan, sProcGas, sSCode, sFCode, sVers, nHead, nHeadL1, iCcCol):

    #> Get paths
    sPthCindi, sPthBlick, sPthL1Blick = GetPths(sDate, sPan, sLoc, sInstituteC, sInstNum, sProcGas, sVers, sSCode, sFCode)

    # > Read Blick data
    sPthDat = glob(ospath.join(par['sL2FitPth'], sPthBlick))
    if len(sPthDat) > 0:
        sPthDat = LatestFile(sPthDat)
        for sPthBlickGlob in sPthDat:
            a2Dat = loadtxt(sPthBlickGlob,
                            converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                            skiprows=nHead)
            # > Read VIS data for EH routine
            sPthDat = glob(sPthBlickGlob.partition('fC')[0] + 'fC07Oc*')
            sPthBlickGlob = LatestFile(sPthDat)[0]
            a2DatEH = loadtxt(sPthBlickGlob,
                              converters={0: lambda s: array([ord(c) ** 2 for c in s]).sum(), 1: datestr2num},
                              skiprows=103)
            bIsEH = a2DatEH[:, 0] == array([ord(c) ** 2 for c in 'EH']).sum()
            a2DatEH = a2DatEH[bIsEH, :a2Dat.shape[1]]  # extract only EH routines and cut by length of original data file
            if len(a2DatEH) > 0:
                # a2DatEH[:, 15:] = -999
                # > Merge UV with VIS data
                iMergPoint = where(a2Dat[:, 1] >= a2DatEH[-1, 1])[0][0]
                a2Dat = vstack((a2Dat[:iMergPoint, :], a2DatEH, a2Dat[iMergPoint:, :]))
                bIsEH = a2Dat[:, 0] == array([ord(c) ** 2 for c in 'EH']).sum()

        #> Get wavelength vector
        GD = GetDataOfRoutine(int(sPan.split('Pandora')[-1][:-2]), int(sPan.split('Pandora')[-1][-1]),
                              sLoc, par['sBlickRootPth'], par['sL0Pth'],
                              par['sOFPth'], par['sCFPth'], [sSCode, -1, -1], [-1, -1, -1])
        _, a1Wvl, _, panPars = GD.GetCalData(num2date(datestr2num(str(sDate))).replace(tzinfo=None))

        #> Read Blick L1 data
        a2Cc = None
        lPthBlickL1Glob = LatestFile(glob(ospath.join(par['sL1Pth'], sPthL1Blick)))
        for sPthBlickL1Glob in lPthBlickL1Glob:
            a2DatL1 = loadtxt(sPthBlickL1Glob,
                              converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                              skiprows=nHeadL1)
            bIdG = zeros((a2DatL1.shape[0]), dtype=bool)
            #> Reduce to L2Fit data
            for iRtnI in xrange(a2Dat.shape[0]):
                idx = where((a2DatL1[:, [3, 4]] == a2Dat[iRtnI, [3, 4]]).sum(axis=1) == 2)[0]
                if len(idx) > 0:
                    if len(idx) == 1:
                        bIdG[idx] = True
                    else:
                        assert True, 'routine doubled!'
                else:
                    assert True, 'corresponding routine not in L1 file!'
            assert (bIdG.shape[0] == a2DatL1.shape[0]) & (bIdG.sum() == a2Dat.shape[0]), 'L1 file filtering issue!'
            a2Cc = a2DatL1[bIdG, iCcCol-1:iCcCol-1+2048]
            a1Scl = a2DatL1[bIdG, iCcCol-1-3]
        fWvlRelInt = 340.
        a1RelInt = ones((a2Cc.shape[0]))
        fWvlColIdx0 = 320.
        fWvlColIdx1 = 340.
        a1ColIdx = ones((a2Cc.shape[0]))
        for iRtnI in xrange(a2Cc.shape[0]):
            #> Relative intensity
            bWvlIdx = (a1Wvl >=  fWvlRelInt - par['fWvlIntAvg']) & (a1Wvl <=  fWvlRelInt + par['fWvlIntAvg'])
            a1RelInt[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx] / a1Scl[iRtnI])
            #> Color index
            bWvlIdx0 = (a1Wvl >= fWvlColIdx0 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx0 + par['fWvlIntAvg'])
            bWvlIdx1 = (a1Wvl >= fWvlColIdx1 - par['fWvlIntAvg']) & (a1Wvl <= fWvlColIdx1 + par['fWvlIntAvg'])
            a1ColIdx[iRtnI] = mean(a2Cc[iRtnI, bWvlIdx0] / a1Scl[iRtnI]) / mean(a2Cc[iRtnI, bWvlIdx1] / a1Scl[iRtnI])

        #> Write data
        ##> Meta info
        a2DatAct = array([]).reshape((a2Dat.shape[0], 0))
        for addCol in [3, 6, 8, 9, 12, 14]:
            if addCol == 3:
                ###> Decimal doy
                delta = (date(int(str(sDate)[:4]), 1, 1) - date(2000, 1, 1)).days
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] - delta))
                ###> Decimal hour of day
                a2DatAct = hstack((a2DatAct, a2Dat[:, [2]] % 1 * 24.))
            elif addCol == 12:
                ###> zenith viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VZA = a2Dat[:, [addCol - 1]]
                a1VZAMode = a2Dat[:, [addCol]]
                a1VZA[a1VZAMode == 0] = a1VZA[a1VZAMode == 0]
                a1VZA[a1VZAMode == 1] = a2Dat[(a1VZAMode == 1).squeeze(), 8 - 1] + a1VZA[a1VZAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVZA = -0.428
                    a1VZA[a1VZA > 5.] += fCorVZA
                ###> Convert to nomial elevation angles
                if par['doOvrwVZA']:
                    a1Elev = ConvertNominalVZA(par, a2Dat[:, 0], a1VZA, int(sPan.split('Pandora')[-1][-1]), doElev=True)
                else:
                    a1Elev = 90. - a1VZA
                a2DatAct = hstack((a2DatAct, a1Elev))
            elif addCol == 14:
                ###> azimuth viewing angle to elevation angle
                ####> discriminate between relative and absolute viewing angles
                a1VAA = a2Dat[:, [addCol - 1]]
                a1VAAMode = a2Dat[:, [addCol]]
                a1VAA[a1VZAMode == 0] = a1VAA[a1VAAMode == 0]
                a1VAA[a1VZAMode == 1] = a2Dat[(a1VAAMode == 1).squeeze(), 9 - 1] + a1VAA[a1VAAMode == 1]
                ###> Correct Pan128s2 before 20160915, correct all angles except zenith measurement
                if (sPan == 'Pandora128s2') and (int(sDate) <= 20160915):
                    print 'VA correction applied'
                    fCorVAA = 0.246
                    a1VAA[a1VZA > 5.] += fCorVAA
                a2DatAct = hstack((a2DatAct, a1VAA))
            else:
                a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))

        ##> Gases
        bOvrw = a2Dat[:, 15] > 2
        if bIsEH.sum() > 0:
            bOvrw |= bIsEH
        ###> O3 DSC with Temp fit to two NO2 DSCs and rescale
        ####> True
        xs1 = 24
        uxs1 = 25
        xst2 = 26
        uxst2 = 27
        # TX = 293.
        TX = 225.
        T1 = 223.
        T2 = 293.
        ####> T = T1
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))
                           * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + ((1 - (a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/(1 - (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-20))
        ####> T = T2
        a2DatAct = hstack((a2DatAct,
                           (a2Dat[:, [xs1 - 1]] * (a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1))
                           * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct, ones((a2DatAct.shape[0], 1)) * -999))
        # a2DatAct = hstack((a2DatAct,
        #                    ((a2Dat[:, [uxs1 - 1]]/a2Dat[:, [xs1 - 1]])**2
        #                     + (((a2Dat[:, [uxst2 - 1]] - TX) / (T2 - T1))/((a2Dat[:, [xst2 - 1]] - TX) / (T2 - T1)))**2)**0.5
        #                    * gp.unitcf[1][1] * 1.e-20))
        ###> NO2 DSC and rescale
        for addCol in [29, 30]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * gp.unitcf[1][1] * 1.e-15))
            a2DatAct[bOvrw, -1] = -999
        ###> HCHO DSC and rescale
        for addCol in [32, 33]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]] * gp.unitcf[1][1] * 1.e-15))
            a2DatAct[bOvrw, -1] = -999
        ###> Ring DSC and rescale
        for addCol in [37, 38]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999

        ##> Aux fit info
        for addCol in [18, 51]:
            a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
            a2DatAct[bOvrw, -1] = -999
        ##> Add relative intensity
        a1RelInt[a1RelInt == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1RelInt.reshape((a1RelInt.shape[0], 1))))
        ##> Add color index
        a1ColIdx[isinf(a1ColIdx)] = -999
        a1ColIdx[a1ColIdx == 0.] = -999
        a2DatAct = hstack((a2DatAct, a1ColIdx.reshape((a1RelInt.shape[0], 1))))
        ##> Add intensity offset
        addCol = 47  # constant term
        a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
        a2DatAct[bOvrw, -1] = -999
        addCol = 49  # linear term
        a2DatAct = hstack((a2DatAct, a2Dat[:, [addCol - 1]]))
        a2DatAct[bOvrw, -1] = -999
        ##> O3 from retrieval output
        xs1 = 24
        uxs1 = 25
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [xs1 - 1]] * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999
        a2DatAct = hstack((a2DatAct,
                           a2Dat[:, [uxs1 - 1]] * gp.unitcf[1][1] * 1.e-20))
        a2DatAct[bOvrw, -1] = -999

        # > Write header
        c = 0
        h = ''
        #> meta info
        h += '* NofHeaderlines: {}'.format(32) + '\n'
        h += '* NofColumns: {}'.format(a2DatAct.shape[1]) + '\n'
        h += '* Instrument identifier: {}_MAXDOAS'.format(sInstituteC) + '\n'
        h += '* Retrieval code: BlickP (v1.3.1, Feb 2017)' + '\n'
        h += '* Created by: Martin Tiefengraber' + '\n'
        h += '* Version: {}_v{}'.format(sProcGas, sVers) + '\n'; c += 1
        h += '* X-Axis (Col {}) = Day of year (DOY) 2016'.format(c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Time of day in hours (UTC)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Total Integration Time(s)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Zenith Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Solar Azimuth Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Elevation Angle (deg)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Viewing Angle (deg) North=0, East=90'.format(c-1, c) + '\n'; c += 1
        #> gases
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223 (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_223_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_293 (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_293_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = NO2_DSCD_298_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = HCHO_DSCD_297 (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = HCHO_DSCD_297_Error (1*10^15 molec/cm2)'.format(c-1, c) + '\n'; c += 1
        #> aux fit info
        h += '* Y{}-Axis (Col {}) = Ring'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Ring_Error'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Fit RMS (in OD)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Spectrum shift (nm, against FRS reference)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Relative Intensity (counts/integration time @ {}nm)'.format(c-1, c, fWvlRelInt) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = Colour index: ({} nm / {} nm)'.format(c-1, c, fWvlColIdx0, fWvlColIdx1) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = intensity offset with normalisation by I, I is the mean intensity in the spectral analysis windows, constant term'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = intensity offset, linear term'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '* Y{}-Axis (Col {}) = O3_DSCD_Error (1*10^20 molecules/cm2)'.format(c-1, c) + '\n'; c += 1
        h += '*DOY UTC Tint SZA SAA Elev Viewing_angle O3_DSCD_223 O3_DSCD_223_Error O3_DSCD_293 O3_DSCD_293_Error NO2_DSCD_298 NO2_DSCD_298_Error HCHO_DSCD_297 HCHO_DSCD_297_Error Ring Ring_Error RMS Spectrum_shift Intens(340) CI(320/340) offset_cst offset_lin O3_DSCD O3_DSCD_Error'

        if par['dProdAna'][sProcGas]:
            f, ax = subplots(10, 2, figsize=(3, 8))
            idxNoon = argmin(a2DatAct[:, 3])
            fltAm = a2DatAct[:idxNoon, 5] > 85.
            fltPm = a2DatAct[idxNoon:, 5] > 85.
            f.suptitle(sPan + ', ' + sProcGas + ', ' + str(sDate))
            ax[0, 0].set_title('AM')
            ax[0, 1].set_title('PM')
            ylab = 'O3 DSCD'
            ylim = ylimO3
            scT1 = 7
            scT2 = 9
            sc = 23
            usc = 24
            row = 0
            ax[row, 0].plot(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, scT1][fltAm], 'o', label='223')
            ax[row, 0].plot(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, scT2][fltAm], 'o', label='293')
            ax[row, 0].errorbar(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], yerr=a2DatAct[:idxNoon, usc][fltAm], fmt='o', label='fit', capsize=0)
            ax[row, 1].plot(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, scT1][fltPm], 'o', label='298')
            ax[row, 1].plot(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, scT2][fltPm], 'o', label='220')
            ax[row, 1].errorbar(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], yerr=a2DatAct[idxNoon:, usc][fltPm], fmt='o', label='fit', capsize=0)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ax[row, 0].set_ylabel(ylab)
            ylab = 'NO2 DSCD'
            ylim = ylimNO2
            sc = 11
            usc = 12
            row = 1
            ax[row, 0].errorbar(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], yerr=a2DatAct[:idxNoon, usc][fltAm], fmt='o', capsize=0)
            ax[row, 1].errorbar(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], yerr=a2DatAct[idxNoon:, usc][fltPm], fmt='o', capsize=0)
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'HCHO DSCD'
            ylim = ylimHCHO
            sc = 13
            usc = 14
            row = 2
            ax[row, 0].errorbar(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], yerr=a2DatAct[:idxNoon, usc][fltAm], fmt='o', capsize=0)
            ax[row, 1].errorbar(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], yerr=a2DatAct[idxNoon:, usc][fltPm], fmt='o', capsize=0)
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'Ring DSCD'
            ylim = ylimRing
            sc = 15
            usc = 16
            row = 4
            ax[row, 0].errorbar(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], yerr=a2DatAct[:idxNoon, usc][fltAm], fmt='o', capsize=0)
            ax[row, 1].errorbar(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], yerr=a2DatAct[idxNoon:, usc][fltPm], fmt='o', capsize=0)
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'RMS'
            ylim = ylimRMS
            sc = 17
            row = 5
            ax[row, 0].semilogy(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], 'o')
            ax[row, 1].semilogy(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], 'o')
            # [ax[row, i].yaxis.set_scale('log') for i in xrange(ax.shape[1])]
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'WVL shift'
            ylim = ylimWVL0
            sc = 18
            row = 6
            ax[row, 0].plot(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], 'o')
            ax[row, 1].plot(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], 'o')
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'Offset'
            ylim = ylimOffs
            sc = 21
            row = 7
            ax[row, 0].plot(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], 'o')
            ax[row, 1].plot(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], 'o')
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'Offset lin'
            ylim = ylimOffsLin
            sc = 22
            row = 7
            ax[row, 0].plot(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], 'o')
            ax[row, 1].plot(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], 'o')
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'RI'
            ylim = ylimRI
            sc = 19
            row = 8
            ax[row, 0].semilogy(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], 'o')
            ax[row, 1].semilogy(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], 'o')
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]
            ylab = 'CI'
            ylim = ylimCI
            sc = 20
            row = 9
            ax[row, 0].plot(a2DatAct[:idxNoon, 3][fltAm], a2DatAct[:idxNoon, sc][fltAm], 'o')
            ax[row, 1].plot(a2DatAct[idxNoon:, 3][fltPm], a2DatAct[idxNoon:, sc][fltPm], 'o')
            ax[row, 0].set_ylabel(ylab)
            [ax[row, i].set_ylim(ylim) for i in xrange(ax.shape[1])]

            [ax[i, 0].set_xlim([100, 40]) for i in xrange(ax.shape[0])]
            [ax[i, 1].set_xlim([40, 100]) for i in xrange(ax.shape[0])]

            [ax[i, 1].get_yaxis().set_visible(False) for i in xrange(ax.shape[0])]
            [ax[i, j].get_xaxis().set_visible(False) for i in xrange(ax.shape[0]-1) for j in xrange(ax.shape[1])]
            # f.subplots_adjust(vspace=0)

            f.savefig(sPan + sProcGas + '_' + str(sDate) + '.png', dpi=300., bbox_inches='tight')
            close()

        #> save file
        sFmtMI = ['%.5f', '%.5f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'] # meta info
        sFmtDat = ['%.4e' for _ in xrange(a2DatAct.shape[1] - len(sFmtMI))]
        savetxt(ospath.join(par['sCampPth'], sPthCindi), a2DatAct, fmt=sFmtMI+sFmtDat, comments='', header=h)
    else:
        print '        ... file "{}" not available!'.format(sPthBlick)
