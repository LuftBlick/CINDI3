__author__ = 'Martin Tiefengraber'
import sys
import os
# Define the desired path
path = os.path.normpath("C:/Blick/src/")
os.chdir(path)
from os import path as ospath
from numpy import array, argmin, hstack, loadtxt, savetxt, where, ones, polyval, unique, asarray, in1d, nan, isnan, \
    nanmean, sum, nansum
from datetime import datetime, timedelta
from matplotlib.dates import num2date, date2num
from params.InputParams import LoadParams
from GetDataOfRoutineBlick import GetDataOfRoutine
from ProcessingWrapper import ProcessingWrapper
from blick_countconverter import regrid_f

from blick_io import blick_io
io = blick_io()

from CINDI3Formats.CINDI3SemiBlindDirSun import \
    CINDI3SemiBlind_NO2vis, \
    CINDI3SemiBlind_NO2visSmall, \
    CINDI3SemiBlind_NO2uv, \
    CINDI3SemiBlind_HCHO, \
    CINDI3SemiBlind_O3vis, \
    CINDI3SemiBlind_O3uv

def ConvertToCINDI3BlindFmt(par):
    print('Start converting to CINDI-3 blind format ...')
    for sDate in par['iDate']:
        # > Loop Pandoras
        for sPanC in par['dPan']:
            # > Loop spectrometers
            for iSpec in par['dPan'][sPanC]:
                sPanNme = 'Pandora{}s{}'.format(sPanC, iSpec)
                sInstituteC, sInstNum = par['dPanIdInst'][sPanNme[7:]]
                for sProcGas in par['dProdFCode'].keys():
                    if iSpec == 109:
                        iHeadL1 = 8
                    elif iSpec == 2:
                        iHeadL1 = 109
                    iHeadL1 = par['dL1NHead'][sPanC+'s'+str(iSpec)][0]
                    iCcCol = par['dL1CcCol'][sPanC+'s'+str(iSpec)][0]
                    iHeadL2Fit = par['dProdNHead'][sProcGas][0]
                    sSCode = par['dProcCode']['s'][0]
                    sProc = par['dProdFCode'][sProcGas][0]
                    iProcVers = par['dProdVers'][sProcGas][0]
                    sProdRefT = par['dProdRefT'][sProcGas][0]
                    sProcType = par['dProcType']

                    print('    ... for Pandora {}, s{}, Prod {}({}), Date {}')\
                        .format(sPanC, iSpec, sProcGas, sProc, sDate)

                    eval(
                        'CINDI3SemiBlind_{}(par, sInstituteC, sDate, par["sLoc"], sInstNum, sPanNme, sProcGas,'
                         'sSCode, sProc, iProcVers, iHeadL2Fit, iHeadL1, iCcCol, sProcType, sProdRefT, sPanC, curr_ref)'.format(sProcGas)
                    )

    print('... finished converting to CINDI-3 blind format.')


def GetExternalReferenceFileName(par, sPanC, iSpec, sRtn, sFuFi, iDate, sDateRef, iSCode):

    sRefNme = ospath.join(par['sRefPth'], 'RefPan{}s{}Rtn{}FuFi{}from{}pm{}min-for{}_s{}.txt'
                          .format(sPanC, iSpec, sRtn, sFuFi, sDateRef,
                                  par['iRefAvgInt'], iDate, iSCode))

    return sRefNme


def PickDataWithTimeRange(par, lVal, lValCmb):

    lValRef = []
    lValRefCmb = []
    lDtUsed = []
    for iDate in xrange(len(lVal)):
        lValRefDate = []
        lValRefCmbDate = []
        #> Generate date time vector
        ##> Datetime and remove time zone
        a1Dt = array([num2date((lVal[iDate][i][0][1][0] + date2num(par['epoch'])*24.*3600.)/(24.*3600.)).replace(tzinfo=None)
                      for i in xrange(len(lVal[iDate]))])
        #> Get datetime range
        ##> Reference date
        a1DtRef = array([datetime.strptime(par['sRefDateTime'][i], '%Y%m%dT%H%M%SZ')
                         for i in xrange(len(par['sRefDateTime']))])
        if a1DtRef.shape[0] == 1:  # use the same reference for all days
            iDate = 0
        #>  Get averaging interval
        ##> Check maximum allowed distance to reference time
        a1DtRefBest = a1Dt[argmin(abs(a1Dt - a1DtRef[iDate]))]
        assert (a1DtRefBest >= a1DtRef[iDate] - timedelta(minutes=int(par['iRefOffs']))) & \
               (a1DtRefBest <= a1DtRef[iDate] + timedelta(minutes=int(par['iRefOffs']))), \
            'Wanted reference time (from {}) not within maximum boundaries (p/m {} min.)!'. \
            format(par['sRefDateTime'][iDate], int(par['iRefOffs']))
        ##> Get averagint interval
        bId = (a1Dt >= a1DtRefBest - timedelta(minutes=int(par['iRefAvgInt']))) & \
              (a1Dt <= a1DtRefBest + timedelta(minutes=int(par['iRefAvgInt'])))
        ##> Check if averaging interval contains data
        assert any(bId), \
            'Wanted closest reference time (from {}, p/m {} min.) not included in the data!'. \
            format(a1DtRefBest, int(par['iRefAvgInt']))
        #> Extract data for reference and average
        for iRtn in where(bId)[0]:
            lValRefDate.append(lVal[iDate][iRtn])
            lValRefCmbDate.append(lValCmb[iDate][iRtn])
        lValRef.append(lValRefDate)
        lValRefCmb.append(lValRefCmbDate)
        lDtUsed.append(a1Dt[bId])

    return lValRef, lValRefCmb, lDtUsed


def ProcessExternalReference(par):

    print('Start processing external reference ...')
    dFuFiAll = {}
    #> Loop Pandoras
    for sPanC in par['dPan']:
        dFuFiAll[sPanC] = {}
        #> Loop spectrometers
        for iSpec in par['dPan'][sPanC]:
            dFuFiAll[sPanC][iSpec] = {}
            #> Loop s-numbers
            for sSCode, iQSCode in zip(par['dSCode']['s'], par['dSCode']['qs']):
                dFuFiAll[sPanC][iSpec][sSCode] = {}
                GD = GetDataOfRoutine(int(sPanC), iSpec, par['sLoc'], par['sBlickRootPth'], par['sL0Pth'],
                                      par['sOFPth'], par['sCFPth'], par['CfSuffixRef'], [sSCode, -1, -1], [iQSCode, -1, -1])
                #> Reformat date vector
                a1Date = array([datetime.strptime(str(par['iDate'][i]), '%Y%m%d')
                                for i in xrange(len(par['iDate']))])
                for sRtn, iRtnCnt in zip(par['sRefRtn'], par['iRtnCnt']):
                    dFuFiAll[sPanC][iSpec][sSCode][sRtn] = {}
                    lDataCmb, lData, _ = GD.GetDataOfRoutine(sRtn, a1Date, iLev='0')
                    #> Delete empty dates (empty when routine not has not been measured)
                    lGoodDates = []
                    for iDateI in xrange(len(lData)):
                        if len(lData[iDateI]) > 0:
                            lGoodDates.append(iDateI)
                        else:
                            print('        ... {} skipped!').format(par['iDate'][iDateI])
                    lData = [lData[i] for i in lGoodDates]
                    lDataCmb = [lDataCmb[i] for i in lGoodDates]
                    a1Date = a1Date[lGoodDates]
                    par['iDate'] = par['iDate'][lGoodDates]
                    par['sRefDateTime'] = par['sRefDateTime'][lGoodDates]
                    if len(lData) > 0:
                        #> Find measurements within wanted temporal range
                        lData, lDataCmb, lDtUsed = PickDataWithTimeRange(par, lData, lDataCmb)
                        #> Loop routines
                        lWvl, lCc, lECc, lCcInfo = GD.DoDataCorrectionOneDate(2048, iRtnCnt, lDataCmb, lData, a1Date)
                        #> Average data
                        for iDateI in xrange(len(a1Date)):
                            #> Load calibration file data
                            calData, _, _, panPars = GD.GetCalData(a1Date[iDateI])
                            # bPixReg = calData[-1][3][0]
                            #> Get functional fiter positions used in measurements
                            lFuFi = []
                            for iRtnI in xrange(len(lDataCmb[iDateI])):
                                for iMeasI in xrange(len(lDataCmb[iDateI][iRtnI][2])):
                                    lFuFi.append(lDataCmb[iDateI][iRtnI][2][iMeasI][3])
                            a1FuFiMea = unique(array(lFuFi))
                            ##> Find name of functional filters used in measurements
                            a1FuFiPosNme = asarray([c.split('-') for c in panPars[1][6]])
                            a1FuFiPos = a1FuFiPosNme[:, 0].astype(int)
                            a1FuFiNme = a1FuFiPosNme[:, 1]
                            a1FuFiNme = a1FuFiNme[in1d(a1FuFiPos, a1FuFiMea)]
                            ##> Remove OPAQUE
                            a1FuFiMea = a1FuFiMea[~(a1FuFiNme=='OPAQUE')]
                            a1FuFiNme = a1FuFiNme[~(a1FuFiNme=='OPAQUE')]
                            ##> Fuse
                            dFuFi = {}
                            for iFuFi, sFuFi in zip(a1FuFiMea, a1FuFiNme):
                                dFuFi[iFuFi] = sFuFi
                            #> Nominel wavelength vector
                            a1Wvl = lWvl[iDateI]
                            #> Initial corrected wavelength vector
                            dWvlCor = {}
                            dCc = {}
                            dECc = {}
                            for sFuFi in a1FuFiNme:
                                dWvlCor[sFuFi] = ones((a1Wvl.shape[0], lCc[iDateI].shape[2],
                                                       lCc[iDateI].shape[0]/len(a1FuFiMea))) * \
                                                       a1Wvl.reshape((a1Wvl.shape[0], 1, 1))
                                dCc[sFuFi] = ones(dWvlCor[sFuFi].shape) * nan
                                dECc[sFuFi] = ones(dWvlCor[sFuFi].shape) * nan
                            #> Loop all routines and measurements within routines
                            for iRtnI in xrange(lCc[iDateI].shape[2]):
                                dFuFiCnt = {}
                                for sFuFi in a1FuFiNme:
                                    dFuFiCnt[sFuFi] = -1
                                lToDel = []
                                for iMeasI in xrange(lCc[iDateI].shape[0]):
                                    #> Only use measuremnts with more than 1 cycle and no saturated data
                                    if (lData[iDateI][iRtnI][iMeasI][4][0] > 1) and (lData[iDateI][iRtnI][iMeasI][4][2] >= 0):
                                        #> Indentify functional filter
                                        sFuFi = dFuFi[lDataCmb[iDateI][iRtnI][2][iMeasI * 2][3]]
                                        dFuFiCnt[sFuFi] += 1
                                        iMeasFuFiI = dFuFiCnt[sFuFi]
                                        #> Calculate corrected wavelength vector
                                        wlcpol = lCcInfo[iDateI][iRtnI][8][iMeasI][1][0]
                                        if len(wlcpol) > 0:  # wavelength shift could be retrieved
                                            a1DeltaWvl = polyval(wlcpol[::-1], a1Wvl - calData[8][0][0])
                                            dWvlCor[sFuFi][:, iRtnI, iMeasFuFiI] = dWvlCor[sFuFi][:, iRtnI, iMeasFuFiI] + a1DeltaWvl
                                        else:
                                            lToDel.append(iRtnI)
                                        #> Resave data for functional filter
                                        dCc[sFuFi][:, iRtnI, iMeasFuFiI] = lCc[iDateI][iMeasI, :, iRtnI]
                                        dECc[sFuFi][:, iRtnI, iMeasFuFiI] = lECc[iDateI][iMeasI, :, iRtnI]
                                    else:
                                        print('         ... measurement has to be skipped!')
                            #> Loop functional filter and create average wavelength vector and spectrum, regrid and save
                            for sFuFi in a1FuFiNme:

                                print('    ... for Pandora {}, s{}, Routine {}, Date {}, FuncFilt {}')\
                                    .format(sPanC, iSpec, sRtn, par['iDate'][iDateI], sFuFi)

                                #> Make reference median wavelength vector and spectrum
                                print(nanmean(dECc[sFuFi][:, :, 0], axis=(0, 1)))
                                ##> Intensity variation filter
                                meanDev = (nanmean(dCc[sFuFi][:, :, 0], axis=0) - nanmean(dCc[sFuFi][:, :, 0], axis=(0, 1))) / \
                                          nanmean(dCc[sFuFi][:, :, 0], axis=0)
                                bUsed = abs(100. * meanDev) < 15.  # in percent
                                a1WvlMean = nanmean(dWvlCor[sFuFi][:, bUsed, :], axis=(1, 2))
                                a1CcMean = nanmean(dCc[sFuFi][:, bUsed, :], axis=(1, 2))
                                a1ECcMean = nansum(dECc[sFuFi][:, bUsed, :]**2, axis=(1, 2))**0.5 / sum(bUsed)
                                # a1ECcMean = nanmean(dECc[sFuFi][:, bUsed, :], axis=(1, 2))
                                for iRtnI in xrange(dWvlCor[sFuFi].shape[1]):
                                    for iMeasI in xrange(dWvlCor[sFuFi].shape[2]):
                                        #> Regrid data to reference median
                                        if not any(isnan(dCc[sFuFi][:, iRtnI, iMeasI])):
                                            dCc[sFuFi][:, iRtnI, iMeasI] = regrid_f(dWvlCor[sFuFi][:, iRtnI, iMeasI],
                                                                                    dCc[sFuFi][:, iRtnI, iMeasI],
                                                                                    a1WvlMean, a1CcMean)
                                            dECc[sFuFi][:, iRtnI, iMeasI] = regrid_f(dWvlCor[sFuFi][:, iRtnI, iMeasI],
                                                                                     dECc[sFuFi][:, iRtnI, iMeasI],
                                                                                     a1WvlMean, a1ECcMean)
                                dCc[sFuFi] = nanmean(dCc[sFuFi], axis=(1, 2))
                                dECc[sFuFi] = nanmean(dECc[sFuFi], axis=(1, 2))
                                #> Save reference to file
                                ##> Make path
                                sRefNme = GetExternalReferenceFileName(par, sPanC, iSpec, sRtn, sFuFi, par['iDate'][iDateI],
                                                                       par['sRefDateTime'][iDateI], sSCode)
                                ##> Save
                                savetxt(sRefNme, hstack((a1WvlMean.reshape(2048,1),
                                                         dCc[sFuFi].reshape(2048,1),
                                                         dECc[sFuFi].reshape(2048,1))))
                            dFuFiAll[sPanC][iSpec][sSCode][sRtn][par['iDate'][iDateI]] = a1FuFiNme
    par['FuFi'] = dFuFiAll
    print('... finished processing external reference.')

    return par

def ProcessCompData(par):
    print('Start processing data ...')
    #> Loop Pandoras
    for sPanC in par['dPan']:
        #> Loop spectrometers
        for iSpec in par['dPan'][sPanC]:
            #> Loop dates
            for iDateC, sDateRefC in zip(par['iDate'], par['sRefDateTime']):
                #> Loop products
                lProds = [prod for prod,do in par['dProdAna'].items() if do]
                for sProd in lProds:
                    #> Loop processing f-codes for each Product:
                    for sFCode in par['dProdFCode'][sProd]:
                        #> Loop reference types:
                        for sRefType in par['dProdRefT'][sProd]:
                            for sCfSuffix in par['dCfSuffix'][sProd]:
                                ProcessingWrapper('Pandora', int(sPanC), iSpec, par['sLoc'], sFCode, par['sOFPth'], par['sCFPth'],
                                                  par['sBlickRootPth'], par['sL0Pth'], par['sL1Pth'], par['sL2Pth'],
                                                  par['sL2FitPth'], par['sPFPth'], par['iDate'][0],
                                                  par, par['sRefRtn'][0], iDateC, sDateRefC, par['dSCode']['s'],
                                                  sRefType, sCfSuffix)
    print('... finished processing data.')

# MAIN
########################################################################################################################
if __name__ == '__main__':

    par = {}
    par['delim'] = ' -> '
    par['epoch'] = datetime(2000,1,1)
    #> Read input file
    par['paramsFile'] = 'processCompDataInput.txt'
    par = LoadParams(par)

    #> Process external reference spectrum
    if par['doExtRef']:
        par = ProcessExternalReference(par)

    #> Process data
    if par['doProcData']:
        ProcessCompData(par)

    #> Convert data to CINDI-3 blind comparison format
    if par['doCINDI2BlindFmt']:
        ConvertToCINDI3BlindFmt(par)