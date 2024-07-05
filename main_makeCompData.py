__author__ = 'Martin Tiefengraber'
import sys
import os
# Define the desired path
path = os.path.normpath("C:/Blick/src/")
os.chdir(path)
from os import path as ospath
from numpy import array, argmin, hstack, savetxt, where, ones, polyval, unique, asarray, in1d, nan, isnan, \
    nanmean, sum, nansum, isin, vstack
from datetime import datetime, timedelta
from matplotlib.dates import num2date, date2num
from params.InputParams import LoadParams
from GetDataOfRoutineBlick import GetDataOfRoutine
from ProcessingWrapper import ProcessingWrapperL1, ProcessingWrapper
from blick_countconverter import regrid_f
from glob import glob
import h5py
from pandas import Timedelta, Timestamp

from matplotlib.pyplot import *

from blick_io import blick_io
io = blick_io()

from LoadData import LoadData

from CINDI3SemiBlind import GetPths
from CINDI3SemiBlind import CINDI3SemiBlind
import CINDI3SemiBlind as CSB

def ConvertToCINDI3BlindFmt(par):

    print('Start converting to CINDI-3 blind format ...')
    dImp = {}
    dImp['loc'] = par['sLoc']
    dImp['sCode'] = par['dSCode']['s'][0]
    for sDate in par['iDate']:
        dImp['date'] = sDate
        # > Loop Pandoras
        for sPanC in par['dPan']:
            # > Loop spectrometers
            for iSpec in par['dPan'][sPanC]:
                dImp['panName'] = 'Pandora{}s{}'.format(sPanC, iSpec)
                dImp['institute'], dImp['instNumber'] = par['dPanIdInst'][dImp['panName'][7:]]
                # get L1 file psath
                _, _, sPthL1 = GetPths(dImp['date'], dImp['panName'], dImp['loc'], dImp['institute'],
                                                       dImp['instNumber'], '', '', '',
                                                       dImp['sCode'], '')
                # get version
                sPthL1 = glob(os.path.join(par['sL1Pth'], sPthL1))[-1]
                sVers = sPthL1.split(dImp['sCode'])[-1].split('p')[0][1:]
                # load l1 data
                l1 = LoadData(par['sL1Pth'], [dImp['date']], dImp['loc'], dImp['panName'].partition('Pandora')[-1],
                              'L1', sVers, [dImp['sCode']], CSB.colAssignL1, doXArray=True)[dImp['sCode']]
                #> Loop products
                lProds = [prod for prod,do in par['dProdAna'].items() if do]
                for sProd in lProds:
                    dImp['prod'] = sProd
                    dImp['fCode'] = par['dProdFCode'][sProd]
                    for sFCode in dImp['fCode']:
                        # > Loop reference types:
                        for sRefType in par['dProdRefT'][sProd]:
                            dImp['refType'] = sRefType
                            # load f-code details
                            dImp['fitPars'] = {}
                            with h5py.File(par['sPFPth'], 'r+') as pssetup:
                                idxFCode = where(pssetup['f_codes'][:, 0]['f-code'] == sFCode)[0][0]
                                FCode = pssetup['f_codes'][idxFCode]
                                for fitPar in ['npol', 'noffs', 'nwlc', 'nresc', 'WL-starts', 'WL-ends', 'Ring']:
                                    dImp['fitPars'][fitPar] = FCode[fitPar][0]
                                for fitPar in ['Fitted gases', 'Gas sources', 'Gas temps']:
                                    dImp['fitPars'][fitPar] = asarray(FCode[fitPar][0].split(','))
                            # overwrite O2O2 with O4, in case
                            dImp['fitPars']['Fitted gases'][dImp['fitPars']['Fitted gases'] == 'O2O2'] = array(['O4'],
                                                                                                               dtype='|S4')
                            dImp['prodVers'] = par['dProdVers'][sProd][0]
                            # select columns to be read from L2Fit file
                            colAssignUsed = {key:CSB.colAssignL2Fit[key] for key in par['dCompCols'][sProd]}
                            ## add temperature
                            for gas in dImp['fitPars']['Fitted gases']:
                                if gas == 'O2O2':
                                    gas = 'O4'
                                if gas == 'GLY':
                                    gas = 'CHOCHO'
                                gasT = '{}_T'.format(gas)
                                colAssignUsed[gasT] = CSB.colAssignL2Fit[gasT]
                            ## add auxilariy information
                            colAssignUsed['RTN'] = 'Two letter code of measurement routine'
                            colAssignUsed['RTNC'] = 'Routine count'
                            colAssignUsed['REPC'] = 'Repetition count'
                            colAssignUsed['ZPM'] = 'Zenith pointing mode'
                            colAssignUsed['APM'] = 'Azimuth pointing mode'
                            colAssignUsed['SZA'] = 'Solar zenith angle for center-time of measurement in degree'
                            colAssignUsed['SAA'] = 'Solar azimuth for center-time of measurement in degree'
                            colAssignUsed['PROCTYPE'] = 'Data processing type index'
                            ## TEMP
                            # colAssignUsed['lam1'] = 'Lower limit used for wavelength scaling'
                            # colAssignUsed['lam2'] = 'Upper limit used for wavelength scaling'
                            colAssignUsed['wvl0'] = 'Wavelength change polynomial coefficient, order 0'
                            colAssignUsed['wvl1'] = 'Wavelength change polynomial coefficient, order 1'
                            colAssignUsed['res0'] = 'Resolution change polynomial coefficient, order 0'
                            colAssignUsed['fitres'] = 'Fitting result index'
                            ##
                            dImp['procType'] = par['dProcType'][sRefType]
                            # get data version
                            _, sPthL2Fit, _ = GetPths(dImp['date'], dImp['panName'], dImp['loc'], dImp['institute'],
                                                      dImp['instNumber'], '', '', '',
                                                      '', sFCode)
                            sPthL2Fit = glob(os.path.join(par['sL2FitPth'], sRefType, sPthL2Fit))[-1]
                            sVers = sPthL2Fit.split(sFCode)[-1].split('p')[0][1:]
                            # load L2Fit data
                            l2fit = LoadData(os.path.join(par['sL2FitPth'], sRefType), [dImp['date']], dImp['loc'],
                                             dImp['panName'].partition('Pandora')[-1], 'L2Fit', sVers,
                                             [sFCode], colAssignUsed, doXArray=True)[sFCode]
                            if l2fit.time.size:
                                # reduce to allowed processing types
                                procTypes = [CSB.proctype2ref[procType] for procType in dImp['procType']]
                                l2fit = l2fit.isel(time=isin(l2fit.PROCTYPE, procTypes))
                                # remove values with fitres != 0
                                l2fit = l2fit.isel(time=isin(l2fit.fitres, [0]))
                                # get name of Cindi format file
                                sPthCindi, _, _ = GetPths(dImp['date'], dImp['panName'], dImp['loc'], dImp['institute'],
                                                          dImp['instNumber'], dImp['prod'], CSB.refTypeSyn[dImp['refType']], dImp['prodVers'], dImp['sCode'], sFCode)
                                # create CINDI-3 file
                                cindi = CINDI3SemiBlind(par['dCompCols'][sProd], par['fWvlINORM'], par['fWvlRef'], l1, l2fit, dImp,
                                                        sPthCindi, CSB.colAssignCindi3, par['dOvrwVZA'], par['dOvrwVAA'],
                                                        par['missValue'])
                                # get converted data columns and respective column description for the header
                                allCols, descrCols = cindi.buildupColumns(CSB.prodMainProd)
                                # get header retrieval settings
                                headRetrSet = cindi.headerRetrSet()
                                # get header for general information
                                headGeneral = cindi.headerGeneral(CSB.prodMainProd, CSB.refTypeSyn, par['comment'])
                                #build full header
                                allHeader = ''.join([line + '\n' for line in headGeneral + headRetrSet + descrCols])

                                # save file
                                savetxt(os.path.join(par['sCampPth'], sPthCindi), allCols, fmt='%.7e',
                                        comments='% ', delimiter=' ', header=allHeader[:-2])
                                print('   ... Pandora {}s{}, product {}, reftype {}, (fcode: {}), date {}'.format(
                                    int(sPanC), iSpec, sProd, sRefType, sFCode, sDate))

    print('... finished converting to CINDI-3 blind format.')

def GetExternalReferenceFileName(par, sPanC, iSpec, sRtn, sFuFi, iDate, sDateRef, iSCode):

    sRefNme = ospath.join(par['sRefPth'], 'RefPan{}s{}Rtn{}FuFi{}from{}pm{}min-for{}_s{}.txt'
                          .format(sPanC, iSpec, sRtn, sFuFi, sDateRef,
                                  par['iRefAvgInt'], iDate, iSCode))

    return sRefNme


def ProcessExternalReference(par):

    print('Start producing external reference ...')
    #> Loop Pandoras
    for sPanC in par['dPan']:
        #> Loop spectrometers
        for iSpec in par['dPan'][sPanC]:
            #> Loop dates
            for iDateC, sDateRefC in zip(par['iDate'], par['sRefDateTime']):
                iDateRefC = int(sDateRefC[:8])
                sCode = par['dSCode']['s'][0]
                print('   ... Pandora {}s{}, scode: {}, for date {}, using ref from date {}'.format(
                    int(sPanC), iSpec, sCode, iDateC, sDateRefC))
                # Process L1 data
                ProcessingWrapperL1('Pandora', int(sPanC), iSpec, par['sLoc'], sCode, par['sOFPth'], par['sCFPth'],
                                    par['sBlickRootPth'], par['sPFPth'], par['sL0Pth'], par['sL1Pth'], par['sL2Pth'],
                                    par['sL2FitPth'], iDateRefC, par['CfSuffixRef'])

                # load data
                # get L1 file psath
                _, _, sPthL1 = GetPths(iDateRefC, 'Pandora{}s{}'.format(sPanC, iSpec), par['sLoc'], '',
                                       '', '', '', '', sCode, '')
                # get version
                sPthL1 = glob(os.path.join(par['sL1Pth'], sPthL1))[-1]
                sVers = sPthL1.split(sCode)[-1].split('p')[0][1:]
                # load l1 data
                l1 = LoadData(par['sL1Pth'], [iDateRefC], par['sLoc'], '{}s{}'.format(int(sPanC), iSpec),
                              'L1', sVers, [sCode], CSB.colAssignL1, doXArray=True)[sCode]
                # restrict to time range
                start_time = Timestamp(sDateRefC) - Timedelta(minutes=float(par['iRefAvgInt']))
                end_time = Timestamp(sDateRefC) + Timedelta(minutes=float(par['iRefAvgInt']))
                l1 = l1.sel(time=slice(start_time, end_time))
                # limit to reference routine
                l1 = l1.isel(time=isin(l1.RTN, par['sRefRtn']))
                # limit to routines where wavelength change could be retrieved
                l1 = l1.isel(time=isin(l1.WVLfi, [0.]))
                # remove saturated data (good data SAT == 0)
                l1 = l1.isel(time=isin(l1.SAT, [0.]))
                # avoid measurements with only one cycle
                l1 = l1.isel(time=l1.NCYB > 1.)

                # crop wavlengths to about 535 nm.
                l1 = l1.sel(lam=slice(280, 535))

                # loop functional filter
                for pos, sFuFi in zip([1, 2], ['OPEN', 'U340']):
                    l1fufi = l1.isel(time=isin(l1.FW2ep, [pos]))
                    sRefNme = GetExternalReferenceFileName(par, sPanC, iSpec, par['sRefRtn'][0], sFuFi, iDateC, sDateRefC, sCode)

                    # measured wavelengths
                    dLam = asarray([polyval(vstack((l1fufi.WVL1.data, l1fufi.WVL0.data))[:, i], l1fufi.lam - 350.) for i in range(l1fufi.WVL0.shape[0])])
                    lamMeas = l1fufi.lam.data + dLam
                    lamMeasMean = lamMeas.mean(0)
                    # mean intensity
                    intMean = l1fufi.INORM.mean('time').data
                    # mean uncertainty
                    intUncMean = (((l1fufi.INORNu**2).sum('time') ** 0.5) / l1fufi.INORNu.shape[0]).data

                    # save reference
                    savetxt(sRefNme, vstack((lamMeasMean, intMean, intUncMean)).T)

    print('... finished producing external reference.')


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
                                print('   ... Pandora {}s{}, product {}, reftype {}, (fcode: {}), date {}'.format(
                                    int(sPanC), iSpec, sProd, sRefType, sFCode, iDateC))
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
        ProcessExternalReference(par)

    #> Process data
    if par['doProcData']:
        ProcessCompData(par)

    #> Convert data to CINDI-3 blind comparison format
    if par['doCINDI2BlindFmt']:
        ConvertToCINDI3BlindFmt(par)