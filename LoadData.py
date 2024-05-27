__author__ = 'Martin Tiefengraber'

from numpy import nan, loadtxt, argmin, datetime64, genfromtxt, unique, polynomial, argmax, asarray, ones, array, \
    where, isnan, arange, hstack, vstack, sum
from datetime import datetime, timedelta
from matplotlib.dates import datestr2num, num2date, date2num
from pandas import read_csv, to_datetime
from pytz import utc
from glob import glob
from os import path as ospath
from os.path import sep as ospathsep
from xarray import Dataset as xarrayDataset
from xarray import merge as xarraymerge
import codecs


def scale_x(x, xlim=[]):
    """
    FROM ALEXANDER CEDE, BLICK SOFTWARE SUITE
    xs=scale_x(x,xlim=[])

    <x> is a vector.
    This function scales <x> (output <xs>) using limits given in <xlim>.
    <xlim> can have 0, 1, or 2 elements.
        <xlim> has 0 elements: the scaled xs=x*scf
        <xlim> has 1 element: the scaled xs=scf*(x/xlim[0]-0.5)
        <xlim> has 2 elements: the scaled xs=scf*((x-xlim[0])/(xlim[1]-xlim[0])-0.5)
    <scf>=3.46 is the scale factor.
    """
    # scale factor
    scf = 3.46
    if xlim == []:
        xs = x * scf
    else:
        if len(xlim) == 1:
            xmin = 0
            xmax = xlim[0]
        else:
            xmin = xlim[0]
            xmax = xlim[1]
        xw = xmax - xmin
        xs = (array(x) - xmin) * scf / xw - scf * 0.5
    return xs

def LatestFile(lPth, bIsL2Tot=False):
    lVersC = []
    lVersP1 = []
    lVersP2 = []
    for sPthC in lPth:
        sPrt = sPthC.partition(ospathsep)[-1].split('_')[-1].partition('.txt')[0]
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

def GetDateArrayFromLimits(iStart, iEnd):
    dtStart = num2date(datestr2num(str(iStart))).date()
    dtEnd = num2date(datestr2num(str(iEnd))).date()
    iNumDays = (dtEnd - dtStart).days
    a1Date = ones(iNumDays + 1) * nan
    for dayI in range(iNumDays + 1):
        a1Date[dayI] = int((dtStart + timedelta(days=dayI)).isoformat().replace('-', ''))
    a1Date = a1Date.astype(int)

    return a1Date

def PickDataWithTimeRange(lVal, lValCmb, lRefDateTime, iRefOffs, sRefAvgInt, epoch=datetime(2000, 1, 1)):
    lValRef = []
    lValRefCmb = []
    lDtUsed = []
    for iDate in range(len(lVal)):
        lValRefDate = []
        lValRefCmbDate = []
        #> Generate date time vector
        ##> Datetime and remove time zone
        a1Dt = array([num2date((lVal[iDate][i][0][1][0] + date2num(epoch)*24.*3600.)/(24.*3600.)).replace(tzinfo=None)
                      for i in range(len(lVal[iDate]))])
        #> Get datetime range
        ##> Reference date
        a1DtRef = array([datetime.strptime(lRefDateTime[i], '%Y%m%dT%H%M%SZ')
                         for i in range(len(lRefDateTime))])
        if a1DtRef.shape[0] == 1:  # use the same reference for all days
            iDate = 0
        #>  Get averaging interval
        ##> Check maximum allowed distance to reference time
        a1DtRefBest = a1Dt[argmin(abs(a1Dt - a1DtRef[iDate]))]
        assert (a1DtRefBest >= a1DtRef[iDate] - timedelta(minutes=float(iRefOffs))) & \
               (a1DtRefBest <= a1DtRef[iDate] + timedelta(minutes=float(iRefOffs))), \
            'Wanted reference time (from {}) not within maximum boundaries (p/m {} min.)!'. \
            format(lRefDateTime[iDate], float(iRefOffs))
        ##> Get averagint interval
        bId = (a1Dt >= a1DtRefBest - timedelta(minutes=float(sRefAvgInt))) & \
              (a1Dt <= a1DtRefBest + timedelta(minutes=float(sRefAvgInt)))
        ##> Check if averaging interval contains data
        assert any(bId), \
            'Wanted closest reference time (from {}, p/m {} min.) not included in the data!'. \
            format(a1DtRefBest, float(sRefAvgInt))
        #> Extract data for reference and average
        for iRtn in where(bId)[0]:
            lValRefDate.append(lVal[iDate][iRtn])
            lValRefCmbDate.append(lValCmb[iDate][iRtn])
        lValRef.append(lValRefDate)
        lValRefCmb.append(lValRefCmbDate)
        lDtUsed.append(a1Dt[bId])

    return lValRef, lValRefCmb, lDtUsed

def GetColumnNumber(sPth, sColumnDescr, sHeadSep='---', iHeadSepCnt=2):
    with codecs.open(sPth, 'r', encoding='ISO-8859-1') as fle:
        iHeadSep = 0
        iColumnNumber = 0
        for sLine in fle:
            if sHeadSep in sLine:
                iHeadSep += 1
                if iHeadSep == iHeadSepCnt:
                    # print 'Data line "{}" not found!'.format()
                    iColumnNumber = nan
                    break
            if sColumnDescr.lower() in sLine.lower():
                sColumnNumber = sLine.partition(':')[0].split()[1]
                if '-' in sColumnNumber:
                    iColumnNumber = [int(sColumnNumber.partition('-')[0]) - 1, int(sColumnNumber.partition('-')[2]) - 1]
                else:
                    iColumnNumber = int(sColumnNumber) - 1
                break
            if iHeadSep:
                iColumnNumber += 1

    return iColumnNumber

def GetNumberHeaderLines(sPth, sHeadSep='---', iHeadSepCnt=2):
    iLine = None
    with codecs.open(sPth, 'r', encoding='ISO-8859-1') as fle:
        iHeadSep = 0
        for iLine, sLine in enumerate(fle):
            if sHeadSep in sLine:
                iHeadSep += 1
                if iHeadSep == iHeadSepCnt:
                    break

    return iLine + 1

def LoadData(sPth, lDates, sLoc, sPan, sLev, sVers, lProcCode, dColAsign, doXArray=False):

    """
    ugly, but works
    """

    #> Append scale factor column for unscaling if L1 spectra are requested
    if sLev == 'L1' and (any([sValue.startswith('Level 1 data') for sValue in dColAsign.values()]) or
                         any([sValue.startswith('L1 data') for sValue in dColAsign.values()])):
        dColAsign['scl'] = 'Scale factor for data'
    if sLev == 'L0' and any([sValue.startswith('Mean over all cycles') for sValue in dColAsign.values()]):
        dColAsign['scl'] = 'Scale factor for data'
    #> Append start and end wavelength for polynomial evaluation for L2Fit files, if requested
    splitCoef = 'polynomial coefficient'
    lPolCoefs = [c for c in dColAsign.values() if splitCoef in c]  # all columns with polynomial coefficients
    sSgnLam1, sSgnLam2 = None, None
    if (sLev == 'L2Fit') and (len(lPolCoefs) > 0):
        sSgnLam1 = 'Lower limit used for wavelength scaling [nm]'
        sSgnLam2 = 'Upper limit used for wavelength scaling [nm]'
        dColAsign['lam1'] = sSgnLam1
        dColAsign['lam2'] = sSgnLam2
        lPolCoefs = sorted(lPolCoefs)
    lCols = sorted(dColAsign.keys())
    #> Get full date vector
    if sLev in ['L2']:
        a1Dates = array([nan])
    else:
        if len(lDates) > 1:
            a1Dates = GetDateArrayFromLimits(*lDates)
        else:
            a1Dates = asarray(lDates)
    #> Initialize output variable
    dxData = None
    if doXArray:
        dxData = {}
    dRtn = {}
    dData = {}
    dMeta = {}
    dGetCols = {}
    #> Loop data
    lGetColsIdx = {}
    for sCode in lProcCode:
        for iDateI, iDateC in enumerate(a1Dates):
            dGetCols[sCode] = []
            #> Read data
            if sLev in ['L2']:
                sPthDat = 'Pandora{}_{}_{}_?{}p?-?.txt'.format(sPan, sLoc, sLev, sCode)
            elif sLev in ['L0']:
                sPthDat = 'Pandora{}_{}_{}_{}.txt'.format(sPan, sLoc, iDateC, sLev)
            else:
                if sVers == 'LATEST':
                    sPthDat = 'Pandora{}_{}_{}_{}_?{}c*p?-?.txt'.format(sPan, sLoc, iDateC, sLev, sCode)
                else:
                    sPthDat = 'Pandora{}_{}_{}_{}_?{}c{}p?-?.txt'.format(sPan, sLoc, iDateC, sLev, sCode, sVers)
            lPthDat = glob(ospath.join(sPth, sPthDat))
            if len(lPthDat) > 0:
                if sLev in ['L0']:
                    sPthDat = lPthDat
                else:
                    sPthDat = LatestFile(lPthDat, bIsL2Tot=sLev in ['L2'])
                for sPthGlob in sPthDat:
                    # > Get columns in data file to be extracted
                    # lGetCols = []
                    for sColC in lCols:
                        iCol = GetColumnNumber(sPthGlob, dColAsign[sColC])
                        if array(isnan(iCol)).any():
                            print('Could not find "{}" in \n{}'.format(dColAsign[sColC], sPthGlob))
                        else:
                            dGetCols[sCode].append(iCol)
                    #> Extract all columns (from column limits)
                    lGetColsFull = []
                    lGetColsIdx[sCode] = []
                    for colC in dGetCols[sCode]:
                        len1 = len(lGetColsFull)
                        if type(colC) == list:
                            newCols = list(range(colC[0], colC[1] + 1))
                            lGetColsFull = lGetColsFull + newCols
                            lGetColsIdx[sCode].append(list(range(len1, len1+len(newCols))))
                        else:
                            lGetColsFull.append(colC)
                            lGetColsIdx[sCode].append([len1])
                    #> Load data
                    # TODO make converters more global!
                    iHeadLines = GetNumberHeaderLines(sPthGlob)
                    iMetaLines = GetNumberHeaderLines(sPthGlob, iHeadSepCnt=1)
                    a1Rtn = None
                    #> Get meta data
                    with codecs.open(sPthGlob, encoding='ISO-8859-1') as fle:
                        lLines = fle.readlines()
                        for iLine in arange(iMetaLines - 1):
                            sLine = lLines[iLine]
                            sKey, _, sValue = sLine.partition(': ')
                            dMeta[sKey] = sValue[:-1]
                    if sLev in ['L2']:
                        a2Dat = loadtxt(sPthGlob,
                                        converters={0: datestr2num},
                                        skiprows=iHeadLines, usecols=lGetColsFull, encoding='latin_1')
                    elif sLev in ['L0']:
                        # Get comment rows to be skipped
                        toskip = []
                        f = open(sPthGlob, encoding='ISO-8859-1')
                        for i, l in enumerate(f.readlines()):
                            if i > (iHeadLines - 1):
                                if l.split(' ')[4] == '#':
                                    toskip.append(i)
                        toskip = list(range(iHeadLines)) + toskip
                        a2Dat = read_csv(sPthGlob, skiprows=toskip, usecols=lGetColsFull, sep=' ',
                                         converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num})
                        a1Rtn = genfromtxt(sPthGlob, dtype='S2', skip_header=iHeadLines, usecols=[0], encoding='ISO-8859-1')
                    else:
                        a2Dat = loadtxt(sPthGlob,
                                        converters={0: lambda s: array([ord(c)**2 for c in s]).sum(), 1: datestr2num},
                                        skiprows=iHeadLines, usecols=lGetColsFull, encoding='latin_1')
                        a1Rtn = genfromtxt(sPthGlob, dtype='S2', skip_header=iHeadLines, usecols=[0], encoding='ISO-8859-1')
                    if iDateI == 0:
                        dRtn[sCode] = array([])
                        dData[sCode] = array([]).reshape((0, len(lGetColsFull)))
                    if len(a2Dat) > 0:
                        dRtn[sCode] = hstack((dRtn[sCode], a1Rtn))
                        dData[sCode] = vstack((dData[sCode], a2Dat))
            else:
                print('{} not available!'.format(ospath.join(sPth, sPthDat)))
    #> Stop of no files have been found
    assert dData, 'No files found!'

    sDtSignal = None
    sLamSignal = None
    if sLev == 'L1':
        sDtSignal = 'UT date and time for beginning of measurement'
        sLamSignal = 'Nominal wavelengths [nm]'
    elif sLev == 'L0':
        sDtSignal = 'UT date and time for beginning of measurement'
        sLamSignal = 'Pixels'
    elif sLev in ['L2Fit']:
        sDtSignal = 'UT date and time for center-time of measurement'
        sLamSignal = 'Nominal wavelengths inside fitting window [nm]'
    elif sLev in ['L2']:
        sDtSignal = 'UT date and time for measurement center'
        sLamSignal = 'Nominal wavelengths inside fitting window [nm]'
    #> Get output variable column assigment
    dCols = {}
    for iColI, sColC in enumerate(lCols):
        dCols[sColC] = iColI
    if doXArray:
        #> Convert output data to dict
        for sCode in dData.keys():
            dCoords = {}
            for sDimsI, sDimsC in enumerate(lCols):
                idx = lGetColsIdx[sCode][sDimsI]
                #> Convert date time
                if dColAsign[sDimsC].startswith(sDtSignal):
                # if sDtSignal in dColAsign[sDimsC]:
                    a1Dt = num2date(dData[sCode][:, idx].squeeze(), tz=utc)
                    dCoords[sDimsC] = array([datetime64(dt) for dt in a1Dt])
                #> Use orignal routine names (strings instead of floats)
                elif dColAsign[sDimsC].startswith('Two letter code of measurement routine'):
                    dCoords[sDimsC] = dRtn[sCode]
                #> Unscale spectra
                elif dColAsign[sDimsC].startswith('Level 1 data for each pixel') or \
                     dColAsign[sDimsC].startswith('L1 data for each pixel') or \
                     dColAsign[sDimsC].startswith('Uncertainty of level 1 data for each pixel') or \
                     dColAsign[sDimsC].startswith('Mean over all cycles of raw counts for each pixel') or \
                     dColAsign[sDimsC].startswith('Uncertainty of raw counts for each pixel') or \
                     dColAsign[sDimsC].startswith('Instrumental uncertainty of level 1 data for each pixel') or \
                     dColAsign[sDimsC].startswith('Independent instrumental uncertainty of L1 data'):
                    a2Spec = dData[sCode][:, idx]
                    a1Scl = dData[sCode][:, lGetColsIdx[sCode][lCols.index('scl')]].squeeze()
                    dCoords[sDimsC] =  (a2Spec.T / a1Scl).T
                #> Unscale residuals
                elif 'slant column residuals' in dColAsign[sDimsC]:
                    a2Res = dData[sCode][:, idx]
                    dCoords[sDimsC] =  a2Res / 1e5
                else:
                    dCoords[sDimsC] = dData[sCode][:, idx]
                dCoords[sDimsC] = dCoords[sDimsC].squeeze()
            dData[sCode] = dCoords
        a1Lam = None
        idxDt = array([val.startswith(sDtSignal) for val in dColAsign.values()])
        if any(idxDt):
            sDtPar = array(list(dColAsign.keys()))[idxDt][0]
            dData[sCode][sDtPar] = to_datetime(dData[sCode][sDtPar]).astype('datetime64[ns]')
            for sCode in dData.keys():
                xDatLam = None
                if (max([dData[sCode][s].ndim for s in dData[sCode]]) > 1) or (len(lPolCoefs) > 0):
                    if sLev in ['L0']:
                        a1Lam = unique(array([c[1]-c[0]+1 for c in dGetCols[code] if type(c) == list if len(c) > 1 for code in dGetCols]))  # extract pixels
                    else:
                        a1Lam = array(list(map(float, dMeta[sLamSignal].split())))
                    xDatLam = xarrayDataset(coords={'time': dData[sCode][sDtPar], 'lam': a1Lam}, attrs=dMeta)
                xDat = xarrayDataset(coords={'time': dData[sCode][sDtPar]}, attrs=dMeta)
                for sPar in dData[sCode].keys():
                    if dData[sCode][sPar].ndim == 1:
                        xDat[sPar] = (['time'], dData[sCode][sPar])
                    else:
                        xDatLam[sPar] = (['time', 'lam'], dData[sCode][sPar])
                # Add evaluated polynomials if requested
                if len(lPolCoefs) > 0:
                    dColAsignInv = {v: k for k, v in dColAsign.items()}
                    # Get scaled wavelength
                    lam1 = float(xDat[dColAsignInv[sSgnLam1]].mean().data)
                    lam2 = float(xDat[dColAsignInv[sSgnLam2]].mean().data)
                    a1LamScl = scale_x(a1Lam, [lam1, lam2])
                    # Get parameters to be evaluated
                    lPolParsUnq = unique(array([c.split(splitCoef)[0] + splitCoef for c in lPolCoefs]))
                    lPolPars = array([c.split(splitCoef)[0] + splitCoef for c in lPolCoefs])
                    lPolParsAll = []
                    for polPar in lPolParsUnq:
                        lPolParsAll.append(array(lPolCoefs)[lPolPars == polPar])
                    # Evaluate polyomials
                    for polPar in lPolParsAll:
                        a1Dat = xDat[dColAsignInv[polPar[0]]].data[None, :]
                        if len(polPar) > 1:
                            for pol in polPar[1:]:
                                a1Dat = vstack((a1Dat, xDat[dColAsignInv[pol]]))
                        anDat = polynomial.polynomial.polyval(a1LamScl, a1Dat)
                        xDatLam[dColAsignInv[polPar[0]][:-1]] = (['time', 'lam'], anDat)
                if xDatLam is not None:
                    dxData[sCode] = xarraymerge([xDat, xDatLam])
                else:
                    dxData[sCode] = xDat
        else:
            dxData = dData

    if doXArray:
        return dxData
    else:
        return dData, dCols