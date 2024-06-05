__author__ = 'Martin_Tiefengraber'

from numpy import array
from params.GetLineEntry import GetLineEntry


def LoadParams(par):

    #> ['<signal name in parameter file>', '<variable name>', '<parameter data type>']
    params = array([
        #> processing parameters
        ['Pandora unit(s) with spectrometers', 'dPan', 'dictint'],
        ['Locations', 'sLoc', 'string'],
        ['Date(s) to process [UTC]', 'iDate', 'int'],
        ['Spectral fitting reference datetime [UTC]', 'sRefDateTime', 'stringarray'],
        ['Average interval around reference datetime [min]', 'iRefAvgInt', 'int'],
        ['Spectral fitting reference from routine', 'sRefRtn', 'stringarray'],
        ['L1 processing s-code', 'dSCode', 'dictstring'],
        ['Calibration file suffix for processing reference', 'CfSuffixRef', 'string'],
        ['Intensity wavelength average pm [nm]', 'fWvlIntAvg', 'float'],
        # > paths
        ['Blick root directory', 'sBlickRootPth', 'string'],
        ['Operation file directory', 'sOFPth', 'string'],
        ['Calibration file directory', 'sCFPth', 'string'],
        ['Pslib directory', 'sPslibPth', 'string'],
        ['Reference data directory', 'sRefPth', 'string'],
        ['Figure directory', 'sPlotPth', 'string'],
        ['L0 data directory', 'sL0Pth', 'string'],
        ['L1 data directory', 'sL1Pth', 'string'],
        ['L2 data directory', 'sL2Pth', 'string'],
        ['L2Fit data directory', 'sL2FitPth', 'string'],
        ['Allowed processing types per reference', 'dProcType', 'dictstring'],
        ['Processing setup file directory', 'sPFPth', 'string'],
        ['Campaign data directory', 'sCampPth', 'string'],
        # > processing steps
        ['Calculate external reference', 'doExtRef', 'bool'],
        ['Process data', 'doProcData', 'bool'],
        ['Convert to CINDI3 blind comparison format', 'doCINDI2BlindFmt', 'bool'],
        ['Overwrite nominal viewing zenith angles', 'doOvrwVZA', 'bool'],
        ['Viewing zenith angles to overwrite', 'dOvrwVZA', 'dictfloat'],
        ['Viewing azimuth angles to overwrite', 'dOvrwVAA', 'dictfloat'],
        ['Product aliases for f codes', 'dProdFCode', 'dictstring'],
        ['Product reference type', 'dProdRefT', 'dictstring'],
        ['Product version number', 'dProdVers', 'dictfloat'],
        ['Calibration file suffixes', 'dCfSuffix', 'dictstring'],
        ['Institution and instrument number assignment', 'dPanIdInst', 'dictstring'],
        ['Analyze product', 'dProdAna', 'dictbool'],
        ['Columns in comparison data', 'dCompCols', 'dictstring'],
        ['Wavelengths for INORM', 'fWvlINORM', 'float'],
        ['Reference wavelength for horizon scans', 'fWvlRef', 'float'],
        ['Missing value', 'missValue', 'float'],
        ['Comment to be added in CINDI3 format', 'comment', 'string'],
    ])
    #> read processing parameters
    for parC in params:
        par[parC[1]] = GetLineEntry(par['paramsFile'], parC[0], par['delim'], parC[2])

    #> reshape parameters
    for parC in ['iDate', 'sRefDateTime', 'sRefRtn']:
        if par[parC].shape == ():
            par[parC] = par[parC].reshape(1, )

    return par
