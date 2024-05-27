__author__ = 'Martin_Tiefengraber'

from numpy import array
from GetLineEntry import GetLineEntry


def LoadParams(par):

    #> ['<signal name in parameter file>', '<variable name>', '<parameter data type>']
    params = array([
        #> processing parameters
        ['Pandora unit(s) with spectrometers', 'dPan', 'dictint'],
        ['Locations', 'sLoc', 'string'],
        ['Date(s) to process [UTC]', 'iDate', 'int'],
        ['Spectral fitting reference datetime [UTC]', 'sRefDateTime', 'stringarray'],
        ['Average interval around reference datetime [min]', 'iRefAvgInt', 'int'],
        ['Maximum allowed time delay of reference time to nearest measured spectrum [min]', 'iRefOffs', 'int'],
        ['Intensity variation filter in percent', 'varFilt', 'float'],
        ['Spectral fitting reference from routine', 'sRefRtn', 'stringarray'],
        ['Spectral fitting reference from routine number of measurements for routine', 'iRtnCnt', 'int'],
        ['L1 processing s-code', 'dSCode', 'dictstring'],
        ['Calibration file suffix for processing reference', 'CfSuffixRef', 'string'],
        ['Intensity wavelength average pm [nm]', 'fWvlIntAvg', 'float'],
        # > paths
        ['Blick root directory', 'sBlickRootPth', 'string'],
        ['Operation file directory', 'sOFPth', 'string'],
        ['Calibration file directory', 'sCFPth', 'string'],
        ['Pslib directory', 'sPslibPth', 'string'],
        ['Reference data directory', 'sRefPth', 'string'],
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
        ['Product version number', 'dProdVers', 'dictint'],
        ['Calibration file suffixes', 'dCfSuffix', 'dictstring'],
        ['L2Fit file number of header lines', 'dProdNHead', 'dictint'],
        ['L1 file number of header lines', 'dL1NHead', 'dictint'],
        ['L1 file start of spectra column', 'dL1CcCol', 'dictint'],
        ['Institution and instrument number assignment', 'dPanIdInst', 'dictstring'],
        ['Analyze product', 'dProdAna', 'dictbool'],
        ['Columns in comparison data', 'dCompCols', 'dictstring'],
        ['Wavelengths for INORM', 'fWvlINORM', 'float'],
    ])
    #> read processing parameters
    for parC in params:
        par[parC[1]] = GetLineEntry(par['paramsFile'], parC[0], par['delim'], parC[2])

    #> reshape parameters
    for parC in ['iDate', 'sRefDateTime', 'sRefRtn', 'iRtnCnt']:
        if par[parC].shape == ():
            par[parC] = par[parC].reshape(1, )

    return par