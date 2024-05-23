# BlickP poor man's wrapper
# Alexander Cede, adapted by Martin Tiefengraber
# This "poor man's" wrapper does not do any error handling!

# IMPORT modules
import os
import h5py

# Define the desired path
path = os.path.normpath("C:/Blick/src")
os.chdir(path)
from glob import glob
from numpy import array, nan, where, argmax, ones, string_
from datetime import datetime, date, timedelta

from blick_psio import *
psio = blick_psio()
from blick_routinereader import *
rr = blick_routinereader()
import blick_processcore as prc
from blickp import create_processor_adv

def ProcessingWrapper(instrnam, instrnum, specnum, locname, fcode,
                    pthOF, pthCF, sBlickRootPth, pthL0, pthL1, pthL2, pthL2Fit, pthPF, currday,
                    par, sRefRtn, iDateC, sDateRefC, sSCode, sRefType, sCfSuffix):

    #Initialize empty dictionary
    dirs = {}

    # Example of setting values, assuming you have specific values to assign
    dirs['blick_dir'] = sBlickRootPth
    dirs['data_dir_L0'] = pthL0
    dirs['data_dir_L1'] = pthL1
    L2FitPth = os.path.join(pthL2Fit, sRefType)
    if not os.path.exists(L2FitPth):
        os.makedirs(L2FitPth)
    dirs['data_dir_L2Fit'] = L2FitPth
    dirs['data_dir_L2'] = pthL2
    dirs['data_dir_iof'] = sBlickRootPth + pthOF
    dirs['data_dir_icf'] = sBlickRootPth + pthCF

    from main_makeCompData import GetExternalReferenceFileName
    with h5py.File(pthPF, 'r+') as pssetup:
        idxFCode = where(pssetup['f_codes'][:,0]['f-code'] == fcode)[0][0]
        FCode = pssetup['f_codes'][idxFCode]
        if sRefType.startswith('Ref'):
            sRefNme = GetExternalReferenceFileName(par, instrnum, specnum, sRefRtn, FCode['Filter types'][0], iDateC, sDateRefC,
                                                   sSCode[0])
            sFCodeRefNme = string_('Ref_' + sRefNme)
        else:
            sFCodeRefNme = sRefType
        FCode['Reference'] = array([sFCodeRefNme], dtype='|S100')
        pssetup['f_codes'][idxFCode] = FCode

    # For the O3VIS product, the processor has to run for both calibration files (Finkenzeller, 2022; ThalmanAndVolkamer, 2013)
    cal_file_name = instrnam + str(instrnum) + 's' + str(specnum) + '_CF_' + sCfSuffix + '.txt'

    processor = create_processor_adv(dirs, instrnum, specnum, locname, None, fcode, None, cal_file=cal_file_name,
                                     proc_setup=pthPF, nuke_l2=True)

    dd = str(currday)
    dd = date(int(dd[:4]), int(dd[4:6]), int(dd[6:]))
    processor.process_day(dd)

    print("Poor man's wrapper finished")

    return curr_ref
