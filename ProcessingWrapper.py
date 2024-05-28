import os
import h5py
import gc
import shutil
from copy import copy

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


def close_h5_files(file_path):
    # Iterate over all objects tracked by the garbage collector
    for obj in gc.get_objects():
        try:
            # Check if the object is an h5py File instance
            if isinstance(obj, h5py.File):
                # In Python 2.7, the filename can be accessed via obj.file.filename
                if obj.file.filename == file_path:
                    obj.close()
        except AttributeError:
            # Handle case where obj.file or obj.file.filename does not exist
            pass
        except ValueError:
            # Handle case where file ID is not valid
            pass


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
    # copy orginal processing setups file to temporary setups file
    pthPForg = copy(pthPF)
    pthPFtmp = ''.join(pthPForg.partition('ProcessingSetups')[:2]) + '_temp.h5'
    shutil.copyfile(pthPForg, pthPFtmp)
    # overwrite h5 file
    with h5py.File(pthPFtmp, 'r+') as pssetup:
        idxFCode = where(pssetup['f_codes'][:,0]['f-code'] == fcode)[0][0]
        FCode = pssetup['f_codes'][idxFCode]
        if sRefType.startswith('Ref'):
            sRefNme = GetExternalReferenceFileName(par, instrnum, specnum, sRefRtn, FCode['Filter types'][0], iDateC, sDateRefC,
                                                   sSCode[0])
            sFCodeRefNme = string_('Ref_' + sRefNme)
        else:
            sFCodeRefNme = sRefType
        FCode['Reference'] = array([sFCodeRefNme], dtype='|S300')
        pssetup['f_codes'][idxFCode] = FCode

    # For the O3VIS product, the processor has to run for both calibration files (Finkenzeller, 2022; ThalmanAndVolkamer, 2013)
    cal_file_name = instrnam + str(instrnum) + 's' + str(specnum) + '_CF_' + sCfSuffix + '.txt'

    processor = create_processor_adv(dirs, instrnum, specnum, locname, None, fcode, None, cal_file=cal_file_name,
                                     proc_setup=pthPFtmp, nuke_l2=True)

    dd = str(currday)
    dd = date(int(dd[:4]), int(dd[4:6]), int(dd[6:]))
    processor.process_day(dd)

    # close processing setups file
    close_h5_files(pthPFtmp)