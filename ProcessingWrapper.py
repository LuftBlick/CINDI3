# BlickP poor man's wrapper
# Alexander Cede, adapted by Martin Tiefengraber
# This "poor man's" wrapper does not do any error handling!

# IMPORT modules
import os
from blickp_processingsetups import FittingSetup
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
                    par, sRefRtn, iDateC, sDateRefC, sSCode, sCfSuffix):

    #Initialize empty dictionary
    dirs = {}

    # Example of setting values, assuming you have specific values to assign
    dirs['blick_dir'] = sBlickRootPth
    dirs['data_dir_L0'] = pthL0
    dirs['data_dir_L1'] = pthL1
    dirs['data_dir_L2Fit'] = pthL2Fit
    dirs['data_dir_L2'] = pthL2
    dirs['data_dir_iof'] = sBlickRootPth + pthOF
    dirs['data_dir_icf'] = sBlickRootPth + pthCF

    # For the O3VIS product, the processor has to run for both calibration files (Finkenzeller, 2022; ThalmanAndVolkamer, 2013)
    cal_file_name = instrnam + str(instrnum) + 's' + str(specnum) + '_CF_' + sCfSuffix + '.txt'
    processor = create_processor_adv(dirs, instrnum, specnum, locname, None, fcode, None, cal_file=cal_file_name,
                                     proc_setup=pthPF, nuke_l2=True)

    reference = processor.config.processing_setups.fcodes[fcode].reference
    # Check if reference starts with REF_ and overwrite it with reference data
    if reference.startswith('REF_'):
        # Overwrite product reference type in f-code
        copy = processor.config.processing_setups.fcodes[fcode]

        from main_makeCompData import GetExternalReferenceFileName
        sRefNme = GetExternalReferenceFileName(par, instrnum, specnum, sRefRtn, copy.filter_types, iDateC, sDateRefC, sSCode[0])

        processor.config.processing_setups.fcodes[fcode] = FittingSetup(
            code=copy.code,
            name=copy.name,
            processing_types=copy.processing_types,
            filter_types=copy.filter_types,
            reference=string_('REF_' + sRefNme), # Overwrite product reference type in f-code
            wlstart=copy.wlstart,
            wlend=copy.wlend,
            npoll=copy.npoll,
            noffs=copy.noffs,
            nwlc=copy.nwlc,
            nresc=copy.nresc,
            fitted_gases=copy.fitted_gases,
            gas_sources=copy.gas_sources,
            gas_od_methods=copy.gas_od_methods,
            gas_temperatures=copy.gas_temperatures,
            ring=copy.ring,
            molecular_scattering=copy.molecular_scattering,
            linear_fitting=copy.linear_fitting,
            uncertainty=copy.uncertainty,
            scodes=copy.scodes,
            diffuse_correction=copy.diffuse_correction,
            time_interploation=copy.time_interploation,
            pixels_to_use=copy.pixels_to_use,
            qf_code=copy.qf_code,
            date=copy.date,
            author=copy.author)

    else:
        print('Reference not REF_')

    dd = str(currday)
    dd = date(int(dd[:4]), int(dd[4:6]), int(dd[6:]))
    processor.process_day(dd)

    curr_ref = str(reference)
    print("Poor man's wrapper finished")

    return curr_ref