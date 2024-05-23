__author__ = 'Martin Tiefengraber'

import sys
import os
# Define the desired path
path = os.path.normpath("C:/Blick/src")
os.chdir(path)
sys.path.insert(0, path)
from numpy import dstack
from blick_psio import *
psio = blick_psio()

#trace gas table and processing setups
psfname = 'C:/Blick/lib/pslib/Blick_ProcessingSetups.h5'
res,gasinfo,sinfo,qsinfo,finfo,qfinfo,rinfo,qrinfo=psio.get_procsetups(psfname)
assert res == 'OK', res

#general calibration parameters
psp.LoadCalibrationParameters(gasinfo)

class GetDataOfRoutine:

    def __init__(self, iPan, iSpec, sLoc, sPthRoot, sPthData, sPthOF, sPthCF, sCfSuffix, lBlickCodes):
        self.iPan = iPan
        self.iSpec = iSpec
        self.sLoc = sLoc
        self.sPthRoot = sPthRoot
        self.sPthData = sPthData
        self.sPthOF = path.join(self.sPthRoot, sPthOF)
        self.sPthCF = path.join(self.sPthRoot, sPthCF)
        self.sCfSuffix = sCfSuffix
        self.lBlickCodes = lBlickCodes  # elem 0 = s number, elem 1 = f number, elem 2 = r number
        #self.lBlickQCodes = lBlickQCodes  # elem 0 = qs number, elem 1 = qf number, elem 2 = qr number


    def GetCalData(self, dtDate):

        instrnam = 'Pandora'
        # read L0 data file header and reset read-info
        sDate = datetime.strftime(dtDate, '%Y%m%d')
        l0fname = instrnam + str(self.iPan) + 's' + str(self.iSpec) + '_' + self.sLoc + '_' + sDate + '_L0.txt'
        res, fileind, L0File = io.open_file(path.join(self.sPthData, l0fname), "r")
        res, l0medata, col0ind, a0line = io.read_dataheader(L0File, gp.lev0col)
        # L0File.close()

        # get corresponding IOF file
        indof = xfus.get_index("Instrument operation file used", l0medata, 0)
        iofname = path.join(self.sPthOF, l0medata[indof][1])
        assert path.isfile(iofname), 'Operationfile "{}" not available!'.format(iofname)

        # read IOF (only read if it has changed relative to the previous day!)
        res, spec_pars, hst_pars, tc_pars, cam_pars, pssys_pars, sbhs_pars, meta_data, checksum = io.get_instrparams(iofname)

        # load ICF
        icfname = os.path.join(self.sPthCF + '{}{}s{}_CF_{}.txt'.format(instrnam, self.iPan, self.iSpec, self.sCfSuffix))
        assert path.isfile(icfname), 'Calibrationfile "{}" not available!'.format(icfname)
        print('   ... CF: {}'.format(icfname))

        # read ICF (only read if it has changed relative to the previous day!)
        ##> change working directory to src directory
        # general calibration parameters
        icfname = path.join(self.sPthCF, icfname)
        res, cal_data = psio.get_calparams(icfname)

        base_name, _ = os.path.splitext(icfname)  # Removes the extension
        ddfname = "{}.ddf".format(path.join(self.sPthCF, base_name))
        if not os.path.exists(ddfname):
            # build last element of cal_data
            res, lCalData = psio.get_lastcalparam(cal_data, spec_pars, hst_pars, tc_pars, cam_pars, pssys_pars, sbhs_pars, gasinfo, changeels=True, mkslm=True)
        else:
            lCalData = io.xload(ddfname)
        # get filterwheel map without counting DIFF
        ispec = self.iSpec - 1
        fwmap = rr.get_fwmap(hst_pars, False)[ispec]

        npix = spec_pars[self.iSpec-1][5][0]
        pix = arange(npix) + 1
        pixs = xfus.scale_x(pix, [npix])   # scaled pixels
        a1Wvl = polyval(lCalData[2][1], pixs)
        if len(spec_pars[self.iSpec-1][5]) == 1:
            nblindpix = 0
        else:
            nblindpix = spec_pars[self.iSpec - 1][5][1]
        a1Wvl = a1Wvl[0 : spec_pars[self.iSpec-1][5][0]-nblindpix]

        return lCalData, a1Wvl, fwmap, [spec_pars, hst_pars, tc_pars, cam_pars, pssys_pars, sbhs_pars]


    def GetDataOfRoutine(self, sRtn, dtDates, iLev):

        #> Load data column description
        lLevCols = None
        sLev = None
        if iLev == "0":
            lLevCols = gp.lev0col
            sLev = 'L0'
        elif iLev == "1":
            lLevCols = psp.lev1col
            sLev = 'L1'
        elif iLev == "2Fit":
            lLevCols = psp.lev3col
            sLev = 'L2Fit'
        sysn = 'Pandora'

        #> Load data of routine
        ##> Initialize output vector
        lDataCmb = []
        lData = []
        lCom = []
        for dtDateC in dtDates:
            sDate = datetime.strftime(dtDateC, '%Y%m%d')
            sPthFle = path.join(self.sPthData, '{}{}s{}_{}_{}_{}.txt'.format(sysn, self.iPan, self.iSpec, self.sLoc, sDate, sLev))

            res, lLMeData, colind, aline = io.read_dataheader(sPthFle, lLevCols)
            assert res == 'OK', res
            res, lA, readinfo = io.read_data(sPthFle, colind, lLevCols, aline, [sRtn])
            if res != 'OK':
                print res

            lDataCmbRtns = []
            lDataRtns = []
            lComRtns = []
            for lAC in lA:
                lDataRtns.append(lAC)
                res, aa, comment = io.combine_data(lAC, lLevCols)
                lDataCmbRtns.append(aa)
                lComRtns.append(comment)
            lDataCmb.append(lDataCmbRtns)
            lData.append(lDataRtns)
            lCom.append(lComRtns)

        return lDataCmb, lData, lCom


    def DoDataCorrectionOneDate(self, iNumPix, iRoutCnt, lDataCmb, lData, dtDates):

        # > Initialize Pandora functions
        import blick_countconverter as coco

        #> s/qs-number index
        #qsind = qsinfo['qs-code'] == self.lBlickQCodes[0]
        sind = sinfo['s-code'] == self.lBlickCodes[0]

        lWvl = []
        lCc = []
        lECc = []
        lCcInfo = []
        for iDayI in xrange(len(dtDates)):
            ##> Load calibration data
            lCalData, a1Wvl, fwmap, [spec_pars, hst_pars, tc_pars, cam_pars, pssys_pars, sbhs_pars] = self.GetCalData(dtDates[iDayI])
            iNumPixReg = lCalData[-1][3][0].sum()
            if iNumPixReg != iNumPix:
                iNumPix = iNumPixReg
            a3Cc = ones((iRoutCnt, iNumPix, 0))
            a3ECc = ones((iRoutCnt, iNumPix, 0))
            lCcInfoDay = []
            #> Loop individual measurement within date
            l1colfilt, sinfoc, snumnotes = psio.get_lev1col(lCalData, spec_pars, hst_pars, tc_pars, cam_pars, pssys_pars, sbhs_pars, sinfo[sind])
            for iRtnI in xrange(len(lData[iDayI])):
                #> Apply data correction steps
                checkwlc = True
                if lData[iDayI][iRtnI][0][2][11] == 8:  # do not check wavelength if a lamp signal was measured
                    checkwlc = False
                ##> correct counts
                # a, cal_data, gasinfo, sinfo, qsinfo, checkwlc=True, islamp=False, flux=None
                a2Cc, a2ECc, ecc_instr, ccinfo = coco.convert_rc(lDataCmb[iDayI][iRtnI], lCalData, gasinfo, sinfo[sind],
                                                                 checkwlc)
                ##> Stack data
                a3Cc = dstack((a3Cc, a2Cc))
                a3ECc = dstack((a3ECc, a2ECc))
                lCcInfoDay.append(ccinfo)
            lWvl.append(a1Wvl)
            lCc.append(a3Cc)
            lECc.append(a3ECc)
            lCcInfo.append(lCcInfoDay)

        return lWvl, lCc, lECc, lCcInfo


if __name__ == '__main__':

    #> For testing
    dtDates = array([datetime(2016, 7, 12)])
    GD = GetDataOfRoutine(110, 1, "LabIBK", "C:/Blick/", "data/L0/", "data/operationfiles/", "data/calibrationfiles/",
                          [4, 0, 0])
    lDataCmb, lData, lCom = GD.GetDataOfRoutine("L5", dtDates, "0")
    lWvl, lCc, lECc, lCcInfo = GD.DoDataCorrectionOneDate(2048, 3, lDataCmb, lData, dtDates)
    print
