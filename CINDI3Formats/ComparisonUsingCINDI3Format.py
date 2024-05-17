from pylab import *
from glob import glob
from os import path as ospath
from linecache import getline
from datetime import timedelta, datetime

from sys import path as syspath
syspath.append("C:/Blick/src/")
syspath.append("C:/Blick/scripts/auxilary/")
from blick_csxfus import blick_csxfus
csxfus = blick_csxfus()
from blick_xfus import blick_xfus
xfus = blick_xfus()
from FunctionsDT import PythonDateNumToBlickDateNum, BlickDateNumToPythonDateNum

def DtConverstion(doy, utc):
    # Ref date
    ## Day
    dtDate = [datetime(2016, 1, 1) + timedelta(days=i + 1 - 1) for i in doy.astype(int)]
    ## Hour
    dtDate = [dtDate[ii] + timedelta(hours=i) for ii, i in enumerate(utc.astype(int))]
    ## Min
    dtDate = [dtDate[ii] + timedelta(minutes=i) for ii, i in enumerate((utc * 60.) % 60.)]
    ## Sec
    dtDate = [dtDate[ii] + timedelta(seconds=i) for ii, i in enumerate((utc * 3600.) % 60.)]
    ## To num
    fDate = date2num(dtDate)
    ## To Blick
    fDateBlick = PythonDateNumToBlickDateNum(fDate)

    return fDate, fDateBlick

def PlotScatSun(d, iPRef, sG, sP):
    a2DatRef = d[iPRef][sG]
    # Ref date
    dtRef, dtBlickRef = DtConverstion(a2DatRef['DOY'], a2DatRef['UTC'])

    h = []
    f, axAll = subplots(5, 2, subplot_kw={'aspect': 'equal'})
    # f.set_size_inches(18 / 2.58, 6. / 2.58)
    axAll = axAll.ravel()

    for iEI, iEC in enumerate(lE):
        ax = axAll[iEI]
        a1Max = array([])
        a1Min = array([])
        for iPI, iPC in enumerate(lPoR):
            a2Dat = d[iPC][sG]
            if any(a2Dat):
                if iEC < 90:
                    bG = (a2Dat['Elev'] >= iEC - 0.5) & (a2Dat['Elev'] < iEC + 0.5)
                    bGRef = (a2DatRef['Elev'] >= iEC - 0.5) & (a2DatRef['Elev'] < iEC + 0.5)
                else:
                    bG = a2Dat['Elev'] > 40
                    bGRef = a2DatRef['Elev'] > 40
                bG &= a2Dat[sP] != -999
                # bG &= (a2Dat[sP] >= -1e20) & (a2Dat[sP] <= 1e20)
                bGRef &= a2DatRef[sP] != -999
                # bGRef &= (a2DatRef[sP] >= -1e20) & (a2DatRef[sP] <= 1e20)

                dt, dtBlick = DtConverstion(d[iPC][sG]['DOY'], d[iPC][sG]['UTC'])
                dtInt, a1Int, a1RefInt, _, _ = \
                    csxfus.timeint(dtBlick[bG], a2Dat[sP][bG], dtBlickRef[bGRef], a2DatRef[sP][bGRef],
                                   x="2", dxmax1=fMaxTimeDif * 60., dxmax2=fMaxTimeDif * 60., meth=1)

                a1Min = append(a1Min, a1Int.min())
                a1Min = append(a1Min, a1RefInt.min())
                a1Max = append(a1Max, a1Int.max())
                a1Max = append(a1Max, a1RefInt.max())

                if iEI == 0:
                    h += ax.plot(a1RefInt, a1Int, 'o', markersize=masi, markeredgecolor='none',
                            markerfacecolor=color[iPI], label=str(iPC))
                else:
                    ax.plot(a1RefInt, a1Int, 'o', markersize=masi, markeredgecolor='none',
                            markerfacecolor=color[iPI])
                ax.set_title('VEA {}'.format(iEC))
            fMin = a1Min.min()
            fMax = a1Max.max()
            ax.set_xlim([fMin, fMax])
            ax.set_ylim([fMin, fMax])
            ax.plot([fMin, fMax], [fMin, fMax], '-', linewidth=1.5, color='k', alpha=0.5, zorder=1)

            ax.locator_params(axis='y', nbins=5)
            ax.locator_params(axis='x', nbins=5)

    f.suptitle('{}_{}'.format(sG, sP), fontsize=12)

    f.subplots_adjust(hspace=0.5)
    f.subplots_adjust(wspace=-0.75)

    legend(handles=h, loc='center left', bbox_to_anchor=(1, 0.5), fancybox = True, framealpha = 0.0)
    # legend(handles=h, ncol=5, loc=2, scatterpoints=1, bbox_to_anchor=(0.125, -0.15, 0.77, .102), mode="expand",
    #        borderaxespad=0.,
    #        bbox_transform=gcf().transFigure, fancybox=True, prop={'size': 8}, framealpha=0.0)
    sFname = '{}_{}.png'.format(sG, sP)
    sFname = sFname.replace('/', '_')
    savefig(sFname, dpi=300., bbox_inches='tight')
    close()


def PlotDirSun(d, iPRef, sG, sP):
    a2DatRef = d[iPRef][sG]
    # Ref date
    dtRef, dtBlickRef = DtConverstion(a2DatRef['DOY'], a2DatRef['UTC'])

    h = []
    f, ax = subplots(1, 1, subplot_kw={'aspect': 'equal'})
    a1Max = array([])
    a1Min = array([])
    for iPI, iPC in enumerate(lPoR):
        a2Dat = d[iPC][sG]
        if any(a2Dat):
            bG = a2Dat[sP] != -999.
            # bG &= a2Dat['RMS'] <= 1.2e-3
            # bG &= a2Dat['RMS'] <= 5.e-3
            bGRef = a2DatRef[sP] != -999.
            # bGRef &= a2DatRef['RMS'] <= 1.2e-3
            # bGRef &= a2DatRef['RMS'] <= 5.e-3

            # a = bG & (a2Dat['RMS'] <= 1.2e-3)
            # b = bG & (a2Dat['RMS'] <= 5.e-3)
            # print sG, iPC, (float(a.sum())/a.shape[0])*100., (float(b.sum())/b.shape[0])*100.

            dt, dtBlick = DtConverstion(d[iPC][sG]['DOY'], d[iPC][sG]['UTC'])
            dtInt, a1Int, a1RefInt, _, _ = \
                csxfus.timeint(dtBlick[bG], a2Dat[sP][bG], dtBlickRef[bGRef], a2DatRef[sP][bGRef],
                               # x="2", dxmax1=-1, dxmax2=-1, meth=1)
                               x="2", dxmax1=fMaxTimeDif * 60., dxmax2=fMaxTimeDif * 60., meth=1)

            a1Min = append(a1Min, a1Int.min())
            a1Min = append(a1Min, a1RefInt.min())
            a1Max = append(a1Max, a1Int.max())
            a1Max = append(a1Max, a1RefInt.max())

            h += ax.plot(a1RefInt, a1Int, 'o', markersize=masi, markeredgecolor=color[iPI],
                    markerfacecolor='none', label='Pan'+str(iPC), alpha=0.6)
            # h += ax.plot(a1RefInt, a1Int, 'o', markersize=masi, markeredgecolor='none',
            #         markerfacecolor=color[iPI], label=str(iPC), alpha=0.6)
        fMin = a1Min.min()
        fMax = a1Max.max()
        ax.set_xlim([fMin, fMax])
        ax.set_ylim([fMin, fMax])
        ax.plot([fMin, fMax], [fMin, fMax], '-', linewidth=1.5, color='k', alpha=0.5, zorder=1)

        ax.set_xlabel('DSCD Pan27')
        ax.set_ylabel('DSCD Pan')

        ax.locator_params(axis='y', nbins=5)
        ax.locator_params(axis='x', nbins=5)

    f.suptitle('{}_{} (RMS <= 1.2e-3)'.format(sG, sP), fontsize=12)

    f.subplots_adjust(hspace=0.5)
    f.subplots_adjust(wspace=-0.75)

    legend(handles=h, loc='center left', bbox_to_anchor=(1, 0.5), fancybox = True, framealpha = 0.0)
    # legend(handles=h, ncol=5, loc=2, scatterpoints=1, bbox_to_anchor=(0.125, -0.15, 0.77, .102), mode="expand",
    #        borderaxespad=0.,
    #        bbox_transform=gcf().transFigure, fancybox=True, prop={'size': 8}, framealpha=0.0)
    sFname = '{}_{}.png'.format(sG, sP)
    sFname = sFname.replace('/', '_')
    savefig(sFname, dpi=300., bbox_inches='tight')
    close()


def PlotDirSun2(d, iPRef, sG, sP):

    h = []
    f, ax = subplots(1, 1)
    for iPI, iPC in enumerate(lP):
        a2Dat = d[iPC][sG]
        if any(a2Dat):
            bG = a2Dat[sP] != -999.
            dt, dtBlick = DtConverstion(d[iPC][sG]['DOY'], d[iPC][sG]['UTC'])
            h += ax.plot(num2date(dt[bG]), a2Dat[sP][bG], 'o', markersize=masi, markeredgecolor=color[iPI],
                    markerfacecolor='none', label='Pan'+str(iPC), alpha=0.6)
        ax.set_xlabel('TIME [UTC]')
        ax.set_ylabel('DSCD Pan')

        # ax.locator_params(axis='y', nbins=5)
        # ax.locator_params(axis='x', nbins=5)

    f.suptitle('{}_{} (RMS <= 1.2e-3)'.format(sG, sP), fontsize=12)

    f.subplots_adjust(hspace=0.5)
    f.subplots_adjust(wspace=-0.75)

    legend(handles=h, loc='center left', bbox_to_anchor=(1, 0.5), fancybox = True, framealpha = 0.0)
    # legend(handles=h, ncol=5, loc=2, scatterpoints=1, bbox_to_anchor=(0.125, -0.15, 0.77, .102), mode="expand",
    #        borderaxespad=0.,
    #        bbox_transform=gcf().transFigure, fancybox=True, prop={'size': 8}, framealpha=0.0)
    sFname = '{}_{}_timeseries.png'.format(sG, sP)
    sFname = sFname.replace('/', '_')
    savefig(sFname, dpi=300., bbox_inches='tight')
    close()



def GetData():

    d = {}
    for iPC in lP:
        d[iPC] = {}
        for sGC in lG:
            d[iPC][sGC] = {}
            a2Full = None
            lH = None
            doInit = True
            for iDC in lD:
                sPthFle = '*_*_{}_{}_*_{}_v{}.asc'.format(iPC, sGC, iDC, iVers)
                lPth = glob(ospath.join(sPthDat, sPthFle))
                del sPthFle
                if len(lPth) == 1:
                    sPth = lPth[0]
                    a2 = loadtxt(sPth, skiprows=dHL[sGC])
                    if doInit:
                        a2Full = array([]).reshape((0, a2.shape[1]))
                        sH = getline(sPth, dHL[sGC])
                        lH = sH[1:].split()
                        doInit = False
                    a2Full = vstack((a2Full, a2))
                else:
                    '{} {} {} skipped'.format(iPC, sGC, iDC)
            if any(a2Full):
                for sH, a1 in zip(lH, a2Full.T):
                    d[iPC][sGC][sH] = a1
                del lH, a2Full

    return d


########################################################################################################################
if __name__ == '__main__':
    sPthDat = 'C:/Blick/data/CINDI2/BlickP_CINDI_2_files_dirsun_Modified_filtered_v4/'
    # sPthDat = 'C:/Blick/data/CINDI2/BlickP_CINDI_2_files_Modified/'
    # sPthDat = 'C:/Blick/data/CINDI2/'
    lD = [20160912, 20160913, 20160914, 20160915, 20160916, 20160917, 20160918, 20160919, 20160920, 20160921, 20160922,
          20160923, 20160924, 20160925, 20160926, 20160927, 20160928]

    # lP = [23, 26, 260, 27, 270, 31, 32]
    # lP = [26, 27]
    lP = [26, 260, 27, 270, 32]

    lE = [1, 2, 3, 4, 5, 6, 8, 15, 30, 90]
    iPRef = 27
    iVers = 4

    # lPoR = [23, 26, 260, 270, 31, 32]
    # lPoR = [26]
    lPoR = [26, 260, 270, 32]

    lG = ['NO2vis', 'NO2visSmall', 'NO2uv', 'HCHO', 'O3vis', 'O3uv']
    dHL = {
    'NO2vis': 33 + 1,
    'NO2visSmall': 33 + 1,
    'NO2uv': 39 + 1,
    'HCHO': 36 + 1,
    'O3vis': 38 + 1,
    'O3uv': 32 + 1,
    }
    # dHL = {
    #     'NO2vis': 33,
    #     'NO2visSmall': 33,
    #     'NO2uv': 39,
    #     'HCHO': 36,
    #     'O3vis': 38,
    #     'O3uv': 32,
    # }

    # fMaxTimeDif = 60.
    fMaxTimeDif = 15.
    masi = 3

    params = {
        'font.size': 6,
    }
    rcParams.update(params)
    colmap = get_cmap('Set1')
    color = [colmap(i) for i in linspace(0, 1, 8)][:-1]

    d = GetData()

    # figure()
    # for i in lP:
    #     plot(d[i]['NO2vis']['DOY'],
    #          d[i]['NO2vis']['NO2_DSCD'],
    #          'o', label='{}'.format(i))
    # legend()

    print
    lProd = [
        'NO2_DSCD',
        'O3_DSCD',
        'O4_DSCD',
        'H2O_DSCD',
        'HCHO_DSCD',
        'BrO_DSCD',
        'NO2_DSCD_220',
        'NO2_DSCD_220', 'NO2_DSCD_298',
        'O3_DSCD_223', 'O3_DSCD_243', 'O3_DSCD_293',
        'O4_DSCD_293',
        'HCHO_DSCD_297',
        'BrO_DSCD_223',
        'H2O_DSCD_296',
        # 'Ring',
        # 'offset_cst',
        # 'offset_lin',
        # 'CI',
        # 'Intens',
    ]

    for sG in lG:
        for sP in lProd:
            lPAll = d[iPRef][sG].keys()
            doPlot = False
            for sPPart in lPAll:
            # for sPPart in lPAll:
                if (sP == sPPart) and ('Error' not in sPPart):
                # if (sP in sPPart) and ('Error' not in sPPart):
                    doPlot = True
                    sP = sPPart
            if doPlot:
            # lPAll = d[iPRef][sG].keys()
            # doPlot = False
            # for sPPart in lPAll:
            #     if (sP == sPPart) and ('Error' not in sPPart):
            #         doPlot = True
            #         sP = sPPart
            # if doPlot:
                print sG, sP
                if iVers <= 2:
                    PlotScatSun(d, iPRef, sG, sP)
                elif iVers == 4:
                    # PlotDirSun2(d, iPRef, sG, sP)
                    PlotDirSun(d, iPRef, sG, sP)