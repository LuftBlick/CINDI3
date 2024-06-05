import pandas as pd
import xarray as xr
import numpy as np
import os
import h5py
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import seaborn as sns

from params.InputParams import LoadParams
from CINDI3SemiBlind import GetPths
import CINDI3SemiBlind as CSB


def process_trace_gas_file(file_path):
    # Read the file content
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Separate the general header, column descriptions, and data
    general_header = []
    column_descriptions = []
    data_start = 0

    # Identify the start of the column header and data based on lines that start with "%"
    data_format_found = False
    for i, line in enumerate(lines):
        if not line.startswith('%'):
            data_start = i
            break
        elif line.startswith('% Data format:'):
            data_format_found = True
        elif data_format_found and line.startswith('% ------------'):
            continue
        elif data_format_found:
            column_descriptions.append(line.strip())
        else:
            general_header.append(line.strip())

    # Parse the general header into attributes
    general_attrs = {}
    for line in general_header:
        line = line.lstrip('%').strip()
        if ':' in line:
            key, value = line.split(':', 1)
            general_attrs[key.strip()] = value.strip()

    # Create a dictionary to hold column attributes and column names
    column_attrs = {}
    column_names = []
    for description in column_descriptions:
        if ':' in description:
            parts = description.split(':', 2)
            col_num = parts[0].strip()
            col_name = parts[1].strip().split()[0]
            long_name = parts[1].strip()+ ' ' + parts[2].strip()
            column_attrs[col_name] = {'long_name': long_name}
            column_names.append(col_name)

    # Read the data into a DataFrame
    data = pd.read_csv(file_path, delim_whitespace=True, skiprows=data_start, names=column_names)

    # Convert DOY to a timestamp
    def doy_to_timestamp(doy, year=2024):
        return datetime(year, 1, 1) + timedelta(days=doy - 1)

    data['time'] = data['DOY'].apply(doy_to_timestamp)
    data = data.set_index('time')

    # Convert DataFrame to xarray Dataset
    dataset = xr.Dataset.from_dataframe(data)

    # Assign attributes to each data variable based on the column descriptions
    for col in column_names:
        if col in column_attrs:
            dataset[col].attrs = column_attrs[col]

    # Assign general attributes to the dataset
    dataset.attrs = general_attrs

    return dataset


def get_data(par):

    data = {}
    dImp = {}
    dImp['loc'] = par['sLoc']
    dImp['sCode'] = par['dSCode']['s'][0]
    for sDate in par['iDate']:
        dImp['date'] = sDate
        data[dImp['date']] = {}
        # > Loop Pandoras
        for sPanC in par['dPan']:
            # > Loop spectrometers
            for iSpec in par['dPan'][sPanC]:
                dImp['panName'] = 'Pandora{}s{}'.format(sPanC, iSpec)
                data[dImp['date']][dImp['panName']] = {}
                dImp['institute'], dImp['instNumber'] = par['dPanIdInst'][dImp['panName'][7:]]
                # > Loop products
                lProds = [prod for prod, do in par['dProdAna'].items() if do]
                for sProd in lProds:
                    dImp['prod'] = sProd
                    data[dImp['date']][dImp['panName']][dImp['prod']] = {}
                    dImp['fCode'] = par['dProdFCode'][sProd]
                    for sFCode in dImp['fCode']:
                        # > Loop reference types:
                        for sRefType in par['dProdRefT'][sProd]:
                            dImp['refType'] = sRefType
                            # load f-code details
                            dImp['fitPars'] = {}
                            with h5py.File(par['sPFPth'], 'r+') as pssetup:
                                idxFCode = np.where(pssetup['f_codes'][:, 0]['f-code'] == sFCode)[0][0]
                                FCode = pssetup['f_codes'][idxFCode]
                                for fitPar in ['npol', 'noffs', 'nwlc', 'nresc', 'WL-starts', 'WL-ends', 'Ring']:
                                    dImp['fitPars'][fitPar] = FCode[fitPar][0]
                                for fitPar in ['Fitted gases', 'Gas sources', 'Gas temps']:
                                    dImp['fitPars'][fitPar] = np.asarray(FCode[fitPar][0].split(','))
                            if FCode['Reference'][0].upper().startswith(sRefType.upper()):
                                dImp['prodVers'] = par['dProdVers'][sProd][0]
                                # get name of Cindi format file
                                sPthCindi, _, _ = GetPths(dImp['date'], dImp['panName'], dImp['loc'], dImp['institute'],
                                                          dImp['instNumber'], dImp['prod'], CSB.refTypeSyn[dImp['refType']],
                                                          dImp['prodVers'], dImp['sCode'], sFCode)
                                # create CINDI-3 file
                                data[dImp['date']][dImp['panName']][dImp['prod']][dImp['refType']] = (
                                    process_trace_gas_file(os.path.join(par['sCampPth'], sPthCindi)))

    return data


def compare_data_timeseries(data, pltpar):

    days = []
    pans = []
    prods = []
    refs = []
    for day in data:
        days.append(day)
        for pan in data[day]:
            pans.append(pan)
            for prod in data[day][pan]:
                prods.append(prod)
                for ref in data[day][pan][prod]:
                    refs.append(ref)
    days = list(np.unique(np.array(days)))
    pans = list(np.unique(np.array(pans)))
    prods = list(np.unique(np.array(prods)))
    refs = list(np.unique(np.array(refs)))

    for day in days:
        for prod in prods:
            for ref in refs:
                if ref in data[day][pans[0]][prod].keys():
                    di = data[day][pans[0]][prod]['Ref']
                    data_vars = [dv for dv in di.data_vars.keys() if dv not in pltpar['skip']]
                    ndv = np.int(np.ceil(float(len(data_vars)) / float(pltpar['cols'])))

                    f, ax = plt.subplots(ndv, pltpar['cols'], sharex='all',
                                         figsize=(pltpar['fpartwidth']*pltpar['cols'], ndv*pltpar['fpartheight']))
                    ax = np.ravel(ax)
                    for pan in pans:
                        di = data[day][pan][prod][ref]
                        # time restriction
                        di = di.isel(time=(di.time.dt.hour >= pltpar['h_start']) & (di.time.dt.hour < pltpar['h_end']))
                        # exclude negative VEA
                        di = di.isel(time=di.VEA >= 0.)
                        for axi, data_var in zip(ax, data_vars):
                            if data_var in pltpar['dolog']:
                                scale = 'log'
                            else:
                                scale = 'linear'
                            di[data_var].plot(ax=axi, label=pan, yscale=scale,
                                              linestyle='none', marker='o', markersize=pltpar['masi'])
                            axi.set_title(data_var)
                            axi.grid()

                    ax[0].legend(loc=2, bbox_to_anchor=(0,1.3), frameon=0, ncol=3)
                    plt.suptitle('{}, {}, {}'.format(day, prod, ref), fontsize=pltpar['fssupt'], y=1.05)
                    plt.tight_layout()

                    figname = 'TimeSeries_{}_{}(v{})_{}.png'.format(day, prod, pltpar['vers'][prod][0], ref)
                    plt.savefig(os.path.join(pltpar['plotpath'], figname), dpi=pltpar['dpi'], bbox_inches='tight')
                    plt.close()


def compare_data_scatter(data, pltpar):

    days = []
    pans = []
    prods = []
    refs = []
    for day in data:
        days.append(day)
        for pan in data[day]:
            pans.append(pan)
            for prod in data[day][pan]:
                prods.append(prod)
                for ref in data[day][pan][prod]:
                    refs.append(ref)
    days = list(np.unique(np.array(days)))
    pans = list(np.unique(np.array(pans)))
    prods = list(np.unique(np.array(prods)))
    refs = list(np.unique(np.array(refs)))

    for day in days:
        for prod in prods:
            for ref in refs:
                if ref in data[day][pans[0]][prod].keys():
                    di = data[day][pans[0]][prod]['Ref']
                    data_vars = [dv for dv in di.data_vars.keys() if dv not in pltpar['skipscatter']]
                    if not pltpar['showErrors']:
                        data_vars = [dv for dv in data_vars if "Error" not in dv]
                    ndv = np.int(np.ceil(float(len(data_vars)) / float(pltpar['colsscatter'])))

                    # f, ax = plt.subplots(ndv, pltpar['colsscatter'], constrained_layout=True)
                    f, ax = plt.subplots(ndv, pltpar['colsscatter'], constrained_layout=True,
                                         figsize=(pltpar['fpartwidthscatter']*pltpar['colsscatter'], ndv*pltpar['fpartheightscatter']))
                    ax = np.ravel(ax)
                    panref = 'Pandora83s1'
                    # prepare reference:
                    diref = data[day][panref][prod][ref]
                    # time restriction
                    diref = diref.isel(
                        time=(diref.time.dt.hour >= pltpar['h_start']) & (diref.time.dt.hour < pltpar['h_end']))
                    # exclude negative VEA
                    diref = diref.isel(time=diref.VEA >= 0.)

                    for i, pan in enumerate(pans):
                        if pan != panref:
                            di = data[day][pan][prod][ref]
                            # time restriction
                            di = di.isel(time=(di.time.dt.hour >= pltpar['h_start']) & (di.time.dt.hour < pltpar['h_end']))
                            # exclude negative VEA
                            di = di.isel(time=di.VEA >= 0.)
                            # interpolate to reference
                            diint = di.interp(time=diref.time)

                            for axi, data_var in zip(ax, data_vars):
                                sns.scatterplot(x=diref[data_var], y=diint[data_var], ax=axi, label=pan, size=diref['VEA'])
                                # axi.scatter(diref[data_var], diint[data_var], label=pan, s=diref['VEA'])
                                # axi.grid(True, zorder=0)
                                axi.set_aspect('equal', 'box')

                                max_tick = max(axi.get_xticks().max(), axi.get_yticks().max())
                                min_tick = min(axi.get_xticks().min(), axi.get_yticks().min())
                                ticks = np.linspace(min_tick, max_tick, num=len(axi.get_xticks()))

                                if i == 0:
                                    axi.plot(ticks, ticks, color='gray', zorder=0)

                                axi.set_xticks(ticks)
                                axi.set_yticks(ticks)

                                axi.set_xlim(min_tick, max_tick)
                                axi.set_ylim(min_tick, max_tick)

                                # Rotate the x-tick labels by 90 degrees
                                axi.tick_params(axis='x', rotation=90)

                                axi.set_title(diref[data_var].attrs['long_name'].split(',')[0])
                                axi.set_xlabel(panref)
                                axi.set_ylabel('Pandora')

                                # Set background to light gray
                                axi.set_facecolor('lightgray')

                                # Remove the legend from each subplot
                                axi.legend_.remove()

                    # Hide unused subplots
                    for axi in ax[len(data_vars):]:
                        f.delaxes(axi)

                    # Move legend outside the subplots
                    handles, labels = ax[0].get_legend_handles_labels()
                    f.legend(handles, labels, loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False,
                             ncol=1, fontsize=pltpar.get('legend_fontsize', 10))
                    # f.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.15), frameon=False,
                    #          ncol=len(pans)-1, fontsize=pltpar.get('legend_fontsize', 10))
                    # f.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.15), frameon=False,
                    #          ncol=len(pans)-1, fontsize=pltpar.get('legend_fontsize', 10))
                    plt.suptitle('{}, {}, {}'.format(day, prod, ref), fontsize=pltpar['fssupt'], y=1.05)

                    figname = 'Scatter_{}_{}(v{})_{}.png'.format(day, prod, pltpar['vers'][prod][0], ref)
                    plt.savefig(os.path.join(pltpar['plotpath'], figname), dpi=pltpar['dpi'], bbox_inches='tight')
                    plt.close()


if __name__ == '__main__':

    par = {}
    par['delim'] = ' -> '
    par['epoch'] = datetime(2000,1,1)
    #> Read input file
    par['paramsFile'] = 'processCompDataInput.txt'
    par = LoadParams(par)


    data = get_data(par)

    pltpar = {
        'masi': 3,
        'cols': 5,
        'colsscatter': 3,
        'fpartwidth': 4,
        'fpartwidthscatter': 3,
        'fpartheight': 2,
        'fpartheightscatter': 3,
        'dpi': 300.,
        'plotpath': par['sPlotPth'],
        'h_start': 5,
        'h_end': 19,
        'fssupt': 15,
        'vers': par['dProdVers'],
        'dolog': ['RMS'],
        'legend_fontsize': 10,
        'showErrors': False,
        'skip': ['DOY', 'UTC', 'VEA', 'VAA', 'INORM_280', 'INORM_290', 'INORM_300', 'INORM_310', 'INORM_320', 'INORM_330', 'INORM_340', 'INORM_350', 'INORM_360', 'INORM_370', 'INORM_380', 'INORM_390', 'INORM_400', 'INORM_410', 'INORM_420', 'INORM_430', 'INORM_440', 'INORM_450', 'INORM_460', 'INORM_470', 'INORM_480', 'INORM_490', 'INORM_500', 'INORM_510', 'INORM_520', 'INORM_530'],
        'skipscatter': ['DOY', 'UTC', 'ACQT', 'SZA', 'SAA', 'RMS', 'VEA', 'VAA', 'NO2_DSCD_220', 'NO2_DSCD_220_Error', 'O3_DSCD_243', 'O3_DSCD_243_Error', 'O4_DSCD_223', 'O4_DSCD_223_Error', 'INORM', 'INORM_280', 'INORM_290', 'INORM_300', 'INORM_310', 'INORM_320', 'INORM_330', 'INORM_340', 'INORM_350', 'INORM_360', 'INORM_370', 'INORM_380', 'INORM_390', 'INORM_400', 'INORM_410', 'INORM_420', 'INORM_430', 'INORM_440', 'INORM_450', 'INORM_460', 'INORM_470', 'INORM_480', 'INORM_490', 'INORM_500', 'INORM_510', 'INORM_520', 'INORM_530'],
    }

    compare_data_timeseries(data, pltpar)
    compare_data_scatter(data, pltpar)
