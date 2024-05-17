__author__ = 'Martin Tiefengraber'

from sys import path as syspath
syspath.append("C:/Blick/src")

from pylab import *

import datetime

from blick_io import blick_io
io = blick_io()

from blick_psio import blick_psio
psio = blick_psio()

from blick_params import blick_params
gp = blick_params()

from blick_xfus import blick_xfus
xfus = blick_xfus()

lLevCols = gp.lev0col

date = datetime.date(2016,8,28)
dataPth = ['C:/Blick/data/L0/Pandora209s1_Izana_20160828_L0.txt']

#> Extract routines and dates
rtn = {}
rtnNme = {}
dt = {}
for dataPthI, dataPthC in enumerate(dataPth):
    res, lLMeData, colind, aline = io.read_dataheader(dataPthC, lLevCols)
    assert res == 'OK', res
    res, lA, readinfo = io.read_data(dataPthC, colind, lLevCols, aline, 'ALL')
    assert res == 'OK', res
    key = dataPthC.split('/')[-1].partition('Pandora')[-1].split('_')[0]
    rtn[key] = []
    rtnNme[key] = []
    dt[key] = []
    for rtnCnt in xrange(len(lA)):
        rtnMeas = lA[rtnCnt][0][0]
        if rtnMeas in ['E0']:
            rtnMeas = 'E1'
        if rtnMeas in ['AB']:
            rtnMeas = 'AC'
        c = array([ord(c) for c in rtnMeas]).sum()
        for rtnCntInt in xrange(len(lA[rtnCnt])):
            rtn[key].append(c)
            rtnNme[key].append(rtnMeas)
            dtP = lA[rtnCnt][rtnCntInt][1][0]
            dtP = xfus.zeit2timestr(dtP, 4)
            dt[key].append(datetime.datetime.strptime(dtP, '%Y%m%dT%H%M%SZ'))
            # dt[key].append(datestr2num(dtP))

#> Extract routine names
rtns = unique(array(rtnNme.values()[0]))
rtnCord = zeros(rtns.shape)
for rtnI, rtnC in enumerate(rtns):
    rtnCord[rtnI] = array([ord(c) for c in rtnC]).sum()

#> Plot timing
figure()
for panC, marker, color in zip(rtn.keys(), ['|', '|', '|', '|', '|', '|'], ['b', 'g', 'r', 'm', 'c', 'y']):
    plot_date(dt[panC], rtn[panC], markersize=15, marker=marker, markerfacecolor='none',  markeredgecolor=color,
              markeredgewidth=4, linestyle='none', label=panC)
#> Plot reference times
axvline(datetime.datetime.combine(date, datetime.time(4,0,0)), color='r')
for hour in [6, 7, 8, 9, 10, 12, 13, 14, 15, 16]:
    for minute in [0, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]:
        axvline(datetime.datetime.combine(date, datetime.time(hour, minute, 0)), color='r')
for hour in [11]:
    for minute in [0, 10, 15, 25, 30, 41, 45, 55]:
        axvline(datetime.datetime.combine(date, datetime.time(hour, minute, 0)), color='r')

xlim([min(dt.values()[0]), max(dt.values()[0])])
yticks(rtnCord, rtns)
legend()
grid(linestyle='-', color='gray')
show()