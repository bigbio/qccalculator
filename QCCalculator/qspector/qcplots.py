import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as robjects
from rpy2.robjects import Formula, Environment
from rpy2.robjects.vectors import IntVector, FloatVector
from rpy2.robjects.lib import grid
from rpy2.robjects.packages import importr, data
from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.robjects.lib.dplyr import DataFrame
from rpy2.robjects.lib.dplyr import (mutate, group_by, summarize)
from rpy2.rinterface import parse
import warnings
import math, datetime
import tempfile
import base64
import re
import numpy as np
from enum import Enum
import os
from statistics import mean, median, stdev


class PlotType(Enum):
    PNG = 1
    SVG = 2
    PLOTLY = 3

def handle_plot_format(pp, plot_type: PlotType):
    if plot_type == PlotType.PLOTLY:
        plotly = importr('plotly')
        ppp = plotly.ggplotly(pp)
        htmlwidgets = importr('htmlwidgets')
        with tempfile.NamedTemporaryFile() as t:
                htmlwidgets.saveWidget(ppp, t.name, libdir='lib', selfcontained = False)
                # start stupid fix to get all the recent libs written in the flask lib directory
                htmlwidgets.saveWidget(ppp, 'bof', libdir='lib', selfcontained = False)
                os.remove('bof')
                # end stupid fix
                with open(t.name, "r") as f:
                    s = f.read()
        return s
    else:
        with tempfile.NamedTemporaryFile() as t:
            if plot_type == PlotType.SVG:
                grdevices = importr('grDevices')
                grdevices.svg(file=t.name)
                pp.plot()
                grdevices.dev_off()

                with open(t.name, "r") as f:
                    s = f.read()
            else:
                grdevices = importr('grDevices')
                grdevices.png(file=t.name, width=512, height=512)
                pp.plot()
                grdevices.dev_off()

                with open(t.name, "rb") as fb:
                    s = base64.b64encode(fb.read()).decode()
        return s

def plot_TIC(tic_table, start_time, plot_type=PlotType.PNG):
    d= {'RT': robjects.POSIXct((tuple([start_time + datetime.timedelta(seconds=i) for i in tic_table.value['RT']]))),
        'int': robjects.FloatVector(tuple(tic_table.value['int']))   }
    dataf = robjects.DataFrame(d)
    scales = importr('scales')
    c0 = robjects.r('c(0,0)')

    lim_maj=int(max(tic_table.value['RT'])//(60*30))
    lim_min=int(max(tic_table.value['RT'])//(60*10))+1
    b_maj=robjects.POSIXct(tuple([start_time + datetime.timedelta(seconds=60*30* i) for i in range(0,lim_maj+1)]))
    b_min=robjects.POSIXct(tuple([start_time + datetime.timedelta(seconds=60*10* i) for i in range(0,lim_min+1)]))

    axislabels = robjects.StrVector(tuple([(datetime.datetime.fromtimestamp(60*30* i)).strftime("%H:%M") for i in range(0,lim_maj+1)]))

    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_line() + \
        ggplot2.aes_string(x='RT', y='int') + \
        ggplot2.scale_x_datetime(breaks=b_maj, minor_breaks=b_min, labels = axislabels, expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Intensity", x="Time") + \
        ggplot2.ggtitle("TIC")
    #does not work: date_minor_breaks=scales.date_breaks("5 minutes")
    # scales.date_format("%H:%M")

    # ltb = robjects.r('theme(plot.margin = unit(c(.1,1,.1,.1), "cm"))')
    # pp = pp + ltb

    return handle_plot_format(pp, plot_type)

def plot_SN(table, mslevel=2, svg_plot=False):
    d= {'SN': robjects.FloatVector(tuple(table.value['SN'])) }
    dataf = robjects.DataFrame(d)
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    m = median(table.value['SN'])
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(binwidth=.5, colour="black", fill="white") + \
        ggplot2.aes_string(x='SN', y='..density..') + \
        ggplot2.geom_density(alpha=.1, fill="green") + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='median(SN, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'median={}'".format(str(round(m,2))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(x="S/N") + \
        ggplot2.ggtitle("S/N distribution in MS{} spectra".format(str(mslevel)))
    
    return handle_plot_format(pp, svg_plot)

def plot_dppm(psm_table, svg_plot=False):
    d= {'deltaPPM': robjects.FloatVector(tuple(psm_table.value['delta_ppm']))   }
    dataf = robjects.DataFrame(d)
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')
    m = median(psm_table.value['delta_ppm'])
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(binwidth=.5, colour="black", fill="white") + \
        ggplot2.aes_string(x='deltaPPM', y='..density..') + \
        ggplot2.geom_density(alpha=.1, fill="green") + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='median(deltaPPM, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'median={}'".format(str(round(m,2))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(expand=c0, limit=c10) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(x=parse('paste(Delta, "ppm")'), y="Frequency density")  + \
        ggplot2.ggtitle("Mass error distribution")
    
    return handle_plot_format(pp, svg_plot)

def plot_lengths(psm_table, svg_plot=False):
    regex_mod = r'(\([^\(]*\))'
    regex_noaa = r'([^A-Za-z])'
    d= {'PeptideSequence': robjects.StrVector(tuple(psm_table.value['peptide_sequence'])),
        'Length': robjects.IntVector(tuple([len(re.sub(regex_noaa, '', re.sub(regex_mod, '', x))) for x in psm_table.value['peptide_sequence']])) }
    dataf = robjects.DataFrame(d)
    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')
    m = mean(d['Length'])
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(binwidth=1, origin = -0.5, colour="black", fill="white") + \
        ggplot2.aes_string(x='Length') + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='mean(Length, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'mean={}'".format(str(round(m,2))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Count")  + \
        ggplot2.ggtitle("Length distribution of identified peptide sequences")
        # parse('paste(Delta, "ppm")' does not work in ggplot2.ggtitle

    return handle_plot_format(pp, svg_plot)

def plot_topn(prec_table, surv_table, svg_plot=False):
    h = np.histogram(prec_table.value['RT'], bins=surv_table.value['RT']+[surv_table.value['RT'][-1]+surv_table.value['RT'][-2]])

    d= {'SN': robjects.FloatVector(tuple(surv_table.value['SN'])), 
        'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in surv_table.value['RT']]))), 
        'TopN': robjects.IntVector(tuple(h[0])) }
    dataf = robjects.DataFrame(d)

    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')
    m = median(d['TopN'])
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(binwidth=1, origin = -0.5, colour="black", fill="white") + \
        ggplot2.aes_string(x='TopN') + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='median(TopN, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'median={}'".format(str(round(m))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(breaks=scales.pretty_breaks(), expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Count")  + \
        ggplot2.ggtitle("TopN sampling")

    return handle_plot_format(pp, svg_plot)

def plot_topn_sn(prec_table, surv_table, svg_plot=False):
    h = np.histogram(prec_table.value['RT'], bins=surv_table.value['RT']+[surv_table.value['RT'][-1]+surv_table.value['RT'][-2]])
    qs =  np.quantile(surv_table.value['SN'], [.25,.5,.75])
    d= {'SN': robjects.FloatVector(tuple(surv_table.value['SN'])), 
        'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in surv_table.value['RT']]))), 
        'TopN': robjects.IntVector(tuple(h[0])),
        'SN.quartile': robjects.FactorVector((tuple( ['q1' if v < qs[0] else 'q2' if qs[0]<v<qs[1] else 'q3' if qs[1]<v<qs[2] else 'q4' for v in surv_table.value['SN']] )))}
    dataf = robjects.DataFrame(d)

    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')
    m = median(d['TopN'])
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.aes_string(x='TopN', fill='SN.quartile') + \
        ggplot2.geom_histogram(binwidth=1, origin = -0.5) + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='median(TopN, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'median={}'".format(str(round(m))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(breaks=scales.pretty_breaks(), expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Count")  + \
        ggplot2.ggtitle("TopN sampling")
        # parse('paste(Delta, "ppm")' does not work in ggplot2.ggtitle

    return handle_plot_format(pp, svg_plot)

def plot_topn_rt(prec_table, surv_table, svg_plot=False):
    h = np.histogram(prec_table.value['RT'], bins=surv_table.value['RT']+[surv_table.value['RT'][-1]+surv_table.value['RT'][-2]])

    d= {'SN': robjects.FloatVector(tuple(surv_table.value['SN'])), 
        'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in surv_table.value['RT']]))), 
        'TopN': robjects.IntVector(tuple(h[0])) }
    dataf = robjects.DataFrame(d)

    b=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10)]))
    lim_maj=int(max(surv_table.value['RT'])//(60*30))
    lim_min=int(max(surv_table.value['RT'])//(60*10))+1
    b_maj=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*30* i) for i in range(0,lim_maj+1)]))
    b_min=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10* i) for i in range(0,lim_min+1)]))

    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')

    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_col(width = 1) + \
        ggplot2.aes_string(x='RT',y='TopN') + \
        ggplot2.scale_x_datetime(breaks=b_maj, minor_breaks=b_min, labels = scales.date_format("%H:%M"), expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Count", x="TopN sampling range")  + \
        ggplot2.ggtitle("TopN utilisation")


    # TODO also plot real histogram with color of SN or target/decoy?
    return handle_plot_format(pp, svg_plot)

def plot_idmap(prec_table, psm_table, svg_plot=False):
    d_psm= {'MZ': robjects.FloatVector(tuple(psm_table.value['MZ'])), 
        'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in psm_table.value['RT']]))), 
        'col': robjects.FactorVector(tuple(["identified"]*len(psm_table.value['MZ']))) }
    dataf_psm = robjects.DataFrame(d_psm)

    d_prc= {'MZ': robjects.FloatVector(tuple(prec_table.value['MZ'])), 
        'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in prec_table.value['RT']]))), 
        'col': robjects.FactorVector(tuple(["recorded"]*len(prec_table.value['MZ']))) }
    dataf_prc = robjects.DataFrame(d_prc)

    b=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10)]))
    lim_maj=int(max(prec_table.value['RT'])//(60*30))
    lim_min=int(max(prec_table.value['RT'])//(60*10))+1
    b_maj=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*30* i) for i in range(0,lim_maj+1)]))
    b_min=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10* i) for i in range(0,lim_min+1)]))

    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')

    pp = ggplot2.ggplot(dataf_prc) + \
        ggplot2.geom_point() + \
        ggplot2.aes_string(x='RT',y='MZ', colour='col') + \
        ggplot2.scale_x_datetime(breaks=b_maj, minor_breaks=b_min, labels = scales.date_format("%H:%M"), expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="m/z", x="Time")  + \
        ggplot2.geom_point(data=dataf_psm) + \
        ggplot2.ggtitle("ID map")

    ltb = robjects.r('theme(legend.title=element_blank())')
    pp = pp + ltb

    # TODO also plot real histogram with color of SN or target/decoy?
    return handle_plot_format(pp, svg_plot)

def plot_charge(prec_table, psm_table=None, svg_plot=False):
    d_prc= {'c': robjects.IntVector(tuple(prec_table.value['c'])), 
        'col': robjects.FactorVector(tuple(["recorded"]*len(prec_table.value['c']))) }
    dataf = robjects.DataFrame(d_prc)

    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')

    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_bar() + \
        ggplot2.aes_string(x='c',fill='col') + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.scale_x_continuous(breaks=scales.pretty_breaks(), expand=c0) + \
        ggplot2.labs(x="Charge state", y="Count")  + \
        ggplot2.ggtitle("Charge states") + \
        robjects.r('theme(legend.title=element_blank())')

    # TODO N/A handling???

    if psm_table:
        d_prc= {'c': robjects.IntVector(tuple(psm_table.value['c'])), 
            'col': robjects.FactorVector(tuple(["identified"]*len(psm_table.value['c']))) }
        dataf_id = robjects.DataFrame(d_prc)
        pp = pp + ggplot2.geom_bar(data=dataf_id)

    # TODO also plot real histogram with color of SN or target/decoy?
    return handle_plot_format(pp, svg_plot)

def plot_peaknum(table, mslevel=2, svg_plot=False):
    d= {'peakcount': robjects.IntVector(tuple(table.value['peakcount']))}
    dataf = robjects.DataFrame(d)
    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')
    m = mean(d['peakcount'])
    binw = round(m/15)+1  #+1 avoids 0 binwidth
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(binwidth=binw, origin = -0.5, colour="black", fill="white") + \
        ggplot2.aes_string(x='peakcount') + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='mean(peakcount, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'mean={}'".format(str(round(m,2))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Count")  + \
        ggplot2.ggtitle("Peak count distribution for MS{} spectra".format(str(mslevel)))

    return handle_plot_format(pp, svg_plot)

def plot_intensities(table, svg_plot=False):
    grdevices = importr('grDevices')
    d= {'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in table.value['RT']]))),
        'int': robjects.FloatVector(tuple(table.value['int']))}
    dataf = robjects.DataFrame(d)
    scales = importr('scales')
    cow = importr('cowplot')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    combiaxis1 = robjects.r('theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin=unit(c(1,1,-0.1,1), "cm"))')
    combiaxis2 = robjects.r('theme(axis.text.y = element_blank(), axis.ticks = element_blank(), plot.margin=unit(c(-0.1,1,1,1), "cm"))')

    ht = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(bins=100) + \
        ggplot2.aes_string(x='int') + \
        combiaxis1 + \
        ggplot2.labs(y="Count", x= robjects.r('element_blank()') )  + \
        ggplot2.ggtitle("Intensity distribution")

    bx = ggplot2.ggplot(dataf) + \
        ggplot2.geom_boxplot() + \
        ggplot2.aes_string(x=1, y='int') + \
        ggplot2.coord_flip() + \
        combiaxis2 + \
        ggplot2.labs(x="", y="Intensities")
        
    # pracma = importr('pracma')
    # AUC <- trapz(QCTIC$MS.1000894_.sec.,QCTIC$MS.1000285)
    # auc = np.trapz(table.value['int'], x=table.value['RT'])

    # Qs <- quantile(QCTIC$MS.1000285,prob = c(0.25, 0.5, 0.75))
    # Qs <- data.frame(Qs)
    # qs =  np.quantile(table.value['int'], [.25,.5,.75])

    pp = cow.plot_grid(ht,bx, ncol = 1, align = 'v', axis = 'l', rel_heights = robjects.r('c(1,.25)'))


    # grdevices.png(file="tests/grid_png_func.png", width=512, height=512)
    # c = cow.plot_grid(ht,bx, ncol = 1, align = 'v', axis = 'l', rel_heights = robjects.r('c(1,.25)'))
    # c.plot()
    # grdevices.dev_off()

    return handle_plot_format(pp, svg_plot)

def plot_events(tic_table, surv_table, prec_table, psm_table=None, svg_plot=False):
    datasources = [("Chromatogram", tic_table),("MS1",surv_table),("MS2", prec_table)]
    if psm_table:
        datasources.append(("Identifications",psm_table))
    
    dataf = robjects.DataFrame({})
    for annot, tab in datasources:
        if annot == "Chromatogram":
            #quartiles of total intensity recorded
            tic = np.sum(tab.value['int'])
            vec = (np.cumsum(tab.value['int'])/tic)
            qs =  [.25,.5,.75]
        else:
            vec = [i[0]+1 for i in enumerate(tab.value['RT'])]
            qs =  np.quantile(vec, [.25,.5,.75])

        d_c= {'Time': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in tab.value['RT']]))),
            'Quartile': robjects.FactorVector((tuple( ['q1' if v < qs[0] else 'q2' if qs[0]<v<qs[1] else 'q3' if qs[1]<v<qs[2] else 'q4' for v in vec] )))
        }
        td = (DataFrame(d_c).
                group_by('Quartile').
                summarize(n='n()', mx='max(Time)', mi='min(Time)').
                mutate(group="'"+annot+"'", dt='hms::as.hms(mx-mi)') )
        dataf = dataf.rbind(td)

    dataf = ( DataFrame(dataf).arrange('Quartile') )
    scales = importr('scales')
    chron = importr('chron')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')

    # TODO time axis fix and Hz annotations
    b=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10)]))
    lim_maj=int(max(tic_table.value['RT'])//(60*30))
    lim_min=int(max(tic_table.value['RT'])//(60*10))+1
    b_maj=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*30* i) for i in range(0,lim_maj+1)]))
    b_min=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10* i) for i in range(0,lim_min+1)]))

    # TODO shitty time scale still does not work with stacked bars
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_bar(stat='identity', position = ggplot2.position_fill(reverse = True)) + \
        ggplot2.aes_string(x='group',y='mx',fill='Quartile') + \
        ggplot2.coord_flip() + \
        chron.scale_y_chron(format ="%H:%M:%S") + \
        ggplot2.labs(x="Events", y="Time")  + \
        ggplot2.ggtitle("Quartiles of Chromatographic, MS1, MS2, and identification events over RT") 

    dataf.to_csvfile('tests/events.csv')
    return handle_plot_format(pp, svg_plot)
        # ggplot2.scale_y_datetime(breaks=b_maj, minor_breaks=b_min,labels = scales.date_format("%H:%M"), expand=c0) + \

