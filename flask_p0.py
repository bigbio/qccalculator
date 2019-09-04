from flask import Flask, request, render_template, jsonify
from MZQC import MZQCFile as mzqc
from typing import List, Dict, Set, Any, Optional, Callable
from collections import defaultdict
from statistics import mean, median, stdev
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
import json
import math, datetime
import tempfile
import base64
import re
import numpy as np
from enum import Enum
import os

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
                htmlwidgets.saveWidget(ppp, t.name, libdir="replaceme", selfcontained = False)
                # start stupid fix to get all the recent libs written in the flask lib directory
                htmlwidgets.saveWidget(ppp, 'bof', libdir='lib', selfcontained = False)
                os.remove('bof')
                # end stupid fix
                with open(t.name, "r") as f:
                    s = f.read()
                    s = s.replace("replaceme", "http://localhost:5000/lib")
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

# https://stackoverflow.com/questions/20646822/how-to-serve-static-files-in-flask 
app = Flask(__name__, template_folder='serving_static', static_url_path='/lib', static_folder='lib')

@app.route('/')
def root():
    message = "Hello there!"
    return render_template('p0.html', message=message)

@app.route('/echo', methods=['POST'])
def hello():
    if request.is_json:
        print("start")
        content = request.json
        m = mzqc.JsonSerialisable.FromJson(json.dumps(content))
        cmpltn = datetime.datetime.strptime("2012-02-03 11:00:41", '%Y-%m-%d %H:%M:%S')
        mxrt = datetime.timedelta(seconds=m.value['RT'][-1])
        pltly = plot_TIC(m, cmpltn-mxrt , PlotType.PLOTLY)
        print("done")
    return jsonify(pltly)

if __name__ == "__main__":
    app.run(debug=True)

#~ python flask_p0.py