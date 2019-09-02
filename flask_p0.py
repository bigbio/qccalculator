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
import os

def handle_plot_widget(ppp):
    s = None
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

def plot_TIC(tic_table):
    d= {'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in tic_table.value['RT']]))),
        'int': robjects.FloatVector(tuple(tic_table.value['int']))   }
    dataf = robjects.DataFrame(d)
    scales = importr('scales')
    c0 = robjects.r('c(0,0)')

    b=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10)]))
    lim_maj=int(max(tic_table.value['RT'])//(60*30))
    lim_min=int(max(tic_table.value['RT'])//(60*10))+1
    b_maj=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*30* i) for i in range(0,lim_maj+1)]))
    b_min=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10* i) for i in range(0,lim_min+1)]))

    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_line() + \
        ggplot2.aes_string(x='RT', y='int') + \
        ggplot2.scale_x_datetime(breaks=b_maj, minor_breaks=b_min, labels = scales.date_format("%H:%M"), expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Intensity", x="Time") + \
        ggplot2.ggtitle("TIC")
    #does not work: date_minor_breaks=scales.date_breaks("5 minutes")
    
    # ltb = robjects.r('theme(plot.margin = unit(c(.1,1,.1,.1), "cm"))')
    # pp = pp + ltb

    plotly = importr('plotly')
    ppp = plotly.ggplotly(pp)
    return handle_plot_widget(ppp)

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
        pltly = plot_TIC(m)
        print("done")
    return jsonify(pltly)

if __name__ == "__main__":
    app.run(debug=True)

#~ python flask_p0.py