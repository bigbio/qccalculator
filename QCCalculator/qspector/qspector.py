from flask import Flask, request, render_template, jsonify
from MZQC import MZQCFile as mzqc
import datetime
import json
import os

from QCCalculator import qcplots

PORT = 5000

# https://stackoverflow.com/questions/20646822/how-to-serve-static-files-in-flask 
app = Flask(__name__, template_folder='serving_static', static_url_path='/lib', static_folder='lib')
#TODO this needs parameterisation in-line with handle_plot_format

@app.route('/')
def root():
    message = "Hello there!"
    return render_template('qspector.html', message=message)

@app.route('/echo', methods=['POST'])
def plotly_echo():
    if request.is_json:
        print("start")
        content = request.json
        m = mzqc.JsonSerialisable.FromJson(json.dumps(content))
        #TODO this needs parameterisation and be fixed in the first place in qccalc
        cmpltn = datetime.datetime.strptime("2012-02-03 11:00:41", '%Y-%m-%d %H:%M:%S')
        mxrt = datetime.timedelta(seconds=m.value['RT'][-1])
        pltly = qcplots.plot_TIC(m, cmpltn-mxrt , qcplots.PlotType.PLOTLY, request.host_url.__str__().rsplit(':', 1)[0], PORT)
        print("done")
    return jsonify(pltly)

if __name__ == "__main__":
    app.run(debug=True, port=PORT)

#~ python qspector.py