import pytest
from MZQC import MZQCFile as mzqc
from QCCalculator.qcplots import *

import tempfile
import datetime
import base64
import re

with open("tests/plot_test_files/test.mzqc", "r") as f:
    tqc = mzqc.JsonSerialisable.FromJson(f.read())['mzQC']
rq = tqc.runQualities[0]

# TODO issues
# * mzqc.JsonSerialisable.FromJson(f.read())['mzQC']
# * inf.fileProperties are quality metrics, not cv params
# * runQualities is list of SetQuality

cmpltn = rq.metadata.inputFiles[0].fileProperties[0].value
maxrt = max(list(filter(lambda x: x.name=="Total ion current chromatogram", rq.qualityMetrics))[0].value['RT'])
strt:datetime.datetime = datetime.datetime.strptime(cmpltn, '%Y-%m-%d %H:%M:%S') - datetime.timedelta(seconds=maxrt)

class TestSvgFunctionality:
    def test_tic_svg(self):
        metric = list(filter(lambda x: x.name=="Total ion current chromatogram", rq.qualityMetrics))
        assert len(metric) == 1
        svg = plot_TIC(metric[0],strt, PlotType.SVG)
        # with open("tests/plot_test_files/tic_svg_func.svg","w") as f:
        #     f.write(svg)
        with open("tests/plot_test_files/tic_svg_func.svg","r") as f:
            cmp_r = f.read()
        assert cmp_r == svg


class TestPlotFunctionality:
    # from here on out the plot code is tested with the PNG format only
    def test_tic_png(self):
        metric = list(filter(lambda x: x.name=="Total ion current chromatogram", rq.qualityMetrics))
        assert len(metric) == 1
        tic = plot_TIC(metric[0],strt)
        # with open("tests/plot_test_files/_tic_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(tic.encode()))
        with open("tests/plot_test_files/tic_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == tic

    def test_sn_png(self):
        metric = list(filter(lambda x: x.name=="Spectrum acquisition metric values - MS2", rq.qualityMetrics))
        assert len(metric) == 1

        sn = plot_SN(metric[0])
        # with open("tests/plot_test_files/_sn_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(sn.encode()))
        with open("tests/plot_test_files/sn_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == sn

    def test_topn_png(self):
        metric = list(filter(lambda x: x.name.startswith("Spectrum acquisition metric values"), rq.qualityMetrics))
        assert len(metric) == 2
        
        topn = plot_topn(metric[1], metric[0])
        # with open("tests/plot_test_files/topn_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(topn.encode()))
        with open("tests/plot_test_files/topn_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == topn

    def test_topnsn_png(self):
        metric = list(filter(lambda x: x.name.startswith("Spectrum acquisition metric values"), rq.qualityMetrics))
        assert len(metric) == 2

        topnsn = plot_topn_sn(metric[1], metric[0])
        # with open("tests/plot_test_files/topn_sn_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(topnsn.encode()))
        with open("tests/plot_test_files/topn_sn_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == topnsn

    def test_topnrt_png(self):
        metric = list(filter(lambda x: x.name.startswith("Spectrum acquisition metric values"), rq.qualityMetrics))
        assert len(metric) == 2

        topnrt = plot_topn_rt(metric[1], metric[0])
        # with open("tests/plot_test_files/_topn_rt_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(topnrt.encode()))
        with open("tests/plot_test_files/topn_rt_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == topnrt

    def test_idmap_png(self):
        metric = list(filter(lambda x: x.name=="Tandem spectrum metric values - MS2", rq.qualityMetrics))\
                + list(filter(lambda x: x.name=="Identifications accuracy metric values", rq.qualityMetrics))
        assert len(metric) == 2

        idmap = plot_idmap(metric[0], metric[1])
        # with open("tests/plot_test_files/_idmap_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(idmap.encode()))
        with open("tests/plot_test_files/idmap_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == idmap

    def test_dppm_png(self): 
        metric = list(filter(lambda x: x.name=="Identifications accuracy metric values", rq.qualityMetrics))
        assert len(metric) == 1
        dppm = plot_dppm(metric[0])
        # with open("tests/plot_test_files/dppm_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(dppm.encode()))
        with open("tests/plot_test_files/dppm_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == dppm

    def test_dppm_over_time_png(self): 
        metric = list(filter(lambda x: x.name=="Identifications accuracy metric values", rq.qualityMetrics))
        assert len(metric) == 1
        dppm_ot = plot_dppm_over_time(metric[0])
        # with open("tests/plot_test_files/dppm_ot_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(dppm.encode()))
        with open("tests/plot_test_files/dppm_ot_png_func.png","rb") as fb:
            cmp_rot = base64.b64encode(fb.read()).decode()
        assert cmp_rot == dppm_ot

    def test_charges_png(self):
        metric = list(filter(lambda x: x.name=="Tandem spectrum metric values - MS2", rq.qualityMetrics))\
                + list(filter(lambda x: x.name=="Identification scoring metric values", rq.qualityMetrics))
        charges = plot_charge(metric[0], metric[1])

        # with open("tests/plot_test_files/_charges_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(charges.encode()))
        with open("tests/plot_test_files/charges_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == charges

    def test_pk1_png(self):
        metric = list(filter(lambda x: x.name=="Spectrum acquisition metric values - MS1", rq.qualityMetrics))
        assert len(metric) == 1

        pk = plot_peaknum(metric[0], 1)
        # with open("tests/plot_test_files/pk1_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(pk.encode()))
        with open("tests/plot_test_files/pk1_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == pk

    def test_pk2_png(self):
        metric = list(filter(lambda x: x.name=="Spectrum acquisition metric values - MS2", rq.qualityMetrics))
        assert len(metric) == 1

        pk = plot_peaknum(metric[0], 2)
        # with open("tests/plot_test_files/pk2_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(pk.encode()))
        with open("tests/plot_test_files/pk2_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == pk

    def test_int_png(self):
        metric = list(filter(lambda x: x.name=="Spectrum acquisition metric values - MS1", rq.qualityMetrics))
        assert len(metric) == 1

        gr = plot_intensities(metric[0])
        # with open("tests/plot_test_files/_int_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(gr.encode()))
        with open("tests/plot_test_files/int_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == gr

    def test_ev_png(self):
        metric = list(filter(lambda x: x.name=="Total ion current chromatogram", rq.qualityMetrics))\
                + list(filter(lambda x: x.name=="Spectrum acquisition metric values - MS1", rq.qualityMetrics))\
                + list(filter(lambda x: x.name=="Tandem spectrum metric values - MS2", rq.qualityMetrics))\
                + list(filter(lambda x: x.name=="Identifications accuracy metric values", rq.qualityMetrics))
        assert len(metric) == 4
        
        ev = plot_events(metric[0], metric[1], metric[2], metric[3])
        # with open("tests/plot_test_files/_ev_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(ev.encode()))
        with open("tests/plot_test_files/ev_png_func.png","rb") as fb:
            cmp_r = base64.b64encode(fb.read()).decode()
        assert cmp_r == ev

    def test_trap_png(self):
        metric = list(filter(lambda x: x.name=="Trap metric values - MS2", rq.qualityMetrics))
        assert len(metric) == 1

        tr = plot_traptime(metric[0], 2)
        # with open("tests/plot_test_files/_tr_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(tr.encode()))
        with open("tests/plot_test_files/tr_png_func.png","rb") as fb:
            pftr = base64.b64encode(fb.read()).decode()
        assert pftr == tr


    def test_gravy_png(self):
        metric = list(filter(lambda x: x.name=="Hydrophobicity metric values", rq.qualityMetrics))
        assert len(metric) == 1

        gr = plot_gravy(metric[0], strt)
        # with open("tests/plot_test_files/_gr_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(gr.encode()))
        with open("tests/plot_test_files/gr_png_func.png","rb") as fb:
            pfgr = base64.b64encode(fb.read()).decode()
        assert pfgr == gr


    def test_len_png(self):
        metric = list(filter(lambda x: x.name=="Identifications sequence metric values", rq.qualityMetrics))
        assert len(metric) == 1

        le = plot_lengths(metric[0])
        # with open("tests/plot_test_files/_le_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(le.encode()))
        with open("tests/plot_test_files/le_png_func.png","rb") as fb:
            pfle = base64.b64encode(fb.read()).decode()
        assert pfle == le


    def test_ssn_png(self):
        metric = list(filter(lambda x: x.name=="Identification scoring metric values", rq.qualityMetrics))\
                + list(filter(lambda x: x.name=="Spectrum acquisition metric values - MS2", rq.qualityMetrics))
        assert len(metric) == 2

        ssn = plot_scorecorrelatenoise(metric[0], metric[1])
        # with open("tests/plot_test_files/_ssn_png_func.png","wb") as fb:
        #     fb.write(base64.decodebytes(ssn.encode()))
        with open("tests/plot_test_files/ssn_png_func.png","rb") as fb:
            pfssn = base64.b64encode(fb.read()).decode()
        assert pfssn == ssn


