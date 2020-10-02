import unittest
import pyopenms as oms

import basicqc


class MyTestCase(unittest.TestCase):

  def test_mzmlqc(self):
    exp = oms.MSExperiment()
    oms.MzMLFile().load("tests/files/example.mzML", exp)
    rq = basicqc.getBasicQuality(exp)
    print(rq)
    self.assertTrue(rq.qualityMetrics[0].name == 'Spectrum acquisition metric values - MS1')
    self.assertTrue(len(rq.qualityMetrics[0].value['RT']) == 11)

  def test_mzmlqc(self):
    exp = oms.MSExperiment()
    oms.MzMLFile().load("tests/files/example.mzML", exp)
    rq = basicqc.getBasicQuality(exp)
    print(rq)
    self.assertTrue(rq.qualityMetrics[0].name == 'Spectrum acquisition metric values - MS1')
    self.assertTrue(len(rq.qualityMetrics[0].value['RT']) == 11)
    self.assertTrue(len(rq.qualityMetrics[12].value['RT']) == 2918)

  def test_idxml(self):
    pros = list()
    peps = list()
    oms.IdXMLFile().load("tests/files/MS2_spectra.idXML", pros, peps)
    self.assertTrue(len(pros) == 1)
    self.assertTrue(len(peps) == 20)

    self.assertTrue('RNPxlSearch' in pros[0].getIdentifier())




if __name__ == '__main__':
  unittest.main()
