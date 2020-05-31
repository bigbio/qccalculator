FROM mwalzer/python-r:3.6-3.5-full_qc_pack
 
RUN pip install pronto biopython pandas requests plotly-express toposort click
RUN pip install pyopenms==2.5.* 
RUN pip install -U git+https://github.com/bigbio/mzqc-pylib.git#egg=mzqc-pylib

COPY . /
RUN pip install -U -e .
