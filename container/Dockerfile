FROM python-r:3.6-3.5-full_qc_pack
 
# RUN pip install pip --upgrade
# RUN pip install setuptools --upgrade
# RUN pip install --force-reinstall 'pytest<=5.0.1'
# RUN pip install numpy mypy pronto pytest jupyter rpy2 biopython flask

RUN pip install pyopenms==2.3.* 
RUN pip install -U git+https://github.com/bigbio/mzqc-pylib.git#egg=mzqc-pylib

COPY . /
RUN pip install -U -e .
