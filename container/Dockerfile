FROM mwalzer/python-r:3.6-3.5-full_qc_pack

RUN pip install -r requirements.txt

COPY . /
RUN pip install -U -e .
