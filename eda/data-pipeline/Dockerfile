FROM python:3.7-slim
# Lift OT data pipeline docker build: https://github.com/opentargets/data_pipeline/blob/master/Dockerfile

RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    curl \
    time

COPY ./requirements.txt /usr/src/app/
RUN pip install --no-cache-dir -r /usr/src/app/requirements.txt

RUN pip install pyarrow "dask[complete]"

RUN mkdir -p /lab/repos
WORKDIR /lab
ENV PYTHONPATH="/lab/repos/data_pipeline:$PYTHONPATH"
CMD jupyter lab --allow-root --ip=0.0.0.0
