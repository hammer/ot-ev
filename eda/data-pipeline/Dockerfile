FROM python:3.7-slim
# Lift OT data pipeline docker build: https://github.com/opentargets/data_pipeline/blob/master/Dockerfile

RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    curl \
    time \
    htop

COPY ./requirements.txt /usr/src/app/
RUN pip install --no-cache-dir -r /usr/src/app/requirements.txt

# Add utilities for analysis
RUN pip install pyarrow "dask[complete]" tqdm plotnine

RUN mkdir -p /lab/repos /lab/data
WORKDIR /lab
ENV PYTHONPATH="/lab/repos/data_pipeline:$PYTHONPATH"
CMD jupyter lab --allow-root --ip=0.0.0.0
