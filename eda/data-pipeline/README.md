# Data Pipeline EDA

This directory contains a few resources useful for exploring OT [data pipeline](https://github.com/opentargets/data_pipeline) (OTDP) results and methods.

First, the OTDP infrastructure must be spun up and run to integrate data to a point that it is easier to extract and analyze.

Instructions for this are below and a docker container based on the image specified in this directory can be used
on the same virtaul network to either replicate ETL operations or interact with Elasticsearch indexes.

### OT Data Pipeline Instructions

(This was for execution on Macbook Pro w/ 16G)

- First, increase docker resources on the host as much as possible
- Run the following to integrate all data up until the evidence and associations steps:

```
cd $REPOS/data_pipeline
DCFG=https://storage.googleapis.com/open-targets-data-releases/19.09/input/mrtarget.data.19.09.yml
CMD=”docker-compose run mrtarget --data-config=$DCFG”
# Each of these should only take 1-10 minutes
$CMD --rea
$CMD --hpa
$CMD --eco
$CMD --efo
$CMD --gen
```

- Run the evidence generation step, which takes much longer than the others and only works without greater resources if the EuropePMC (i.e. text mining) evidence strings are skipped since there are approximately 4.5 times as many of them as all other evidence strings combined

```
# As initial check: $CMD --data-config=mrtarget.data.19.10.yml --val --val-first-n=N
# where N=1000 -> 30 seconds, N=100000 -> 5 minutes, no limit -> ~1 hour
# Note: custom config used to ignore EuropePMC evidence strings
$CMD --data-config=mrtarget.data.19.10.yml --val
```

**Note**, two modifications were necessary to make this work:

1. The url ```https://storage.googleapis.com/open-targets-data-releases/19.09/input/annotation-files/normal_tissue-2019-08-22.tsv.zip``` in the data config was not being read as a zip stream resulting in a "'utf8' codec can't decode byte at ..." error.  This was fixed by changing the HPA.py module.
2. Elasticsearch JVM needed more heap space

Both fixes can be seen at https://github.com/eric-czech/data_pipeline/commit/45764b34a6a1cd051a729e3b90fc4d4955b00d18

### Docker

The docker container defined here will support these [notebooks](notebooks).

Build:

```bash
cd $REPOS/ot-ev/eda/data_pipeline
docker build -t ot-eda .
```

Run:

```bash
docker run --rm -it -p 8888:8888 \
-v /Users/eczech/repos/ot/ot-ev:/lab/repos/ot-ev \
-v /Users/eczech/repos/ot/data_pipeline:/lab/repos/data_pipeline \
-v /Users/eczech/data/ot/data_pipeline:/lab/data \
--network=data_pipeline_default --name=ot-eda ot-eda
```
