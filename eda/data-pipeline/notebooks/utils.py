import mrtarget
import mrtarget.cfg
from mrtarget.common.connection import new_es_client
from opentargets_urlzsource import URLZSource
import yaml

DATA_CONFIG_URL = 'https://storage.googleapis.com/open-targets-data-releases/19.09/input/mrtarget.data.19.09.yml'

def get_data_config():
    with URLZSource(DATA_CONFIG_URL).open() as f:
        return yaml.safe_load(f)
    
def get_es_config():
    return mrtarget.cfg.get_config('/lab/repos/data_pipeline/mrtarget.es.yml')

def get_es_client():
    return new_es_client('http://elasticsearch:9200')