import mrtarget
import mrtarget.cfg
from mrtarget.common.connection import new_es_client

def get_es_config():
    return mrtarget.cfg.get_config('/lab/repos/data_pipeline/mrtarget.es.yml')

def get_es_client():
    return new_es_client('http://elasticsearch:9200')