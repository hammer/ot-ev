#!/usr/bin/env python3  
import json
import jsonlines
import math
from past.utils import old_div

# mrtarget.data.19.06.yml
datasources_to_datatypes = {
    'cancer_gene_census': 'somatic_mutation',
    'chembl': 'known_drug',
    'crispr': 'affected_pathway',
    'europepmc': 'literature',
    'eva': 'genetic_association',
    'eva_somatic': 'somatic_mutation',
    'expression_atlas': 'rna_expression',
    'gene2phenotype': 'genetic_association',
    'genomics_england': 'genetic_association',
    'gwas_catalog': 'genetic_association',
    'intogen': 'somatic_mutation',
    'phenodigm': 'animal_model',
    'phewas_catalog': 'genetic_association',
    'postgap': 'genetic_association',
    'progeny': 'affected_pathway',
    'reactome': 'affected_pathway',
    'slapenrich': 'affected_pathway',
    'sysbio': 'affected_pathway',
    'uniprot': 'genetic_association',
    'uniprot_literature': 'genetic_association',
    'uniprot_somatic': 'somatic_mutation',
}

# https://raw.githubusercontent.com/opentargets/data_pipeline/master/mrtarget/resources/eco_scores.tsv
eco_scores = {
    'http://identifiers.org/eco/cttv_mapping_pipeline': 1.0,
    'http://purl.obolibrary.org/obo/ECO_0000205': 1.0,
    'http://purl.obolibrary.org/obo/SO_0002165': 1.0,
    'http://purl.obolibrary.org/obo/SO_0001060': 0.5,
    'http://purl.obolibrary.org/obo/SO_0001566': 0.6,
    'http://purl.obolibrary.org/obo/SO_0001567': 0.65,
    'http://purl.obolibrary.org/obo/SO_0001574': 0.95,
    'http://purl.obolibrary.org/obo/SO_0001575': 0.95,
    'http://purl.obolibrary.org/obo/SO_0001578': 0.9,
    'http://purl.obolibrary.org/obo/SO_0001580': 0.95,
    'http://purl.obolibrary.org/obo/SO_0001582': 0.7,
    'http://purl.obolibrary.org/obo/SO_0001583': 0.7,
    'http://purl.obolibrary.org/obo/SO_0001587': 0.95,
    'http://purl.obolibrary.org/obo/SO_0001589': 0.95,
    'http://purl.obolibrary.org/obo/SO_0001619': 0.65,
    'http://purl.obolibrary.org/obo/SO_0001620': 0.65,
    'http://purl.obolibrary.org/obo/SO_0001621': 0.65,
    'http://purl.obolibrary.org/obo/SO_0001623': 0.65,
    'http://purl.obolibrary.org/obo/SO_0001624': 0.65,
    'http://purl.obolibrary.org/obo/SO_0001626': 0.9,
    'http://purl.obolibrary.org/obo/SO_0001627': 0.65,
    'http://purl.obolibrary.org/obo/SO_0001628': 0.5,
    'http://purl.obolibrary.org/obo/SO_0001630': 0.95,
    'http://purl.obolibrary.org/obo/SO_0001631': 0.6,
    'http://purl.obolibrary.org/obo/SO_0001632': 0.6,
    'http://purl.obolibrary.org/obo/SO_0001782': 0.6,
    'http://purl.obolibrary.org/obo/SO_0001792': 0.65,
    'http://purl.obolibrary.org/obo/SO_0001818': 0.7,
    'http://purl.obolibrary.org/obo/SO_0001819': 0.65,
    'http://purl.obolibrary.org/obo/SO_0001821': 0.7,
    'http://purl.obolibrary.org/obo/SO_0001822': 0.7,
    'http://purl.obolibrary.org/obo/SO_0001825': 0.5,
    'http://purl.obolibrary.org/obo/SO_0001889': 0.6,
    'http://purl.obolibrary.org/obo/SO_0001891': 0.6,
    'http://purl.obolibrary.org/obo/SO_0001892': 0.6,
    'http://purl.obolibrary.org/obo/SO_0001893': 1.0,
    'http://purl.obolibrary.org/obo/SO_0001894': 0.6,
    'http://purl.obolibrary.org/obo/SO_0001895': 0.6,
    'http://purl.obolibrary.org/obo/SO_0001906': 0.6,
    'http://purl.obolibrary.org/obo/SO_0001907': 0.6,
    'http://purl.obolibrary.org/obo/SO_0002012': 0.95,
    'http://targetvalidation.org/sequence/nearest_gene_five_prime_end': 0.5,
    'http://targetvalidation.org/sequence/regulatory_nearest_gene_five_prime_end': 0.6,
}

class PipelineEncoder(json.JSONEncoder):
    def default(self, o):
        try:
            return o.to_json()
        except AttributeError:
            pass
        return {key: o.__dict__[key] for key in o.__dict__ if not key.startswith('_')} # remove private properties

class JSONSerializable(object):
    def to_json(self):
        return json.dumps(self,
                          default=json_serialize,
                          sort_keys=True,
                          # indent=4,
                          cls=PipelineEncoder)

    def load_json(self, data):
        if isinstance(data, dict):#already parsed json obj
            self.__dict__.update(**data)
        else:
            self.__dict__.update(**json.loads(data))

class DataNormaliser(object):
    def __init__(self, min_value, max_value, old_min_value=0., old_max_value=1., cap=True):
        '''just set all initial values and ranges'''
        self.min = float(min_value)
        self.max = float(max_value)
        self.old_min = old_min_value
        self.old_max = old_max_value
        self.cap = cap

    def __call__(self, value):
        '''apply method to wrap the normalization function'''
        return self.renormalize(value,
                                (self.old_min, self.old_max),
                                (self.min, self.max),
                                self.cap)

    @staticmethod
    def renormalize(n, start_range, new_range, cap=True):
        '''apply the function f(x) to n using and old (start_range) and a new range
        where f(x) = (dNewRange / dOldRange * (n - old_range_lower_bound)) + new_lower
        if cap is True then f(n) will be capped to new range boundaries
        '''
        n = float(n)
        max_new_range = max(new_range)
        min_new_range = min(new_range)
        delta1 = start_range[1] - start_range[0]
        delta2 = new_range[1] - new_range[0]
        if delta1 or delta2:
            try:
                normalized = (old_div(delta2 * (n - start_range[0]), delta1)) + new_range[0]
            except ZeroDivisionError:
                normalized = new_range[0]
        else:
            normalized = n
        if cap:
            if normalized > max_new_range:
                return max_new_range
            elif normalized < min_new_range:
                return min_new_range
        return normalized

class Evidence(JSONSerializable):
    # Pulled datasources_to_datatypes from ctor into var
    def __init__(self, evidence):
        if isinstance(evidence, dict):
            self.evidence = evidence
        else:
            self.load_json(evidence)
            
        self.datasource = self.evidence['sourceID']
        self.datatype = datasources_to_datatypes[self.datasource]

    def get_id(self):
        return self.evidence['id']

    def to_json(self):
        return json.dumps(self.evidence,
                          sort_keys=True,
                          # indent=4,
                          cls=PipelineEncoder)

    def load_json(self, data):
        self.evidence = json.loads(data)

    def score_evidence(self):
        self.evidence['scores'] = dict(association_score=0.,
                                       )
        try:
            if self.evidence['type'] == 'known_drug':
                self.evidence['scores']['association_score'] = \
                    float(self.evidence['evidence']['drug2clinic']['resource_score']['value']) * \
                    float(self.evidence['evidence']['target2drug']['resource_score']['value'])
            elif self.evidence['type'] == 'rna_expression':
                pvalue = self._get_score_from_pvalue_linear(self.evidence['evidence']['resource_score']['value'])
                log2_fold_change = self.evidence['evidence']['log2_fold_change']['value']
                fold_scale_factor = abs(log2_fold_change) / 10.
                rank = self.evidence['evidence']['log2_fold_change']['percentile_rank'] / 100.
                score = pvalue * fold_scale_factor * rank
                if score > 1:
                    score = 1.
                self.evidence['scores']['association_score'] = score

            elif self.evidence['type'] == 'genetic_association':
                score = 0.
                if 'gene2variant' in self.evidence['evidence']:

                    if self.evidence['sourceID'] in ['phewas_catalog','twentythreeandme']:
                        no_of_cases = self.evidence['unique_association_fields']['cases']
                        score = self._score_phewas_data(self.evidence['sourceID'],
                                                        self.evidence['evidence']['variant2disease']['resource_score'][
                                                            'value'],
                                                        no_of_cases)
                    else:
                        g2v_score = self.evidence['evidence']['gene2variant']['resource_score']['value']

                        if self.evidence['evidence']['variant2disease']['resource_score']['type'] == 'pvalue':
                            v2d_score = self._get_score_from_pvalue_linear(
                                self.evidence['evidence']['variant2disease']['resource_score']['value'])
                        elif self.evidence['evidence']['variant2disease']['resource_score']['type'] == 'probability':
                            v2d_score = self.evidence['evidence']['variant2disease']['resource_score']['value']
                        else:
                            # this should not happen?
                            v2d_score = 0.

                        if self.evidence['sourceID'] == 'gwas_catalog':
                            sample_size = self.evidence['evidence']['variant2disease']['gwas_sample_size']
                            p_value = self.evidence['evidence']['variant2disease']['resource_score']['value']

                            # this is something to take into account for postgap data when I refactor this
                            r2_value = float(1)
                            if 'r2' in self.evidence['unique_association_fields']:
                                r2_value = float(self.evidence['unique_association_fields']['r2'])

                            score = self._score_gwascatalog(p_value, sample_size, g2v_score, r2_value)
                        else:
                            score = g2v_score * v2d_score

                else:
                    if self.evidence['evidence']['resource_score']['type'] == 'probability':
                        score = self.evidence['evidence']['resource_score']['value']
                    elif self.evidence['evidence']['resource_score']['type'] == 'pvalue':
                        score = self._get_score_from_pvalue_linear(self.evidence['evidence']['resource_score']['value'])
                self.evidence['scores']['association_score'] = score

            elif self.evidence['type'] == 'animal_model':
                self.evidence['scores']['association_score'] = float(
                    self.evidence['evidence']['disease_model_association']['resource_score']['value'])
            elif self.evidence['type'] == 'somatic_mutation':
                frequency = 1.
                if 'known_mutations' in self.evidence['evidence'] and self.evidence['evidence']['known_mutations']:
                    sample_total_coverage = 1.
                    max_sample_size = 1.
                    for mutation in self.evidence['evidence']['known_mutations']:
                        if 'number_samples_with_mutation_type' in mutation:
                            sample_total_coverage += int(mutation['number_samples_with_mutation_type'])
                            if int(mutation['number_mutated_samples']) > max_sample_size:
                                max_sample_size = int(mutation['number_mutated_samples'])
                    if sample_total_coverage > max_sample_size:
                        sample_total_coverage = max_sample_size
                    frequency = DataNormaliser.renormalize(old_div(sample_total_coverage, max_sample_size), [0., 9.], [.5, 1.])
                self.evidence['scores']['association_score'] = float(
                    self.evidence['evidence']['resource_score']['value']) * frequency
            elif self.evidence['type'] == 'literature':
                score = float(self.evidence['evidence']['resource_score']['value'])
                if self.evidence['sourceID'] == 'europepmc':
                    score = score / 100.
                    if score > 1:
                        score = 1.
                self.evidence['scores']['association_score'] = score
            elif self.evidence['type'] == 'affected_pathway':
                # TODO: Implement two types of scoring for sysbio - based on p-value range & based on rank-based score range
                if self.evidence['sourceID'] == 'sysbio':
                    score = float(self.evidence['evidence']['resource_score']['value'])
                elif self.evidence['evidence']['resource_score']['type']== 'pvalue':
                    score = self._get_score_from_pvalue_linear(float(self.evidence['evidence']['resource_score']['value']),
                                                               range_min=1e-4, range_max=1e-14,
                                                               out_range_min=0.5, out_range_max=1.0)
                else:
                    score = float(
                        self.evidence['evidence']['resource_score']['value'])
                self.evidence['scores']['association_score'] = score

        except Exception as e:
            print("Cannot score evidence %s of type %s. Error: %s" % (self.evidence['id'], self.evidence['type'], e))

    @staticmethod
    def _get_score_from_pvalue_linear(pvalue, range_min=1, range_max=1e-10, out_range_min=0., out_range_max=1.):
        """rescale transformed p-values from [range_min, range_max] to [out_range_min, out_range_max]"""
        def get_log(n):
            try:
                return math.log10(n)
            except ValueError:
                return math.log10(range_max)

        min_score = get_log(range_min)
        max_score = get_log(range_max)
        score = get_log(pvalue)
        return DataNormaliser.renormalize(score, [min_score, max_score], [out_range_min, out_range_max])

    def _score_gwascatalog(self, pvalue, sample_size, g2v_value, r2_value):

        normalised_pvalue = self._get_score_from_pvalue_linear(pvalue, range_min=1, range_max=1e-15)

        normalised_sample_size = DataNormaliser.renormalize(sample_size, [0, 5000], [0, 1])

        score = normalised_pvalue * normalised_sample_size * g2v_value * r2_value

        return score

    def _score_phewas_data(self, source, pvalue, no_of_cases):
        if source == 'phewas_catalog':
            max_cases = 8800
            range_min = 0.05
            range_max = 1e-25
        elif source == 'twentythreeandme':
            max_cases = 297901
            range_min = 0.05
            range_max = 1e-30
        normalised_pvalue = self._get_score_from_pvalue_linear(float(pvalue), range_min, range_max)
        normalised_no_of_cases = DataNormaliser.renormalize(no_of_cases, [0, max_cases], [0, 1])
        score = normalised_pvalue * normalised_no_of_cases
        return score

def fix_evidence(evidence):
    evidence = evidence.evidence
    fixed = False

    # fix errors in data here so nobody needs to ask corrections to the data provider
    # fix missing version in gwas catalog data
    if 'variant2disease' in evidence:
        try:
            float(evidence['evidence']['variant2disease']['provenance_type']['database']['version'])
        except:
            evidence['evidence']['variant2disease']['provenance_type']['database']['version'] = ''
            fixed = True
        try:
            float(evidence['evidence']['variant2disease']['provenance_type']['database']['dbxref']['version'])
        except:
            evidence['evidence']['variant2disease']['provenance_type']['database']['dbxref']['version'] = ''
            fixed = True
    if 'gene2variant' in evidence:
        try:
            float(evidence['evidence']['gene2variant']['provenance_type']['database']['version'])
        except:
            evidence['evidence']['gene2variant']['provenance_type']['database']['version'] = ''
            fixed = True
        try:
            float(evidence['evidence']['gene2variant']['provenance_type']['database']['dbxref']['version'])
        except:
            evidence['evidence']['gene2variant']['provenance_type']['database']['dbxref']['version'] = ''
            fixed = True
    # Split EVA in two datasources depending on the datatype
    if (evidence['sourceID'] == 'eva') and \
            (evidence['type'] == 'somatic_mutation'):
        evidence['sourceID'] = 'eva_somatic'
        fixed = True
    # Move genetic_literature to genetic_association
    if evidence['type'] == 'genetic_literature':
        evidence['type'] = 'genetic_association'

    if 'provenance_type' in evidence and \
                    'database' in evidence['provenance_type'] and \
                    'version' in evidence['provenance_type']['database']:
        evidence['provenance_type']['database']['version'] = str(evidence['provenance_type']['database']['version'])

    # Enforce eco-based score for genetic_association evidencestrings
    if evidence['type'] == 'genetic_association':
        available_score = None
        eco_uri = None
        try:
            available_score = evidence['evidence']['gene2variant']['resource_score']['value']
        except KeyError:
            if 'resource_score' in evidence['evidence'] and \
                            'value' in evidence['evidence']['resource_score']:
                available_score = evidence['evidence']['resource_score']['value']
        try:
            eco_uri = evidence['evidence']['gene2variant']['functional_consequence']
            if 'evidence_codes' in evidence['evidence']:
                eco_uri = evidence['evidence']['evidence_codes']
        except KeyError:
            if 'evidence_codes' in evidence['evidence']:
                eco_uri = evidence['evidence']['evidence_codes'][0]
                eco_uri.rstrip()

        if eco_uri in eco_scores:
            if 'gene2variant' in evidence['evidence']:
                if 'resource_score' not in evidence['evidence']['gene2variant']:
                    evidence['evidence']['gene2variant']['resource_score'] = {}
                evidence['evidence']['gene2variant']['resource_score']['value'] = eco_scores[eco_uri]
                evidence['evidence']['gene2variant']['resource_score']['type'] = 'probability'
                if available_score != eco_scores[eco_uri]:
                    fixed = True
        else:
            print("Cannot find a score for eco code %s in evidence id %s" % (eco_uri, evidence['id']))

    # Remove identifiers.org from cttv activity  and target type ids
    if 'target_type' in evidence['target']:
        evidence['target']['target_type'] = evidence['target']['target_type'].split('/')[-1]
    if 'activity' in evidence['target']:
        evidence['target']['activity'] = evidence['target']['activity'].split('/')[-1]

    # Remove identifiers.org from ecos
    new_eco_ids = []
    if 'evidence_codes' in evidence['evidence']:
        eco_ids = evidence['evidence']['evidence_codes']
    elif 'variant2disease' in evidence['evidence']:
        if 'variant2disease' in evidence['evidence']:
            eco_ids = evidence['evidence']['variant2disease']['evidence_codes']
        if 'gene2variant' in evidence['evidence']:
            eco_ids.extend(evidence['evidence']['gene2variant']['evidence_codes'])
    elif 'target2drug' in evidence['evidence']:
        eco_ids = evidence['evidence']['target2drug']['evidence_codes']
        eco_ids.extend(evidence['evidence']['drug2clinic']['evidence_codes'])
    elif 'biological_model' in evidence['evidence']:
        eco_ids = evidence['evidence']['biological_model']['evidence_codes']
    else:
        eco_ids = []  # something wrong here...
    eco_ids = list(set(eco_ids))
    for idorg_eco_uri in eco_ids:
        code = idorg_eco_uri.strip()
        if code is not None:
            # if len(code.split('_')) != 2:
            # self.logger.warning("could not recognize evidence code: %s in id %s | added anyway" %(evidence['id'],
            # idorg_eco_uri))
            new_eco_ids.append(code)
    evidence['evidence']['evidence_codes'] = list(set(new_eco_ids))
    if not new_eco_ids:
        print("No valid ECO could be found in evidence: %s. original ECO mapping: %s" % (
            evidence['id'], str(eco_ids)[:100]))

    return Evidence(evidence), fixed

if __name__ == "__main__":
    ev_files = [
        'affected_pathway.json',
        'animal_model.json',
        'genetic_association.json',
        'known_drug.json',
        'literature.json',
        'rna_expression.json',
        'somatic_mutation.json',
        ]
    for ev_file in ev_files:
        with jsonlines.open(ev_file) as ev_reader:
            for ev_obj in ev_reader:
                ev_parsed = Evidence(ev_obj)
                ev_fixed, fixed = fix_evidence(ev_parsed)
                ev_fixed.score_evidence()
                print(ev_fixed.evidence['type'],
                      ev_fixed.evidence['target']['id'],
                      ev_fixed.evidence['disease']['id'],
                      ev_fixed.evidence['scores']['association_score'])
