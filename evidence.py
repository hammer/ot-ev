import json

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
            'uniprot_somatic': 'somatic_mutation'}
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
