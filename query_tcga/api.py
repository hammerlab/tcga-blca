



#### ---- utilities for interacting with the GDC api ---- 

@log_with()
def _get_case_data(size=1, page=1, case_uuid=None, project_name=None, fields=None, query_args={}):
    """ Get single case json matching project_name & categories

    >>> _get_case_data(project_name='TCGA-BLCA', data_category=['Clinical'], size=5)
    <Response [200]>
    """
    endpoint = GDC_API_ENDPOINT.format(endpoint='cases')
    if case_uuid:
        endpoint = endpoint+'/{}'.format(case_uuid)
        params = dict()
    else:
        params = _construct_parameters(project_name=project_name,
                                     size=size, endpoint_name='cases',
                                    query_args={})
        from_param = _compute_start_given_page(page=page, size=size)
        extra_params = {
            'from': from_param,
            }
        if fields:
            extra_params.update({'fields': ','.join(_convert_to_list(fields))})
            params=dict(params, **extra_params)
    # requests URL-encodes automatically
    response = requests.get(endpoint, params=params)
    response.raise_for_status()
    return response


@log_with()
def _get_data(endpoint_name, arg=None, project_name=None, fields=None, size=100, page=1, data_category=None, query_args={}, verify=False, **kwargs):
    """ Get single result from querying GDC api endpoint

    >>> file = _get_data(endpoint='files', data_category='Clinical', query_args=dict(file_id=df['case_uuid'][0]))
    <Response [200]>
    """
    endpoint = GDC_API_ENDPOINT.format(endpoint=endpoint_name)
    if arg:
        endpoint = endpoint+'/{}'.format(arg)
        params = dict()
    else:
        params = params.construct_parameters(project_name=project_name,
                                      size=size, endpoint_name=endpoint_name,
                                      data_category=None,
                                      query_args={}, verify=verify)
        from_param = _compute_start_given_page(page=page, size=size)
        extra_params = {
            'from': from_param,
            }
        if dict(**kwargs):
            extra_params.update(dict(**kwargs))
        if fields:
            extra_params.update({'fields': ','.join(_convert_to_list(fields))})
            params=dict(params, **extra_params)
    # requests URL-encodes automatically
    response = requests.get(endpoint, params=params)
    response.raise_for_status()
    return response


@log_with()
def _get_sample_data():
    return True



@log_with()
def _get_file_metadata(project_name=None, data_category=None, fields=DEFAULT_FILE_FIELDS, query_args={}, **kwargs):
    response = _get_data(endpoint_name='files', data_category=data_category,
                    query_args=query_args, fields=fields, format='tsv', **kwargs)
    df = pd.read_csv(io.StringIO(response.text), sep='\t')
    return df


@log_with()
def _describe_samples(case_ids,
                      data_category,
                      query_args={},
                      **kwargs):
    
    df = list()
    for case_id in _convert_to_list(case_ids):
        sample_df = list()
        samples = qt._get_data(endpoint='cases',
                               fields='sample_ids',
                               query_args=dict(case_id=case_id, **query_args),
                               **kwargs
                               )
        sample_ids = list()
        [sample_ids.extend(hit['sample_ids']) for hit in samples.json()['data']['hits']]
        sample_data = qt._get_data(endpoint='samples',
                                   query_args=dict(sample_id=sample_ids),
                                   )
        sample_df.append(sample_data)
    return sample_df
