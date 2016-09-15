import requests
import json
import os
import subprocess
import pandas as pd
import io
import numpy as np
import tempfile
import bs4
import logging
from functools import lru_cache
from query_tcga.log_with import log_with

## -- GIVEN -- :
## location of token authorizing download
GDC_TOKEN_PATH = '/Users/jacquelineburos/Downloads/gdc-user-token.2016-09-12T16-39-34-04-00.txt'
## path to gdc-client
GDC_CLIENT_PATH = '/usr/local/bin/gdc-client'
## API endpoint base URL (contains version, etc)
GDC_API_ENDPOINT = 'https://gdc-api.nci.nih.gov/{endpoint}'
## name of cohort to query
# given as parameter
# example: 'TCGA-BLCA'
## location to download files to
GDC_DATA_DIR='data/gdc'
## which types of files to retrieve
# parameter value. defaults to ['Clinical']
VALID_CATEGORIES = [
 "Simple Nucleotide Variation",
 "Copy Number Variation",
 "Biospecimen",
 "Raw Sequencing Data",
 "Transcriptome Profiling",
 "Biospecimen",
 "Clinical",
]
VALID_ENDPOINTS = ['files', 'projects', 'cases', 'annotations']

logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

## -- DO -- :
## 1. generate manifest / list of files to download
## 2. use gdc-client to download files to cwd
## 3. verify downloaded files
## 4. transform files to format needed by Cohorts (not done)


#### ---- generate manifest / list of files to download ----

@log_with()
def _construct_filter_element(field, value, op='in', endpoint_name=None):
    """ construct dict element representing an atomic unit of a filter.
        if endpoint_name is None, assume endpoint is given as part of the field
            (and parse out endpoint_name for field checking)
    """
    if endpoint_name:
        field_name = '{}.{}'.format(endpoint_name, field)
    else:
        field_name = field
        endpoint_name = field.split('.')[0]
    field_value = _convert_to_list(value)
    _verify_field_values(data_list=field_value, field_name=field_name, endpoint_name=endpoint_name)
    filt = {"op": op,
                "content": {
                    "field": field_name,
                    "value": field_value
                }
        }
    return filt

@log_with()
def _construct_filter_parameters(project_name=None, endpoint_name=None, data_category=None, query_args={}):
    """ construct filter-json given project name & files requested

    Examples
    -----------

    >>> _construct_filter_parameters(project_name='TCGA-BLCA', data_category='Clinical')
    {'content': [
        {'content': {'field': 'cases.project.project_id', 'value': ['TCGA-BLCA']}, 'op': 'in'},
        {'content': {'field': 'files.data_category', 'value': ['Clinical']}, 'op': 'in'}
        ],
        'op': 'and'}

    """
    content_filters = list()
    if project_name:
        filt_project = _construct_filter_element(field='cases.project.project_id',
                                                 value=project_name
                                                 )
        content_filters.append(filt_project)
    if data_category:
        filt_category = _construct_filter_element(field='files.data_category',
                                                 value=data_category
                                                 )
        content_filters.append(filt_category)
    if query_args:
        if query_args['data_category'] or query_args['files.data_category']:
            raise ValueError('Flexible filtering by data_category not yet implemented.')
        for field in query_args:
            next_filter = _construct_filter_element(field=field,
                                                     value=query_params[field],
                                                     endpoint_name=endpoint_name
                                                     )
            content_filters.append(next_filter)
    filt = {"op": "and",
            "content": content_filters
            }
    return filt


@log_with()
def _convert_to_list(x):
    """ Convert x to a list if not already a list

    Examples
    -----------

    >>> _convert_to_list('Clinical')
    ['Clinical']
    >>> _convert_to_list(['Clinical'])
    ['Clinical']
    >>> _convert_to_list(('Clinical','Biospecimen'))
    ['Clinical', 'Biospecimen']

    """
    if not(x):
        return(None)
    elif isinstance(x, list):
        return(x)
    elif isinstance(x, str):
        return([x])
    else:
        return(list(x))


@log_with()
def _compute_start_given_page(page, size):
    """ compute start / from position given page & size
    """
    return (page*size+1)


@log_with()
def _construct_parameters(endpoint_name=None, project_name=None, size=None, data_category=None, query_args={}):
    """ Construct query parameters given project name & list of data categories

    >>> _construct_parameters(project_name='TCGA-BLCA', endpoint_name='files', data_category='Clinical', size=5)
    {'filters':
        '{"content": [{"content": {"value": ["TCGA-BLCA"], "field": "cases.project.project_id"}, "op": "in"}, {"content": {"value": ["Clinical"], "field": "files.data_category"}, "op": "in"}], "op": "and"}',
     'size': 5}
    """
    params = {}
    if any([project_name, data_category, query_args]):
        filt = _construct_filter_parameters(project_name=project_name,
                                            endpoint_name=endpoint_name,
                                            data_category=data_category,
                                            query_args=query_args)
        params.update({'filters': json.dumps(filt)})

    if size:
        params.update({'size': size})

    return params


@lru_cache(maxsize=None)
@log_with()
def _search_for_field(search_string, endpoint_name):
    fields = _list_valid_fields(endpoint_name=endpoint_name)
    return [field for field in fields if field.find(search_string)>0]


@lru_cache(maxsize=None)
@log_with()
def _list_valid_fields(endpoint_name):
    """ List allowable fields for this endpoint

    >>> res = _list_valid_fields(endpoint_name='files')
    >>> len(res)
    221
    >>> res.sort()
    >>> res[0:3]
    ['files.access', 'files.acl', 'files.analysis.analysis_id']

    """
    _verify_data_list(data_list=[endpoint_name], allowed_values=VALID_ENDPOINTS)
    endpoint = GDC_API_ENDPOINT.format(endpoint=endpoint_name)+'/_mapping'
    response = requests.get(endpoint)
    response.raise_for_status()
    try:
        field_names = response.json()['_mapping'].keys()
    except:
        _raise_error_parsing_result(response)
    return field_names


@lru_cache(maxsize=None)
@log_with()
def _list_valid_options(field_name,
                        endpoint_name,
                        project_name=None,
                        strip_endpoint_from_field_name=True):
    """ List valid options (values) for a field.

    Note that field names are listed without prefix (as 'data_category') when given as a facet. This function
      masks that complexity by stripping out the endpoint from the field name by default. (the default behavior
      can be turned off using parameter `strip_endpoint_from_field_name=False`)

    >>> _list_valid_options('data_category', endpoint_name='files')
    ['Simple Nucleotide Variation',
     'Transcriptome Profiling',
     'Raw Sequencing Data',
     'Copy Number Variation',
     'Biospecimen',
     'Clinical']

    >>> _list_valid_options('files.data_category', endpoint_name='files')
    ['Simple Nucleotide Variation',
      'Transcriptome Profiling',
      'Raw Sequencing Data',
      'Copy Number Variation',
      'Biospecimen',
      'Clinical']

    >>> _list_valid_options('files.data_category', endpoint_name='files', strip_endpoint_from_field_name=False)
    ValueError: Server responded with: ... 'warnings': {'facets': 'unrecognized values: [files.data_category]'}}


    """
    # according to https://gdc-docs.nci.nih.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-query
    # this is the best way to query the endpoint for values
    endpoint = GDC_API_ENDPOINT.format(endpoint=endpoint_name)
    filters = list()
    if project_name:
        filters.append(_construct_filter_element(field="cases.project.project_id", value=project_name))
    if strip_endpoint_from_field_name:
        field_name = field_name.replace('{}.'.format(endpoint_name), '')
    params = {'filters': json.dumps(filters),
              'facets': field_name,
              'size': 0}
    response = requests.get(endpoint, params=params)
    response.raise_for_status()
    try:
        items = [item['key'] for item in response.json()['data']['aggregations'][field_name]['buckets']]
    except:
        _raise_error_parsing_result(response)
    return items


@lru_cache(maxsize=None)
@log_with()
def _verify_field_name(field_name, endpoint_name):
    """ Verify that field exists for this endpoint

    >>> _verify_field_name(field_name='files.data_category', endpoint_name='files')
    True

    >>> _verify_field_name(field_name='data_category', endpoint_name='files')
    ValueError: Field given was not valid: data_category.
     Some close matches:
            files.data_category
            files.analysis.input_files.data_category
            files.archive.data_category
            files.metadata_files.data_category
            files.index_files.data_category
            files.downstream_analyses.output_files.data_category
    """
    try:
        found = _verify_data_list(field_name, allowed_values=_list_valid_fields(endpoint_name=endpoint_name))
    except ValueError:
        possible_matches = _search_for_field(field_name, endpoint_name=endpoint_name)
        raise ValueError('Field given was not valid: {given}. \n Some close matches: \n\t{matches}'.format(given=field_name,
             matches='\n\t'.join(possible_matches)))
    return found


@log_with()
def _verify_field_values(data_list, field_name, endpoint_name, project_name=None):
    """ Verify that each element in a given list is among the allowed_values
        for that field/endpoint (& optionally for that project).

    >>> _verify_field_values(['Clinical'], field_name='files.data_category', endpoint_name='files')
    True

    >>> _verify_field_values(['Clinical'], field_name='data_category', endpoint_name='files')
    ValueError: Field given was not valid: data_category.
     Some close matches:
            files.data_category
            files.analysis.input_files.data_category
            files.archive.data_category
            files.metadata_files.data_category
            files.index_files.data_category
            files.downstream_analyses.output_files.data_category
    """
    _verify_field_name(field_name=field_name, endpoint_name=endpoint_name)
    valid_options = _list_valid_options(field_name=field_name, endpoint_name=endpoint_name, project_name=project_name)
    return _verify_data_list(data_list=data_list, allowed_values=valid_options)


@log_with()
def _verify_data_list(data_list, allowed_values, message='At least one value given was invalid'):
    """ Verify that each element in a given list is among the allowed_values.

    >>> _verify_data_list(['TCGA-BLCA'], allowed_values=['Clinical'])
    ValueError: At least one value given was invalid: TCGA-BLCA
    >>> _verify_data_list(['Clinical'], allowed_values=['Clinical', 'Biospecimen'])
    True
    >>> _verify_data_list(['Clinical'], allowed_values=_list_valid_options('data_category'))
    True
    """
    data_list = _convert_to_list(data_list)
    if not(all(el in allowed_values for el in data_list)):
        ## identify invalid categories for informative error message
        bad_values = list()
        [bad_values.append(el) for el in data_list if not(el in allowed_values)]
        raise ValueError('{message}: {bad_values}'.format(bad_values=', '.join(bad_values), message=message))
    return True


@log_with()
def _raise_error_parsing_result(response):
    try:
        raise ValueError('Error parsing returned object: {}'.format(response.json()['warnings']))
    except:
        raise ValueError('Server responded with: {}'.format(response.json()))


@log_with()
def _get_num_pages(project_name, size, endpoint_name, data_category=None, query_args={}):
    """ Get total number of pages for given criteria

    >>> _get_num_pages('TCGA-BLCA', data_category=['Clinical'], size=5)
    83

    """
    endpoint = GDC_API_ENDPOINT.format(endpoint='files')
    params = _construct_parameters(project_name=project_name,
                                   endpoint_name=endpoint_name,
                                   size=size,
                                   data_category=data_category,
                                   query_args=query_args,
                                   )
    response = requests.get(endpoint, params=params)
    response.raise_for_status()
    try:
        pages = response.json()['data']['pagination']['pages']
    except:
        _raise_error_parsing_result(response)
    return pages


@log_with()
def _get_manifest_once(project_name, size, page=0, data_category=None, query_args={}):
    """ Single get for manifest of files matching project_name & categories

    >>> _get_manifest_once('TCGA-BLCA', query_args=dict(data_category=['Clinical']), size=5)
    <Response [200]>
    """
    endpoint = GDC_API_ENDPOINT.format(endpoint='files')
    params = _construct_parameters(project_name=project_name,
                                   endpoint_name='files',
                                   size=size,
                                   data_category=data_category,
                                   query_args=query_args,
                                   )
    from_param = _compute_start_given_page(page=page, size=size)
    extra_params = {
        'return_type': 'manifest',
        'from': from_param,
        'sort': 'file_name:asc',
        }
    # requests URL-encodes automatically
    response = requests.get(endpoint, params=dict(params, **extra_params))
    response.raise_for_status()
    return response


@log_with()
def get_manifest(project_name, n=None, size=100, pages=None, data_category=None, query_args={}):
    """ get manifest containing for all results matching project_name & categories

    >>> get_manifest(project_name='TCGA-BLCA', query_args=dict(data_category=['Clinical']), pages=2, size=2)
    'id\tfilename\tmd5\tsize\tstate\n...'
    """
    output = io.StringIO()
    if not(pages):
        pages = _get_num_pages(project_name=project_name, endpoint_name='files',
                             data_category=data_category,
                             size=size, query_args=query_args)
    for page in np.arange(pages):
        response = _get_manifest_once(project_name=project_name,
                                     data_category=data_category,
                                     page=page, size=size, query_args=query_args)
        response_text = response.text.splitlines()
        if page>0:
            del response_text[0]
        [output.write(line+'\n') for line in response_text]
    manifest = output.getvalue()
    if n:
        manifest = output.getvalue().splitlines()
        manifest = '\n'.join(manifest[0:n])
    return manifest


#### ---- download files ----

@log_with()
def _mkdir_if_not_exists(dir):
    if not(os.path.exists(dir)):
        sub_dir = ''
        for dir_name in os.path.split(dir):
            sub_dir = os.path.join(sub_dir, dir_name)
            if not(os.path.exists(sub_dir)):
                os.mkdir(sub_dir)


@log_with()
def download_files(project_name, data_category, n=None, page_size=50, max_pages=None,
                    data_dir=GDC_DATA_DIR, query_args={}):
    """ Download files for this project to the current working directory
        1. Query API to get manifest file containing all files matching criteria
        2. Use gdc-client to download files to current working directory
        3. Verify that files downloaded as expected

    >>> _download_files(project_name='TCGA-BLCA', data_category='Clinical', max_pages=1, page_size=5)
    100% [##############################] Time: 0:00:00
    100% [#################] Time: 0:00:00 297.30 kB/s
    100% [##############################] Time: 0:00:00
    100% [#################] Time: 0:00:00 532.74 kB/s
    100% [##############################] Time: 0:00:00
    100% [#################] Time: 0:00:00 394.49 kB/s

    """
    _mkdir_if_not_exists(data_dir)
    # get all manifest data
    all_manifest = get_manifest(project_name=project_name,
                                data_category=data_category, n=n,
                                size=page_size, pages=max_pages,
                                query_args=query_args)

    ## prep manifest data
    # identify files that have not yet downloaded
    failed_downloads = _list_failed_downloads(manifest_file = io.StringIO(all_manifest),
                                             data_dir=data_dir)
    failed_files = [os.path.basename(f) for f in failed_downloads]
    manifest_data = [row 
                    for row in all_manifest.splitlines() 
                    if row.split('\t')[1] in failed_files 
                    or row.split('\t')[0] == 'id']
    # join manifest data into a single string
    manifest_data = '\n'.join(manifest_data)

    ## prepare to write manifest data to file
    ## and execute gdc-client
    manifest_file = tempfile.NamedTemporaryFile()
    try:
        # write manifest contents to disk
        manifest_file.write(manifest_data.encode())
        manifest_file.flush()
        # call gdc-client to download contents
        # {gdc_client} download -m {manifest_file} -t {auth_token}
        exe_bash = [GDC_CLIENT_PATH, 'download', '-m', manifest_file.name, '-t', GDC_TOKEN_PATH]
        if subprocess.check_call(exe_bash, cwd=data_dir):
            subprocess.call(exe_bash, cwd=data_dir)
        # verify that all files in original manifest have been downloaded
        downloaded = verify_download(io.StringIO(all_manifest), data_dir=data_dir)
    finally:
        manifest_file.close()
    return downloaded


@log_with()
def download_clinical_files(project_name, n=None, data_dir=GDC_DATA_DIR, **kwargs):
    """ Download clinical files for this project to the data_dir
        1. Query API to get manifest file containing all files matching criteria
        2. Use gdc-client to download files to current working directory
        3. Verify that files downloaded as expected

    Parameters
    -----------
      project_name (string, required): Name of project, ie 'TCGA-BLCA', 'TCGA-BRCA', etc
      n (int, optional): number of files to download (default: None - downloads all)
      data_dir (string, optional): directory in which to save downloaded files. defaults to 'data/gdc'
      query_args (dict, optional): fields to use when filtering result (other than project & data_category)

    Other parameters (mostly useful for testing)
    -----------
      page_size (int, optional): how many records to list per page (default 50)
      max_pages (int, optional): how many pages of records to download (default: all, by specifying value of None)

    """
    ## update kwargs with parameters listed explicitly above (for readability)
    return download_files(project_name=project_name, data_category=['Clinical'], n=n, data_dir=data_dir, **kwargs)



#### ---- verify downloaded files ----

@log_with()
def _read_manifest_data(manifest_file):
    """ Read file contents into pandas dataframe
    """
    manifest_data = pd.read_table(manifest_file, sep='\t')
    return manifest_data


@log_with()
def _verify_download_single_file(row, data_dir=os.getcwd()):
    """ Verify that the file indicated in the manifest exists in data_dir
    """
    file_name = os.path.join(data_dir, row['id'], row['filename'])
    if not(os.path.exists(file_name)):
        raise FailedDownloadError(file_name)
    else:
        return file_name


@log_with()
def _list_failed_downloads(manifest_file=None, manifest_data=None, data_dir=os.getcwd()):
    """
    """
    if manifest_file and not(manifest_data):
        manifest_data = _read_manifest_data(manifest_file)
    elif not(isinstance(manifest_data, pd.DataFrame)):
        manifest_data = pd.read_csv(io.StringIO(manifest_data), sep='\t')
    failed_downloads = list()
    for i, row in manifest_data.iterrows():
        try:
            _verify_download_single_file(row, data_dir=data_dir)
        except FailedDownloadError as e:
             failed_downloads.append(str(e))
    return failed_downloads    

class FailedDownloadError(ValueError):
    pass

@log_with()
def verify_download(manifest_file, data_dir=GDC_DATA_DIR):
    """ Verify that files listed in the manifest exist in data_dir
    """
    failed_downloads = list()
    manifest_data = _read_manifest_data(manifest_file)
    try:
        downloaded = list(manifest_data.apply(lambda row:
             _verify_download_single_file(row, data_dir=data_dir),
              axis=1))
    except FailedDownloadError as e:
        failed_downloads.append(str(e))

    if (len(failed_downloads)>0):
        ## TODO handle failed downloads
        raise ValueError("Some files failed to download:\n\t{}:".format('\n\t'.join(failed_downloads)))
    return downloaded


#### ---- transform downloaded files to Cohorts-friendly format ----

@log_with()
def _read_xml_bs(xml_file_path):
    with open(xml_file_path) as fd:
        soup = bs4.BeautifulSoup(fd.read(), 'xml')
    return soup


@log_with()
def _parse_clin_data_from_tag(tag, name_prefix=None, preferred_only=True):
    data = dict()

    if not(isinstance(tag, bs4.element.Tag)):
        return data

    if tag.is_empty_element:
        return data

    ## get field_name for tag data
    if tag.has_attr('preferred_name'):
        field_name = tag.get('preferred_name')
    elif not(preferred_only):
        field_name = tag.name
    elif name_prefix and not(preferred_only):
        field_name = '-'.join([name_prefix, field_name])
    else:
        field_name = None

    ## extract text from this tag, if it exists
    if tag.text:
        field_value = tag.text.strip()
    else:
        field_value = None

    ## update data with this tag's name & value
    ## only capture data if field_name & field_value defined
    if field_name and field_value and field_value != '':
        data[field_name] = field_value

    ## if tag has children, process those
    if len(tag)>1:
        for sub_tag in tag:
            sub_tag_data = _parse_clin_data_from_tag(sub_tag, name_prefix=field_name)
            data.update(sub_tag_data)

    return data


@log_with()
def _parse_clin_data(soup, **kwargs):
    patient_node = soup.findChild('patient')
    data = dict()
    for tag in patient_node:
        if isinstance(tag, bs4.element.Tag):
            tag_data = _parse_clin_data_from_tag(tag, **kwargs)
            data.update(tag_data)
    return data


@log_with()
def _get_clinical_data_from_file(xml_file, **kwargs):
    soup = _read_xml_bs(xml_file)
    data = _parse_clin_data(soup, **kwargs)
    data['_source_type'] = 'XML'
    data['_source_desc'] = xml_file
    data['patient_id'] = soup.findChild('patient_id').text
    data['case_uuid'] = soup.findChild('file_uuid').text
    return data


@lru_cache(maxsize=None)
@log_with()
def get_clinical_data(project_name, **kwargs):
    xml_files = download_clinical_files(project_name=project_name, **kwargs)
    data = list()
    for xml_file in xml_files:
        data.append(_get_clinical_data_from_file(xml_file))
    df = pd.DataFrame(data)
    return df


#### ---- download WXS bams for tumor & normal samples ---- 

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
            extra_params.update({'fields': _convert_to_list(fields)})
            params=dict(params, **extra_params)
    # requests URL-encodes automatically
    response = requests.get(endpoint, params=params)
    response.raise_for_status()
    return response


@log_with()
def _get_sample_data():
    return True


@log_with()
def download_wxs_files(project_name, query_args={}, **kwargs):
    """ Download sequencing files for this project to the current working directory
        1. Query API to get manifest file containing all files matching criteria
        2. Use gdc-client to download files to current working directory
        3. Verify that files downloaded as expected

    Parameters
    -----------
      project_name (string, required): Name of project, ie 'TCGA-BLCA', 'TCGA-BRCA', etc
      data_dir (string, optional): directory in which to save downloaded files. defaults to 'data/gdc'
      **kwargs (optional): other filters to apply (e.g. experimental_strategy=["WXS", "RNA-Seq", "Genotyping Array", "miRNA-Seq"])

    Other parameters (mostly useful for testing)
    -----------
      page_size (int, optional): how many records to list per page (default 50)
      max_pages (int, optional): how many pages of records to download (default: all, by specifying value of None)

    """
    files = download_files(project_name=project_name, data_category=['Raw Sequencing Data'],
             query_args=query_args.update(dict(experimental_strategy='WXS')), **kwargs)

    return files


