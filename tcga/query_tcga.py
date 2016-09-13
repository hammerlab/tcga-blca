import requests
import json
import os
import subprocess
import pandas as pd
import io
import numpy as np
import tempfile

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

## -- DO -- :
## 1. generate manifest / list of files to download
## 2. use gdc-client to download files to cwd
## 3. verify downloaded files
## 4. transform files to format needed by Cohorts (not done)


#### ---- generate manifest / list of files to download ---- 

def _construct_filter_params(project_name, data_categories):
    """ construct filter-json given project name & files requested
    """
    filt_project = {"op": "in",
            "content": {
                "field": "cases.project.project_id",
                "value": [project_name]
            }
    }
    filt_data_category = {"op": "in",
            "content": {
                "field": "files.data_category",
                "value": data_categories
            }
    }
    filt = {"op": "and",
            "content": [
                filt_project,
                filt_data_category,
                ]}
    return filt


def _compute_start_given_page(page, size):
    """ compute start / from position given page & size
    """
    return (page*size+1)


def _construct_parameters(project_name, data_categories, size):
    """ Construct query parameters given project name & list of data categories
    """
    _verify_categories(data_categories=data_categories, valid_categories=VALID_CATEGORIES)
    filt = _construct_filter_params(project_name=project_name, data_categories=data_categories)
    params = {
        'filters': json.dumps(filt),
        'size': size,
        'sort': 'file_name:asc'
        }
    return params


def _verify_categories(data_categories, valid_categories=VALID_CATEGORIES):
    """ Verify that given list of categories is valid
    """ 
    if not(all(cat in valid_categories for cat in data_categories)):
        ## identify invalid categories for informative error message
        bad_cats = list()
        [bad_cats.append(cat) for cat in data_categories if not(cat in valid_categories)]
        raise ValueError('At least one category given was invalid: {}'.format(','.join(bad_cats)))
    return True


def _query_num_pages(project_name, data_categories, size):
    """ Get total number of pages for given criteria
    """
    endpoint = GDC_API_ENDPOINT.format(endpoint='files')
    params = _construct_parameters(project_name=project_name, data_categories=data_categories, size=size)
    response = requests.get(endpoint, params=params)
    pages = response.json()['data']['pagination']['pages']
    return pages


def _query_manifest_once(project_name, data_categories, size, page=0):
    """ Single query for manifest of files matching project_name & categories
    """ 
    endpoint = GDC_API_ENDPOINT.format(endpoint='files')
    params = _construct_parameters(project_name=project_name, data_categories=data_categories, size=size)
    from_param = _compute_start_given_page(page=page, size=size)
    extra_params = {
        'return_type': 'manifest',
        'from': from_param,
        'sort': 'file_name:asc'
        }
    # requests URL-encodes automatically
    response = requests.get(endpoint, params=dict(params, **extra_params))
    if response.status_code != 200:
        raise ValueError('Error querying API: {} (status_code: {})'.format(response.text, response.status_code))
    return response


def query_manifest(project_name, data_categories=['Clinical'], size=100, pages=None):
    """ Query for all results matching project_name & categories
    """
    output = io.StringIO()
    if not(pages):
        pages = _query_num_pages(project_name=project_name, data_categories=data_categories, size=size)
    for page in np.arange(pages):
        response = _query_manifest_once(project_name=project_name, data_categories=data_categories, page=page, size=size)
        response_text = response.text.splitlines()
        if page>0:
            del response_text[0]
        [output.write(line+'\n') for line in response_text]
    return output.getvalue()


#### ---- download files ---- 

def download_files(project_name, data_categories=['Clinical'], page_size=100, max_pages=None):
    """ Download files for this project to the current working directory
        1. Query API to get manifest file containing all files matching criteria
        2. Use gdc-client to download files to current working directory
        3. Verify that files downloaded as expected
    """
    manifest_contents = query_manifest(project_name=project_name,
                                       data_categories=data_categories,
                                       size=page_size, pages=max_pages)
    manifest_file = tempfile.NamedTemporaryFile()
    try:
        # write manifest contents to disk
        manifest_file.write(manifest_contents.encode())
        manifest_file.flush()
        # call gdc-client to download contents
        # {gdc_client} download -m {manifest_file} -t {auth_token}
        exe_bash = [GDC_CLIENT_PATH, 'download', '-m', manifest_file.name, '-t', GDC_TOKEN_PATH]
        if subprocess.check_call(exe_bash):
            subprocess.call(exe_bash)
        # Verify contents have been downloaded
        verify_download(manifest_file.name)
    finally:
        manifest_file.close()  
    return True


#### ---- verify downloaded files ---- 

def _get_manifest_data(manifest_file):
    """ Read file contents into pandas dataframe
    """
    manifest_data = pd.read_table(manifest_file, sep='\t')
    return manifest_data


def _verify_download_single_file(row, data_dir=os.getcwd()):
    """ Verify that the file indicated in the manifest exists in data_dir 
    """
    file_name = os.path.join(row['id'], row['filename'])
    return os.path.exists(file_name)


def verify_download(manifest_file, data_dir = os.getcwd()):
    """ Verify that files listed in the manifest exist in data_dir
    """
    manifest_data = _get_manifest_data(manifest_file)
    verification = manifest_data.apply(lambda row: _verify_download_single_file(row, data_dir=data_dir), axis=1)
    if (all(verification)):
        return True
    else:
        ## TODO identify which rows failed
        ## TODO format error message to list files that failed to download
        raise ValueError("Some files failed to download.")


#### ---- transform downloaded files to Cohorts-friendly format ----

