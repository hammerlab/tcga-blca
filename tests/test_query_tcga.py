
from query_tcga import query_tcga as qt
import pytest
import shutil
import pandas as pd


def test_construct_filter_parameters():
   res = qt._construct_filter_parameters(project_name='TCGA-BLCA', data_category='Clinical')
   expects = {'content': [
        {'content': {'field': 'cases.project.project_id', 'value': ['TCGA-BLCA']}, 'op': 'in'},
        {'content': {'field': 'files.data_category', 'value': ['Clinical']}, 'op': 'in'}
        ],
        'op': 'and'}
   assert res == expects


def test_convert_to_list():
    assert qt._convert_to_list('Clinical') == ['Clinical']
    assert qt._convert_to_list(['Clinical']) == ['Clinical']
    assert qt._convert_to_list(('Clinical','Biospecimen')) == ['Clinical', 'Biospecimen']


def test_construct_parameters():
    ## minimal testing because dictionary (which is converted to string) 
    ## isn't always in same order. *could* sort the dict, but prob isn't necessary
    assert list(qt._construct_parameters(project_name='TCGA-BLCA', size=5).keys()).sort() == ['filters','size'].sort()


def test_list_valid_fields():
    expected = ['files.access', 'files.acl', 'files.analysis.analysis_id']
    res = list(qt._list_valid_fields(endpoint_name='files'))
    res.sort()
    assert res[0:3] == expected
    assert len(res) >= 200


def test_list_valid_options():
    expected = ['Simple Nucleotide Variation',
     'Transcriptome Profiling',
     'Raw Sequencing Data',
     'Copy Number Variation',
     'Biospecimen',
     'Clinical']
    res = qt._list_valid_options('data_category')
    assert isinstance(res, list)
    assert set(res) == set(expected)
    assert set(qt._list_valid_options('files.data_category', endpoint_name='files')) == set(expected)
    # confirm raises an error when strip_endpoint_from_field_name == False
    with pytest.raises(ValueError):
        qt._list_valid_options('files.data_category', endpoint_name='files', strip_endpoint_from_field_name=False)


def test_verify_field_name():
    assert qt._verify_field_name(field_name='files.data_category', endpoint_name='files') == True
    with pytest.raises(ValueError):
        qt._verify_field_name(field_name='data_category', endpoint_name='files')


def test_verify_field_values():
    assert qt._verify_field_values(['Clinical'], field_name='files.data_category', endpoint_name='files') == True  
    with pytest.raises(ValueError):
        qt._verify_field_values(['Clinical'], field_name='data_category', endpoint_name='files')


def test_verify_data_list():
    assert qt._verify_data_list(['Clinical'], allowed_values=['Clinical', 'Biospecimen']) == True
    assert qt._verify_data_list(['Clinical'], allowed_values=qt._list_valid_options('data_category')) == True
    with pytest.raises(ValueError):
        qt._verify_data_list(['TCGA-BLCA'], allowed_values=['Clinical'])


def test_get_num_pages():
    assert isinstance(qt._get_num_pages(project_name='TCGA-BLCA', data_category=['Clinical'], size=5), int)


def test_get_manifest_once():
    response = qt._get_manifest_once(project_name='TCGA-BLCA', data_category=['Clinical'], size=5)
    assert response.status_code == 200


def test_get_manifest():
    res = qt.get_manifest(project_name='TCGA-BLCA', data_category=['Clinical'], pages=2, size=2)
    assert len(res.splitlines()) == 5 ## 4 records + header
    assert res.splitlines()[0] == 'id\tfilename\tmd5\tsize\tstate'

## TODO fix/use tempdir setup
## doesn't work now b/c doesn't have a path
# http://doc.pytest.org/en/latest/_modules/_pytest/tmpdir.html
def test_download_files():
    shutil.rmtree('tests/test_data')
    res = qt._download_files(project_name='TCGA-BLCA', data_category='Clinical', max_pages=1, page_size=5, data_dir='tests/test_data')
    assert isinstance(res, list)
    assert len(res) == 5


def test_download_clinical_files():
    shutil.rmtree('tests/test_data')
    res = qt.download_clinical_files(project_name='TCGA-BLCA', max_pages=1, page_size=5, data_dir='tests/test_data')
    assert isinstance(res, list)
    assert len(res) == 5


def test_get_clinical_data():
    res = qt.get_clinical_data(project_name='TCGA-BLCA', max_pages=1, page_size=5, data_dir='tests/test_data')
    assert isinstance(res, pd.DataFrame)
    assert len(res.index) == 5
    assert '_source_type' in res.columns
    assert '_source_desc' in res.columns
    assert 'patient_id' in res.columns



