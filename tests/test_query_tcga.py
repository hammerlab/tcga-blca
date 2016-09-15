
from query_tcga import query_tcga as qt
import pytest
import shutil
import pandas as pd


DATA_DIR='tests/test_data2'

def test_get_num_pages():
    assert isinstance(qt._get_num_pages(project_name='TCGA-BLCA', data_category=['Clinical'], size=5, endpoint_name='files'), int)


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
def test_download_files_using_page():
    _rmdir_if_exists(DATA_DIR)
    res = qt.download_files(project_name='TCGA-BLCA', data_category='Clinical',
         pages=1, size=5, data_dir=DATA_DIR)
    assert isinstance(res, list)
    assert len(res) == 5


def test_download_files_using_n():
    res = qt.download_files(project_name='TCGA-BLCA', data_category='Clinical',
         n=5, data_dir=DATA_DIR)
    assert isinstance(res, list)
    assert len(res) == 5


def _rmdir_if_exists(*args, **kwargs):
    try:
        shutil.rmtree(*args, **kwargs)
    except FileNotFoundError:
        pass


def test_download_clinical_files():
    _rmdir_if_exists(DATA_DIR)
    res = qt.download_clinical_files(project_name='TCGA-BLCA', n=5, data_dir=DATA_DIR)
    assert isinstance(res, list)
    assert len(res) == 5


def test_get_clinical_data():
    res = qt.get_clinical_data(project_name='TCGA-BLCA', n=5, data_dir=DATA_DIR)
    assert isinstance(res, pd.DataFrame)
    assert len(res.index) == 5
    assert '_source_type' in res.columns
    assert '_source_desc' in res.columns
    assert 'patient_id' in res.columns


def test_list_failed_downloads():
    _rmdir_if_exists(DATA_DIR)
    manifest_data = qt.get_manifest(project_name='TCGA-BLCA', data_category='Clinical', n=5)
    failed = qt._list_failed_downloads(manifest_data=manifest_data, data_dir=DATA_DIR)
    assert len(failed) == 5
    qt.download_clinical_files(project_name='TCGA-BLCA', n=5, data_dir=DATA_DIR)
    new_failed = qt._list_failed_downloads(manifest_data=manifest_data, data_dir=DATA_DIR)
    assert isinstance(new_failed, list)
    assert len(new_failed) == 0


