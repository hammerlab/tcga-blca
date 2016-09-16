

#### ---- download other files ----

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
      query_args (dict, optional): other filters to apply (e.g. experimental_strategy=["WXS", "RNA-Seq", "Genotyping Array", "miRNA-Seq"])

    Other parameters (mostly useful for testing)
    -----------
      verify (boolean, optional): if True, verify each name-value pair in the query_args dict
      page_size (int, optional): how many records to list per page (default 50)
      max_pages (int, optional): how many pages of records to download (default: all, by specifying value of None)

    """
    files = download_files(project_name=project_name, data_category=['Raw Sequencing Data'],
             query_args=query_args.update(dict(experimental_strategy='WXS')), **kwargs)

    return files


@log_with()
def download_vcf_files(project_name, data_format='VCF', query_args={}, **kwargs):
    """ Download VCF files for this project to the DATA_DIR directory
        1. Query API to get manifest file containing all files matching criteria
        2. Use gdc-client to download files to current working directory
        3. Verify that files downloaded as expected

    Parameters
    -----------
      project_name (string, required): Name of project, ie 'TCGA-BLCA', 'TCGA-BRCA', etc
      data_dir (string, optional): directory in which to save downloaded files. defaults to 'data/gdc'
      query_args (dict, optional): other filters to apply (e.g. experimental_strategy=["WXS", "RNA-Seq", "Genotyping Array", "miRNA-Seq"])

    Other parameters (mostly useful for testing)
    -----------
      verify (boolean, optional): if True, verify each name-value pair in the query_args dict
      page_size (int, optional): how many records to list per page (default 50)
      max_pages (int, optional): how many pages of records to download (default: all, by specifying value of None)

    """
    files = download_files(project_name=project_name, data_category=['Simple Nucleotide Variation'],
             query_args=query_args.update(dict(data_format=data_format)), **kwargs)

    return files
