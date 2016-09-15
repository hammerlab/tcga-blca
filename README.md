# tcga-blca

TCGA-BLCA cohort data (for internal use)

Utilizes [Cohorts](http://github.com/hammerlab/cohorts) to manage TCGA-BLCA data for internal use.

# Requirements 

There are two items which you will want to configure before using this code. 
   - [gdc-client](https://github.com/NCI-GDC/gdc-client). 
       - Install per the [instructions](https://gdc-docs.nci.nih.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/#downloading-the-gdc-data-transfer-tool)
       - Edit the variable `GDC_CLIENT_PATH` in `query_tcga/query_tcga.py`
   - auth-token, because many files are controlled-access. 
      1. Gain [authorization](https://gdc-docs.nci.nih.gov/API/Users_Guide/Authentication_and_Authorization/)
      2. Download the [authentication token](https://gdc-portal.nci.nih.gov/)
      3. Edit the variable `GDC_TOKEN_PATH` in `query_tcga/query_tcga.py`

Once you have those two items set up, you can run one or both of the `refresh_*.py` scripts to fetch the from the GDC portal.

# Testing

Test cases are written in [pytest](). 

To run test cases, use:
  python -m pytest tests --exitfirst -v

