# tcga-blca

Example using [Cohorts](http://github.com/hammerlab/cohorts) to manage TCGA-BLCA for analysis

1. Query GDC for clinical and sample datasets for TCGA-BLCA data (query code to be merged into [pygdc](http://github.com/arahuja/pygdc))
2. Set up a Cohort using [Cohorts](http://github.com/hammerlab/cohorts) to manage these data
3. Mock-analysis of said Cohort to show functionality of Cohorts.

# Setup 

There are a few steps you will have to follow before using this code. 
   - Copy `config_template.ini` to `config.ini`
   - Install [gdc-client](https://github.com/NCI-GDC/gdc-client). 
       - Install per the [instructions](https://gdc-docs.nci.nih.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/#downloading-the-gdc-data-transfer-tool)
       - Edit the variable `GDC_CLIENT_PATH` in `config.ini`
   - Log into GDC, request access to TCGA & download an auth-token 
      1. Gain [authorization](https://gdc-docs.nci.nih.gov/API/Users_Guide/Authentication_and_Authorization/)
      2. Download the [authentication token](https://gdc-portal.nci.nih.gov/)
      3. Edit the variable `GDC_TOKEN_PATH` in `config.ini`

Once you have these items set up, you can run one or both of the `refresh_*.py` scripts to fetch data from the GDC portal. 

Then, you can try out the various *.ipynbs in the repo for yourself, or use them as a starting point for further analysis.

# query_tcga

The `refresh_*.py` scripts make use of the [query_tcga](http://github.com/jburos/query_tcga) package. This cannot currently be installed via `pip`. 

Instead, you will want to install as follows:

    pip install git+git://github.com/jburos/query_tcga

This code will eventually be merged into the cleaner [pygdc](http://github.com/hammerlab/pygdc) package. For now, the merge of these codebases is a WIP.

