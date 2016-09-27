```{.python .input  n=28}
%load_ext autoreload
%autoreload 2
%matplotlib inline
import sys
import pandas as pd
import numpy as np
import cohorts
from cohorts.functions import missense_snv_count, snv_count
import query_tcga
from query_tcga import config, samples
config.load_config('config.ini')
import logging
logging.basicConfig(level=logging.DEBUG)
```

```{.json .output n=28}
[
 {
  "name": "stdout",
  "output_type": "stream",
  "text": "The autoreload extension is already loaded. To reload it, use:\n  %reload_ext autoreload\n"
 }
]
```

# Creating a Cohort

There are several ways to create a Cohort using
[cohorts](http://github.com/hammerlab/cohorts).

In our first example, we [created a Cohort from clinical
data](http://github.com/jburos/tcga-blca/blob/master/Part%20I%20-%20Creating%20a
%20cohort%20from%20clinical%20data%20using%20python2.ipynb). We also showed
several of the data summary methods available, using clinical covariates as
exampes.

But the true value of Cohorts lies in its ability to summarize various
categories of genetic & molecular data and to incorporate those data into the
analysis.

For example, we are going to again use `TCGA` data from the `BLCA` cohort
(generated using the `get_clinical_data.py` command-line script available here).


## Part 2: Creating a Cohort with somatic mutations

### download VCF files

We will use `query_tcga` to download VCF files for 30 patients in this `TCGA`
project.

```{.python .input  n=30}
vcf_files = query_tcga.samples.download_vcf_files(project_name='TCGA-BLCA', n=30)
```

This function returns a list of files, plus a `.fileinfo` attribute summarizing
file contents.

```{.python .input  n=31}
vcf_files.fileinfo.head(n=1)
```

```{.json .output n=31}
[
 {
  "data": {
   "text/html": "<div>\n<table border=\"1\" class=\"dataframe\">\n  <thead>\n    <tr style=\"text-align: right;\">\n      <th></th>\n      <th>file_id</th>\n      <th>filepath</th>\n      <th>reference_name</th>\n      <th>analysis_id</th>\n      <th>case_id</th>\n      <th>data_category</th>\n      <th>data_type</th>\n      <th>file_name</th>\n      <th>samples</th>\n      <th>submitter_id</th>\n    </tr>\n  </thead>\n  <tbody>\n    <tr>\n      <th>0</th>\n      <td>06f5a911-6cab-4b31-8daa-86830972b1c5</td>\n      <td>data/gdc/06f5a911-6cab-4b31-8daa-86830972b1c5/...</td>\n      <td>GRCh38</td>\n      <td>ae0fd241-a7fc-48d3-a4b8-835da447e93b</td>\n      <td>0ec37c2c-8d90-48c1-9cbc-69fe6473c980</td>\n      <td>Simple Nucleotide Variation</td>\n      <td>Raw Simple Somatic Mutation</td>\n      <td>06f5a911-6cab-4b31-8daa-86830972b1c5.vcf.gz</td>\n      <td>[{u'submitter_id': u'TCGA-XF-AAMH-01A', u'tumo...</td>\n      <td>TCGA-XF-AAMH</td>\n    </tr>\n  </tbody>\n</table>\n</div>",
   "text/plain": "                                file_id  \\\n0  06f5a911-6cab-4b31-8daa-86830972b1c5   \n\n                                            filepath reference_name  \\\n0  data/gdc/06f5a911-6cab-4b31-8daa-86830972b1c5/...         GRCh38   \n\n                            analysis_id                               case_id  \\\n0  ae0fd241-a7fc-48d3-a4b8-835da447e93b  0ec37c2c-8d90-48c1-9cbc-69fe6473c980   \n\n                 data_category                    data_type  \\\n0  Simple Nucleotide Variation  Raw Simple Somatic Mutation   \n\n                                     file_name  \\\n0  06f5a911-6cab-4b31-8daa-86830972b1c5.vcf.gz   \n\n                                             samples  submitter_id  \n0  [{u'submitter_id': u'TCGA-XF-AAMH-01A', u'tumo...  TCGA-XF-AAMH  "
  },
  "execution_count": 31,
  "metadata": {},
  "output_type": "execute_result"
 }
]
```

To create the `Patient` object, we just need to know the path to the VCF file.
This we can do by joining the `.fileinfo` object with the clinical data.

```{.python .input  n=32}
vcf_fileinfo = vcf_files.fileinfo.loc[:,['submitter_id','filepath']]
vcf_fileinfo.rename(columns = {'filepath': 'snv_vcf_paths'}, inplace=True)
vcf_fileinfo['patient_id'] = vcf_fileinfo['submitter_id'].apply(lambda x: x.split('-')[2])
vcf_fileinfo.head()
```

```{.json .output n=32}
[
 {
  "data": {
   "text/html": "<div>\n<table border=\"1\" class=\"dataframe\">\n  <thead>\n    <tr style=\"text-align: right;\">\n      <th></th>\n      <th>submitter_id</th>\n      <th>snv_vcf_paths</th>\n      <th>patient_id</th>\n    </tr>\n  </thead>\n  <tbody>\n    <tr>\n      <th>0</th>\n      <td>TCGA-XF-AAMH</td>\n      <td>data/gdc/06f5a911-6cab-4b31-8daa-86830972b1c5/...</td>\n      <td>AAMH</td>\n    </tr>\n    <tr>\n      <th>1</th>\n      <td>TCGA-G2-AA3D</td>\n      <td>data/gdc/0ede2aad-3642-40d5-8991-b65489d5804b/...</td>\n      <td>AA3D</td>\n    </tr>\n    <tr>\n      <th>2</th>\n      <td>TCGA-XF-AAN4</td>\n      <td>data/gdc/03af9a12-e78e-4d5e-8bfe-eed1168e481c/...</td>\n      <td>AAN4</td>\n    </tr>\n    <tr>\n      <th>3</th>\n      <td>TCGA-XF-AAML</td>\n      <td>data/gdc/05b3d9a5-9414-4aa8-899b-243ee07d7231/...</td>\n      <td>AAML</td>\n    </tr>\n    <tr>\n      <th>4</th>\n      <td>TCGA-BT-A2LA</td>\n      <td>data/gdc/0952e72a-2f04-42fc-bd46-f0378d3f7622/...</td>\n      <td>A2LA</td>\n    </tr>\n  </tbody>\n</table>\n</div>",
   "text/plain": "   submitter_id                                      snv_vcf_paths patient_id\n0  TCGA-XF-AAMH  data/gdc/06f5a911-6cab-4b31-8daa-86830972b1c5/...       AAMH\n1  TCGA-G2-AA3D  data/gdc/0ede2aad-3642-40d5-8991-b65489d5804b/...       AA3D\n2  TCGA-XF-AAN4  data/gdc/03af9a12-e78e-4d5e-8bfe-eed1168e481c/...       AAN4\n3  TCGA-XF-AAML  data/gdc/05b3d9a5-9414-4aa8-899b-243ee07d7231/...       AAML\n4  TCGA-BT-A2LA  data/gdc/0952e72a-2f04-42fc-bd46-f0378d3f7622/...       A2LA"
  },
  "execution_count": 32,
  "metadata": {},
  "output_type": "execute_result"
 }
]
```

### load clinical data

Here, we will modify the `prep_patient_data` code we used in our earlier script
to include processing the extra parameter: `snv_vcf_paths`.

In practice, we would include this code in a utility script (we typically name
this `analysis/utils/data.py`) so that it can be easily reused, but for this
exercise we will code this manually.

```{.python .input  n=33}
clinical_data = pd.read_csv('data/clinical.csv', sep='|')
clinical_data = clinical_data.merge(vcf_fileinfo, on='patient_id', how='left')
assert clinical_data['snv_vcf_paths'].count()>0
clinical_data.dropna(subset=['snv_vcf_paths'], inplace=True, axis=0)
#clinical_data.head()
```

```{.python .input  n=34}
def prep_patient_data(row):
    # capture key outcome data elements
    patient_id = row['patient_id']
    deceased = row['vital_status'] != 'Alive'
    progressed = row['treatment_outcome_at_tcga_followup'] != 'Complete Response'
    censor_time = float(row['last_contact_days_to'])
    deceased_time = float(row['death_days_to'])
    progressed_time = float(row['new_tumor_event_dx_days_to'])
    
    # compute age at diagnosis
    row['age'] = (-1*row['birth_days_to'])/365.25
    
    # save back in 'row' as-is so that we can see raw values
    row['progressed_time'] = progressed_time
    row['deceased_time'] = deceased_time
    row['censor_time'] = censor_time
    row['progressed'] = progressed
    row['deceased'] = deceased

    # clean up censor time - a number of obs have NaN values 
    if np.isnan(censor_time):
        censor_time = max(progressed_time, deceased_time, censor_time)
    if censor_time > progressed_time:
        censor_time = progressed_time
    if censor_time > deceased_time:
        censor_time = deceased_time

    # save time-to-event-or-censor data elements
    os = deceased_time if deceased else censor_time
    pfs = progressed_time if progressed else os
    
    # again, make sure outcomes aren't NaN
    if np.isnan(os):
        os = censor_time

    if np.isnan(pfs):
        pfs = os
    
    # save transformed versions of outcome back to 'row' object for inspection
    row['pfs'] = pfs
    row['os'] = os
    row['censor_time'] = censor_time
    
    # force progressed time to be < os 
    pfs = min(pfs, os) 
    
    # definition of benefit for this cohort
    benefit = pfs <= 365.25
    
    # these conditions are required by Cohorts
    assert(not np.isnan(pfs))
    assert(not np.isnan(os))
    assert pfs <= os, 'PFS {pfs} is not <= OS {os} for Patient {patid}'.format(pfs=pfs, os=os, patid=patient_id)
    
    # capture snv_vcf_paths, if they exist
    if 'snv_vcf_paths' in row.keys() and isinstance(row['snv_vcf_paths'], str):
        snv_vcf_paths = query_tcga.helpers.convert_to_list(row['snv_vcf_paths'])
    else:
        snv_vcf_paths = None
    
    # create our patient object
    patient = cohorts.Patient(
        id=str(patient_id),
        deceased=deceased,
        progressed=progressed,
        os=os,
        pfs=pfs,
        benefit=benefit,
        additional_data=row,
        snv_vcf_paths=snv_vcf_paths,
    )
    return(patient)

```

As in the earlier example, we will process the rows of our clinical data into a
list of Patients

```{.python .input  n=35}
patients = []
for (i, row) in clinical_data.iterrows():
    patients.append(prep_patient_data(row))
```

And, create a Cohort object from the list

```{.python .input  n=36}
blca_cohort2 = cohorts.Cohort(patients, cache_dir='data-cache2')
```

```{.json .output n=36}
[
 {
  "name": "stdout",
  "output_type": "stream",
  "text": "{'dataframe_hash': 463748958046215442,\n 'provenance_file_summary': {u'cohorts': u'0.2.0+6.g205397b',\n                             u'isovar': u'0.0.6',\n                             u'mhctools': u'0.2.3',\n                             u'numpy': u'1.11.1',\n                             u'pandas': u'0.18.1',\n                             u'pyensembl': u'0.9.7',\n                             u'scipy': u'0.17.1',\n                             u'topiary': u'0.0.21',\n                             u'varcode': u'0.4.14'}}\n"
 }
]
```

## plot survival

As before, we can plot survival for this cohort.

This time, however, instead of passing the name of a field in our dataframe, we
will pass a function which summarizes the variants (VCFs) for this patient.

```{.python .input  n=37}
blca_cohort2.plot_survival(on=cohorts.functions.snv_count)
```

```{.json .output n=37}
[
 {
  "name": "stdout",
  "output_type": "stream",
  "text": "# no condition 5\n# with condition 5\n"
 },
 {
  "data": {
   "text/plain": "<lifelines.StatisticalResult: \nResults\n   df: 1\n   alpha: 0.95\n   t 0: -1\n   test: logrank\n   null distribution: chi squared\n\n   __ p-value ___|__ test statistic __|____ test result ____|__ is significant __\n         0.32222 |              0.980 |  Cannot Reject Null |       False       \n>"
  },
  "execution_count": 37,
  "metadata": {},
  "output_type": "execute_result"
 },
 {
  "data": {
   "image/png": "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\nxCYKEhERsYmCREREbKIgERERmyhIRETEJgoSERGxiYJERERsoiARERGbKEhERMQmChIREbGJgkRE\nRGyiIBEREZsoSERExCYKEhERsYmCREREbKIgERERmyhIRETEJgoSERGxiYJERERsoiARERGbKEhE\nRMQmChIREbGJgkRERGyiIBEREZsoSERExCYKEhERsYmCREREbKIgERERmyhIRETEJgoSERGxiYJE\nRERsoiARERGbKEhERMQmChIREbGJgkRERGyiIBEREZsoSERExCYKEhERsYmCREREbKIgERERmyhI\nRETEJgU+SJYsWULz5s0JDAykc+fO7N+/P79LEhFxKAU6SJYvX05kZCTt2rVjypQpeHp60qtXL86c\nOZPfpYmIOIwCHSRTpkyhc+fOvPbaazz11FNMnz6d4sWLM2/evPwuTUTEYRTYIDl58iR//vknISEh\n5mUuLi40adKEHTt25GNlIiKOpcAGyYkTJzCZTPj6+mZZXrFiRU6fPo1hGPlUmYiIYymwQZKYmAiA\nm5tbluVubm5kZmaSnJycH2WJiDgcl/wu4J+6fcRhMpnuut7J6e4ZmZGRAcC5c+dypzCxC5cTbv03\nLi4ufwsReUDc/s68/R16pwIbJB4eHgAkJSVRokQJ8/KkpCScnZ0pWrToXbe7ePEiABEREblfpOSj\nwrf+Mzc0f8sQecBcvHgx2yWFAhskvr6+GIbB6dOn8fb2Ni+Pi4vDz8/P4nb+/v4sWLCAUqVK4ezs\nnAeViogUfBkZGVy8eBF/f/9s6wpskPj5+VGuXDk2bdpEgwYNAEhLS+P777/PcifXXxUpUoQnn3wy\nr8oUEXlg/PVI5LYCGyQAvXv3ZuTIkXh4eBAUFMT8+fO5du0aL7/8cn6XJiLiMExGAb9Pdt68eXzx\nxRdcvXqVqlWrMnToUAICAvK7LBERh1Hgg0RERPJXgX2ORERE7IOCREREbOJQQaIh57O6du0aVatW\nzfbTv39/c5sZM2YQEhJCzZo16dGjB8eOHcuyj9TUVEaNGkWjRo0ICgrijTfe4MKFC3n9VvLM5s2b\nCQoKyrY8J/opISGBIUOGULduXerUqcOwYcPMIzg8CO7WdwcPHsz2+atWrRpjx441t3HEvsvMzGTu\n3Lm0atWKWrVq0bp1axYsWJCljV195gwHsWzZMqNatWrGtGnTjG3bthm9e/c2ateubcTFxeV3aflm\n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\nRU5ErkJggg==\n",
   "text/plain": "<matplotlib.figure.Figure at 0x116879a50>"
  },
  "metadata": {},
  "output_type": "display_data"
 }
]
```

This is not dissimilar from the `high_grade_tumor` function we created in the
earlier example:

```{.python .input  n=38}
def high_grade_tumor(row):
    return row['tumor_grade'] == 'High Grade'
```

However, this function is a little more intricate.

Let's look, for example, at the `snv_count` function ([as of v0.2.0](https://git
hub.com/hammerlab/cohorts/blob/0.2.0/cohorts/functions.py#L24-L38)):

```{python}
def snv_count(row, cohort, filter_fn=None, normalized_per_mb=None, **kwargs):
    filter_fn = first_not_none_param([filter_fn, cohort.filter_fn], no_filter)
    normalized_per_mb = first_not_none_param([normalized_per_mb,
cohort.normalized_per_mb], False)
    patient_id = row["patient_id"]
    patient_variants = cohort.load_variants(
        patients=[cohort.patient_from_id(patient_id)],
        filter_fn=filter_fn,
        **kwargs)
    if patient_id in patient_variants:
        count = len(patient_variants[patient_id])
        if normalized_per_mb:
            count /= float(get_patient_to_mb(cohort)[patient_id])
        return count
    return np.nan
```

There are essentially three steps to this function, other than processing input
params:
1. Identify the patient, given the row
2. Load variants for this patient (compute them, or load from cache using the
union of VCFs).
    - Optionally apply any QC filters defined at Cohort or function level (more
on this to come)
3. Count the variants associated with this patient

## recreate snv_count for one patient

If we were to recreate these steps semi-manually, they would look something like
the following:

```{.python .input  n=39}
row = blca_cohort2.as_dataframe().iloc[0,:]
filter_fn = cohorts.variant_filters.no_filter
patient_id = row['patient_id']
```

```{.python .input  n=40}
patient_variants = blca_cohort2.load_variants(patients=[blca_cohort2.patient_from_id(patient_id)], filter_fn=filter_fn)
```

```{.python .input  n=41}
count = len(patient_variants[patient_id])
```

```{.python .input  n=42}
count
```

```{.json .output n=42}
[
 {
  "data": {
   "text/plain": "3436"
  },
  "execution_count": 42,
  "metadata": {},
  "output_type": "execute_result"
 }
]
```

## plot snv_count by benefit

```{.python .input  n=43}
blca_cohort2.plot_benefit(cohorts.functions.snv_count)
```

```{.json .output n=43}
[
 {
  "name": "stdout",
  "output_type": "stream",
  "text": "Mann-Whitney test: U=14.0, p-value=0.834531622711 (two-sided)\n"
 },
 {
  "data": {
   "text/plain": "MannWhitneyResults(U=14.0, p_value=0.834531622711, sided_str='two-sided')"
  },
  "execution_count": 43,
  "metadata": {},
  "output_type": "execute_result"
 },
 {
  "data": {
   "image/png": "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\nRK5CYII=\n",
   "text/plain": "<matplotlib.figure.Figure at 0x122cbf590>"
  },
  "metadata": {},
  "output_type": "display_data"
 }
]
```

## retrieve dataframe for a cohort

As before, we can retrieve the same summary statistics as a dataframe.

Here it is possible to provide a list of functions so that multiple summary
statistics can be computed in the same dataframe.

Also note that **the function now returns a tuple** when a list is passed to the
parameter `on`.
    (By comparison, when no `on` parameter is given, `as_dataframe()` simply
returns a dataframe.)

```{.python .input  n=44}
extra_cols, df = blca_cohort2.as_dataframe(on=[cohorts.functions.snv_count])
print(extra_cols)
```

```{.json .output n=44}
[
 {
  "name": "stdout",
  "output_type": "stream",
  "text": "['snv_count']\n"
 }
]
```

## cached effects

There's a fair amount of processing required to load the VCFs into varcode &
summarize patient effects.

To speed up this process, `Cohorts` caches intermediate results of these
computations on the filesystem.

We now have a new folder in our `data-cache` directory:

```{.python .input  n=45}
ls -1 data-cache2/
```

```{.json .output n=45}
[
 {
  "name": "stdout",
  "output_type": "stream",
  "text": "\u001b[34mcached-variants\u001b[m\u001b[m/\r\n"
 }
]
```

Within that directory, we have one folder per patient containing cached
variants:

```{.python .input  n=46}
ls -1 data-cache2/cached-variants/ | head
```

```{.json .output n=46}
[
 {
  "name": "stdout",
  "output_type": "stream",
  "text": "A0S7/\r\nA1AB/\r\nA2LA/\r\nA2OE/\r\nA2PC/\r\nA3I6/\r\nA3IV/\r\nA3MH/\r\nA3SL/\r\nA3WX/\r\n"
 }
]
```

We can also now check our data-source summary to ensure that all cached data
were prepared in the same environment.

```{.python .input  n=47}
blca_cohort2.summarize_data_sources()
```

```{.json .output n=47}
[
 {
  "data": {
   "text/plain": "{'dataframe_hash': 463748958046215442,\n 'provenance_file_summary': {u'cohorts': u'0.2.0+6.g205397b',\n  u'isovar': u'0.0.6',\n  u'mhctools': u'0.2.3',\n  u'numpy': u'1.11.1',\n  u'pandas': u'0.18.1',\n  u'pyensembl': u'0.9.7',\n  u'scipy': u'0.17.1',\n  u'topiary': u'0.0.21',\n  u'varcode': u'0.4.14'}}"
  },
  "execution_count": 47,
  "metadata": {},
  "output_type": "execute_result"
 }
]
```

The next time these data are used, the provenance_file_summary will be updated.

## Summary and next steps

We can now start to see some of the utility of `Cohorts` for managing clinical +
somatic mutation data.

But, there is more to come.

Next up: advanced features (filters) & the dataframe_loader

