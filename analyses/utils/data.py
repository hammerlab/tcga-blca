import pandas as pd
import numpy as np
import cohorts
import os
import query_tcga
from query_tcga import helpers

def _load_clinical_data(filepath='data/clinical.csv'):
    clinical_data = pd.read_csv(filepath, sep='|')
    return clinical_data


def _load_vcf_fileinfo(filepath='data/vcf_fileinfo.csv'):
    all_vcf_files = pd.read_csv(filepath, sep='|')
    vcf_fileinfo = all_vcf_files.loc[:,['submitter_id','filepath']]
    vcf_fileinfo.rename(columns = {'filepath': 'snv_vcf_paths'}, inplace=True)
    vcf_fileinfo['patient_id'] = vcf_fileinfo['submitter_id'].apply(lambda x: x.split('-')[2])
    vcf_fileinfo['snv_vcf_paths'] = vcf_fileinfo['snv_vcf_paths'].apply(query_tcga.helpers.convert_to_list)
    vcf_fileinfo_agg = vcf_fileinfo.groupby('patient_id').agg({'snv_vcf_paths': 'sum'}).reset_index()
    return vcf_fileinfo_agg


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
    if 'snv_vcf_paths' in row.keys() and isinstance(row['snv_vcf_paths'], list):
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
        additional_data=row.to_dict(),
        snv_vcf_paths=snv_vcf_paths,
        #hla_alleles=row["hla_allele_list"],
        #indel_vcf_paths=indel_vcf_paths,
        #normal_sample=normal_sample,
        #tumor_sample=tumor_sample,
    )
    return(patient)


def init_cohort(project_data_dir='data', cache_dir='data-cache', **kwargs):
    """ Construct cohort from data downloaded from TCGA using `query_tcga`

        Assumes that clinical & VCF files have been downloaded to GDC_DATA_DIR, and 
        that summary files have been stored in `project_data_dir`.

        cache_dir & **kwargs are passed to cohorts.Cohort()/ See ?cohorts.Cohort for details

    """
    ## start with clinical data, construct list of Patients
    vcf_fileinfo = _load_vcf_fileinfo(filepath=os.path.join(project_data_dir, 'vcf_fileinfo.csv'))
    clinical_data = _load_clinical_data(filepath=os.path.join(project_data_dir, 'clinical.csv'))

    # merge clinical & vcf data
    clinical_data = clinical_data.merge(vcf_fileinfo, on='patient_id', how='left')
    assert clinical_data['snv_vcf_paths'].count()>0
    clinical_data.dropna(subset=['snv_vcf_paths'], inplace=True, axis=0)
    assert clinical_data.duplicated('patient_id').any() == False, 'Duplicates by patient_id'

    patients = []
    for (i, row) in clinical_data.iterrows():
        patients.append(prep_patient_data(row))

    cohort = cohorts.Cohort(patients, cache_dir=cache_dir, **kwargs)

    return cohort


