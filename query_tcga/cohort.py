import cohorts
from query_tcga import query_tcga as qt
import numpy as np

def prep_patient_data(row):
    patient_id = row['case_id']
    deceased = row['vital_status'] != 'Alive'
    progressed = row['treatment_outcome_at_tcga_followup'] != 'Complete Response'
    censor_time = float(row['last_contact_days_to'])
    deceased_time = float(row['death_days_to'])
    progressed_time = min(float(row['new_tumor_event_dx_days_to']),
                         float(row['new_tumor_event_surgery_days_to']),
                         float(row['new_tumor_event_surgery_days_to_met']),
                         )

    row['progressed_time'] = progressed_time
    row['deceased_time'] = deceased_time
    row['censor_time'] = censor_time
    row['progressed'] = progressed
    row['deceased'] = deceased

    censor_time = censor_time or 0
    if np.isnan(censor_time):
        censor_time = max(progressed_time, deceased_time, censor_time)

    if censor_time > progressed_time:
        censor_time = progressed_time
    if censor_time > deceased_time:
        censor_time = deceased_time

    os = deceased_time if deceased else censor_time
    pfs = progressed_time if progressed else os

    if np.isnan(os):
        os = censor_time

    if np.isnan(pfs):
        pfs = os

    row['pfs'] = pfs
    row['os'] = os
    row['censor_time'] = censor_time

    pfs = min(pfs, os) ## force progressed time to be < os 

    benefit = pfs <= 365.25

    assert(not np.isnan(pfs))
    assert(not np.isnan(os))
    assert pfs <= os, 'PFS {pfs} is not <= OS {os} for Patient {patid}'.format(pfs=pfs, os=os, patid=patient_id)

    patient = cohorts.Patient(
        id=patient_id,
        deceased=deceased,
        progressed=progressed,
        os=os,
        pfs=pfs,
        benefit=benefit,
        additional_data=row,
    )
    return(patient)


def build_cohort(project_name, **kwargs):

    clin = qt.get_clinical_data(project_name=project_name, **kwargs)


    patients = []
    for (i, row) in clin.iterrows():
        patients.append(prep_patient_data(row))

    cohort = cohorts.Cohort(
            patients=patients,
            cache_dir='data-cache'
        )
    return(cohort)


