import pandas as pd
import cohorts

def init_cohort(project_name, data_dir):
    ## start with clinical data, construct list of Patients
    patients = list()
    df = pd.read_csv('data/clinical.csv', sep='|')
    for i, row in df.iterrows():
        days_to_death = row['death_days_to']
        is_deceased = row['vital_status'] == 'Dead'
        days_to_progression = min(row['new_tumor_event_dx_days_to'], row['death_days_to'])
        is_progressed = row['treatment_outcome_at_tcga_followup'] == 'Progressive Disease'
        is_progressed_or_deceased = is_deceased or is_progressed
        patient = cohorts.Patient(id=row['patient_id'],
                                  os=days_to_death,
                                  pfs=days_to_progression, 
                                  deceased=is_deceased,
                                  progressed=is_progressed,
                                  progressed_or_deceased=is_progressed_or_deceased,
                                  #hla_alleles=row["hla_allele_list"],
                                  #snv_vcf_paths=snv_vcf_paths,
                                  #indel_vcf_paths=indel_vcf_paths,
                                  #normal_sample=normal_sample,
                                  #tumor_sample=tumor_sample,
                                  additional_data=row.to_dict())
        patients.append(patient)

    cohort = cohorts.Cohort(patients)


    return cohort


