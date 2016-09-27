from query_tcga import samples, config
import numpy as np
import logging

config.load_config('config.ini')
logging.basicConfig(level=logging.DEBUG)

def refresh_vcf_data(project_name, data_dir, filename):
	vcf_files = samples.download_vcf_files(project_name=project_name, data_dir=data_dir)
	logging.info("VCF data downloaded to {}".format(data_dir))

	## check for non-uniqueness in reference names
	vcf_file_summary = vcf_files.fileinfo
	reference_names = np.unique(vcf_file_summary['reference_name'])
	if len(reference_names)==1:
		logging.info('All VCFs have reference genome: {}'.format(reference_names[0]))
	else:
		reference_summary = vcf_file_summary.groupby('reference_name').agg(len)
		logging.warn('Not all VCFs have the same inferred reference! \n {}'.format(reference_summary))
	vcf_file_summary.to_csv(filename, sep='|', index=False)
	print("VCF fileinfo written to {file} for project {project}".format(file=filename, project=project_name))


if __name__ == '__main__':
	refresh_vcf_data('TCGA-BLCA', data_dir='data/gdc', filename='data/vcf_fileinfo.csv')
