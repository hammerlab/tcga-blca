from query_tcga import query_tcga as qt
from query_tcga import config
import logging

config.load_config('config.ini')
logging.basicConfig(level=logging.DEBUG)

def refresh_clin_data(project_name, filename):
	clin = qt.get_clinical_data(project_name=project_name)
	clin.to_csv(filename, sep='|', index=False)
	print("Clinical data written to {file} for project {project}".format(file=filename, project=project_name))

if __name__ == '__main__':
	refresh_clin_data(project_name='TCGA-BLCA', filename='data/clinical.csv')

