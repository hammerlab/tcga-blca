from query_tcga import query_tcga as qt
import logging

logging.basicConfig(level=logging.INFO)

def refresh_clin_data(project_name, filename):
	clin = qt.get_clinical_data(project_name=project_name)
	clin.to_csv(filename, sep='|', index=False)

if __name__ == '__main__':
	refresh_clin_data(project_name='TCGA-BLCA', filename='data/clinical.csv')
	print("Clinical data written to {file} for project {project}".format(file='data/clinical.csv', project='TCGA-BLCA'))

