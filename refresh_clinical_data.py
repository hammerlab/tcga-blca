from query_tcga import query_tcga as qt
from query_tcga import config
import logging

config.set_value(USE_CACHE=True)
config.set_value(GDC_TOKEN_PATH='/Users/jacquelineburos/Downloads/gdc-user-token.2016-09-26T12-23-27-04-00.txt')
config.set_value(GDC_DATA_DIR='/Users/jacquelineburos/projects/tcga-blca/data/gdc/')

logging.basicConfig(level=logging.DEBUG)

def refresh_clin_data(project_name, filename):
	clin = qt.get_clinical_data(project_name=project_name)
	clin.to_csv(filename, sep='|', index=False)

if __name__ == '__main__':
	refresh_clin_data(project_name='TCGA-BLCA', filename='data/clinical.csv')
	print("Clinical data written to {file} for project {project}".format(file='data/clinical.csv', project='TCGA-BLCA'))

