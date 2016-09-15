from query_tcga import query_tcga as qt

def refresh_clin_data(project_name, filename):
	clin = qt.get_clinical_data(project_name=project_name)
	clin.to_csv(filename, sep='|', index=False)

if __name__ == '__main__':
	refresh_clin_data('TCGA-BLCA', filename='data/clinical.csv')
	print("Clinical data written to {}".format('data/clinical.csv'))

