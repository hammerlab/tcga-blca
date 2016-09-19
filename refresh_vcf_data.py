from query_tcga import query_tcga as qt

def refresh_vcf_data(project_name):
	qt.download_vcf_data(project_name=project_name, data_dir='data/gdc')

if __name__ == '__main__':
	refresh_vcf_data('TCGA-BLCA')
	print("VCF data downloaded to {}".format('/data/gdc'))
