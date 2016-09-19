from query_tcga import query_tcga as qt

def refresh_wxs_data(project_name):
	qt.download_wxs_data(project_name=project_name, data_dir='data/gdc')

if __name__ == '__main__':
	refresh_wxs_data('TCGA-BLCA')
	print("WXS data downloaded to {}".format('/data/gdc'))
