from query_tcga import samples

def refresh_wxs_data(project_name):
	samples.download_wxs_data(project_name=project_name, data_dir='data/gdc')

if __name__ == '__main__':
	refresh_wxs_data('TCGA-BLCA')
	print("WXS data downloaded to {}".format('/data/gdc'))
