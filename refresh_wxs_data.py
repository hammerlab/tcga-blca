from query_tcga import samples, config
import logging

config.load_config('config.ini')
logging.basicConfig(level=logging.DEBUG)

def refresh_wxs_data(project_name, filename, data_dir):
	wxs_files = samples.download_wxs_files(project_name=project_name, data_dir=data_dir)
	wxs_files.fileinfo.to_csv(filename, sep='|', index=False)
	print("WXS fileinfo written to {file} for project {project}".format(file=filename, project=project_name))


if __name__ == '__main__':
	refresh_wxs_data('TCGA-BLCA', filename='data/wxs_fileinfo.csv', data_dir=config.get_setting_value('GDC_DATA_DIR'))
