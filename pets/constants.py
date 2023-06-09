import os
COMMON_OPTPARSE = os.path.abspath(os.path.join(ROOT_PATH, '..', 'pets', 'common_optparse.py'))
REPORT_FOLDER = os.path.abspath(os.path.join(ROOT_PATH, '..', 'templates'))
EXTERNAL_DATA = os.path.abspath(os.path.join(ROOT_PATH, '..', 'external_data'))
EXTERNAL_CODE = os.path.abspath(os.path.join(ROOT_PATH, '..', 'external_code'))
HPO_FILE = os.path.join(EXTERNAL_DATA, 'hp.json')
MONDO_FILE = os.path.join(EXTERNAL_DATA, 'mondo.obo')
IC_FILE = os.path.join(EXTERNAL_DATA, 'uniq_hpo_with_CI.txt')