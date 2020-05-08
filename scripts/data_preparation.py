import yaml
import pandas as pd

stream = open('config.yaml', 'r')
config = yaml.load(stream, Loader=yaml.SafeLoader)
config['project'] = 'PRJNA336352'

runTable_file = config['params']['general'][config['project']]['runTable']['file']
selected_columns = config['params']['general'][config['project']]['runTable']['selected_columns']
aggregate = config['params']['general'][config['project']]['runTable']['aggregate_column']
chip_input_value = config['params']['general'][config['project']]['runTable']['chip_input_value']

runTable = pd.read_csv(runTable_file, sep = ',')

if 'aggregate_column' not in runTable.columns:
   agg = runTable[aggregate].str.split(expand = True)[0]
   agg.replace(to_replace='none', value='Input', inplace=True)
   runTable['aggregate_column'] = agg

runTable.to_csv(runTable_file, sep = ',')    
