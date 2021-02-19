Please follow the format as shown in the two sample input datasets.  The formatting is quite strict to that standard.  More work could be done to accept a variety of formatting for datasets.  

Please contact us at <dan.hoeflinger@gmail.com> for help troubleshooting problems.

To run the code with the example datasets run the following:

import GrowthCurveModeler as gcm
gcm.GrowthCurveModeler('Growth Curve Sample Data Set 1.xlsx', threshold=0.2)
gcm.GrowthCurveModeler('Growth Curve Sample Data Set 2.xlsx', threshold=0.2)






To use the new metadata + raw data mode, please use the following:

import GrowthCurveModeler as gcm
gcm.GrowthCurveModeler('raw_data_file.xlsx', MetaDataFile='meta_data_file.xlsx')



Please read the docstring for explanation on how to use the optional arguments.  Example raw and metadata files will be added soon as an example.
