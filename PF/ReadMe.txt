%% The power flow code used here is downloaded from 
https://www.mathworks.com/matlabcentral/fileexchange/56074-linear-load-flow-in-power-distribution-systems-unbalanced-case?s_tid=prof_contriblnk


%% This folder uses the *xlsx file as the input, runs the three-phase unbalanced power flow, and save the feeder data and power flow results as *mat file, which will be used as the input data in DSE Code.



%% Simply run the following command in MATLAB command window. 

>> SEData('IEEE37')
>> SEData('IEEE123')