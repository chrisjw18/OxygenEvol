#make generic calibration data available to OxygenEvol package

generic_cali<-read.csv('PSI_F3500_Calibration.csv')
devtools::use_data(generic_cali)
