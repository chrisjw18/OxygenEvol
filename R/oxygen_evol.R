#OxygenEvol

##NB: This relies on time stamps of measurements made on the FL3500
#     and the FireSting for matching data, thus both need to be run
#     on the same computer, or computers whose date/time are in sync.

# - read in FL3500 output (.txt)
# - extract fluoroescence measurement times
# - extract the PAR values for each light step from FL3500 data
# - read in FireSting output (.txt)
# - plot whole oxygen trace showing timing of F measurements
# - extract oxygen traces per light step, fit lm, take slope as rate of O2 evolution
# - combine data into data frame
# - make plots of Oxygen evolution ~ light step and ~ PAR (Scatter plot)

#Return:

## License information
#
# OxygenEvol - Data processing from the FireSting optode in R
#     Copyright (C) 2018 by Bruno Jesus, Douglas A. Campbell, Christopher J. Williamson.
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact Chris Williamson at c.williamson@bristol.ac.uk

#'@title oxygen_evol.
#'@description extract and match firesting oxygen data to FL3500 PSI fluorescence data.
#'@param fluorwin_filename name of the fluorwin .txt file to process.
#'@param firesting_filename name of the firesting oxygen data file to process.
#'@param no_light_steps the number of light steps in the curve.
#'@param time_step length of the window (in seconds) of oxygen data to extract per
#'    light step.
#'@param calibration_file flashlet calibration file used to calculate PAR from
#'    voltages for FL3500 data. Generic calibration file is built into package
#'    so leave as NA if no specific calibration available.
#'@param data_output logical (default FALSE), if TRUE dataframes produced are
#'    written to CSV and plots written to single .pdf file (see Values below)
#'@param end_respiration logical (default FALSE), if TRUE function will perform
#'    an end respiration calculation using oxygen data from immediately after the
#'    final light step of the light curve over the duration set by time_step
#' @return oxy_results list containing the following two dataframes:
#'   1. raw_oxygen dataframe containing the tidy outputs from the Firesting
#'   optode (i.e. complete oxygen trace with formatted time columns).
#'   2. matched_oxygen dataframe containing the FL3500 light step with matched
#'   calculated rate of oxygen evolution (units of umol s-1), with columns
#'   showing voltage and PAR levels per light step
#'
#'   If data_output = TRUE:
#'   x2 CSV files of raw_oxygen and matched_oxygen datasets appended with oxygen
#'   datafile name
#'
#'   AND
#'
#'   .pdf document appended with oxygen datafile name showing:
#'   1. complete oxygen trace with grey areas highlighting which sections have
#'   been clipped and matched to fluorescence data (As determined by the
#'   time_step selected)
#'   2. individual plots of each linear regression for oxygen evolution
#'   calculations
#'   3. calculated oxygen evolution ~ light step over the course of the light
#'   curve with PAR levels shown in blue trace in brackground
#'   4. scatter plot of oxygen evolution ~ PAR
#'@examples To import example data, run the following code:
#'
#'example_fluorwin_data<-system.file("extdata", "example_fluorwin_data.txt", package = "OxygenEvol")
#'example_oxygen_data<-system.file("extdata", "example_oxygen_data.txt", package = "OxygenEvol")
#'oxygen_evol(example_fluorwin_data, example_oxygen_data, no_light_steps=11, time_step=60, calibration_file=NA, data_output=F, end_respiration = T)
#'@export

oxygen_evol<-function(fluorwin_filename, firesting_filename, no_light_steps, time_step=60, calibration_file=NA, data_output=F, end_respiration=F){


  ###############read in fluorwin .txt output to extract timings
  #of measurements and format them
  f.times<-read.csv(fluorwin_filename, skip=9, sep='\t', header=F)[1:no_light_steps,c(1,3)]
  a<-grepl("AM|PM", f.times[,2])
  if (a[1]==FALSE){
    f.times[,2] <- as.POSIXct(f.times[,2], format="%d/%m/%Y %H:%M:%S")
  }else{
    f.times[,2] <- as.POSIXct(f.times[,2], format="%m/%d/%Y %I:%M:%S %p")
  }
  names(f.times)<-c('dataset', 'meas_time')

  #add catch for end dark respiration assessment
  if(end_respiration==T){
    end.resp<-data.frame(dataset='end_resp_measurement', meas_time=f.times$meas_time[nrow(f.times)]+time_step)
    f.times<-rbind(f.times, end.resp)
  }


  #calculate time window for viewing oxygen based on user input
  #this is the time window backwards from the time of F measurements
  #over which we will regress oxygen evolution per light step
  f.times$start.times<-f.times$meas_time - time_step

  #######add in blank columns to house regression coefficients
  f.times$light_step<-1:nrow(f.times)
  f.times$slope<-NA
  f.times$r.sq<-NA
  f.times$p.val<-NA

  ########get the PAR values for each light step in the FL3500 output
  #generate the meta.data
  data.raw<-scan(fluorwin_filename,character(0),sep='\n',quiet=TRUE);
  data.sep <- which(data.raw=="-----------------------------------");
  # Extract Meta-Data
  meta.data.names<-unlist(strsplit(data.raw[data.sep[2]+2],split='\t'),use.names=FALSE);
  tst<-unlist(strsplit(data.raw[data.sep[2]+3],split='\t'));
  meta.data<-as.data.frame(as.list(as.numeric(tst[2:length(tst)])));
  names(meta.data)<-meta.data.names[2:length(meta.data.names)];


  #use generic calibration table or read in calibration if available
  if(is.na(calibration_file)){
    cali<-generic_cali
  } else{
    cali<-read.csv(calibration_file)
  }


  #generate the voltages table and look up PAR using calibration
  blue<-meta.data[1,"Blue_light_curves_enable"]
  red<-meta.data[1,"Red_light_curves_enable"]

  if (blue==1){
      lc_voltage<-names(meta.data)[which(grepl('Blue',names(meta.data)))]
      lc_voltage<-lc_voltage[which(lc_voltage!='Blue_light_curves_enable')]
      lc_voltage<-lc_voltage[which(lc_voltage!='Blue_light_curves')]

      voltages<-data.frame(step=rep(NA,no_light_steps), voltage=NA, par=NA)
      for(i in 1:(no_light_steps)){
        voltages$step[i]<-lc_voltage[which(lc_voltage==paste('Blue_light_curves',i-1,sep='_'))]
        voltages$voltage[i]<-meta.data[1, voltages$step[i]]
        #look up PAR using voltage in calibration file - must already be loaded under 'cali'
        voltages$par[i]<-cali[cali$Power.Level==voltages$voltage[i], paste('Blue','.Actinic',sep='')]
        #print(voltages)
      }
  } else {
      lc_voltage<-names(meta.data)[which(grepl('Red',names(meta.data)))]
      lc_voltage<-lc_voltage[which(lc_voltage!='Red_light_curves_enable')]
      lc_voltage<-lc_voltage[which(lc_voltage!='Red_light_curves')]

      voltages<-data.frame(step=rep(NA,no_light_steps), voltage=NA, par=NA)
      for(i in 1:(no_light_steps)){
        voltages$step[i]<-lc_voltage[which(lc_voltage==paste('Red_light_curves',i-1,sep='_'))]
        voltages$voltage[i]<-meta.data[1, voltages$step[i]]
        #look up PAR using voltage in calibration file - must already be loaded under 'cali'
        voltages$par[i]<-cali[cali$Power.Level==voltages$voltage[i], paste('Red','.Actinic',sep='')]
      }
  }#end of else

  #input into f.timees dataframe
  if(end_respiration==F){
  f.times$voltage<-voltages$voltage
  f.times$par<-voltages$par
  } else {
    volt_append<-data.frame(step='end_resp',voltage=0, par=0)
    voltages<-rbind(voltages, volt_append)
    f.times$voltage<-voltages$voltage
    f.times$par<-voltages$par
  }

  #add a column of fluorwin output data filename into f.times for traceability.
  f.nam<-strsplit(fluorwin_filename, split='[.]')[[1]][1]
  f.times$fluorwin_filename<-f.nam

  ###############read in optode data, tidy up data frame, match to F data
  m1<-read.table(firesting_filename,skip=13,header=T,sep="\t")

  #replace names for something more obvious
  names(m1)<-c("date","time","time_s","comment","ch1_O2","ch2_O2","ch3_O2","ch4_O2",
                "ch1_temp","ch2_temp","ch3_temp","ch4_temp","pressure","humidity",
                "temp_probe","int_temp","analogue_in","ch1_raw","ch2_raw","ch3_raw"
                ,"ch4_raw",
                "ch1_signal","ch2_signal","ch3_signal","ch4_signal",
                "ch1_light","ch2_light","ch3_light","ch4_light")
  #remove last column
  m1<-m1[,-30]

  #establish an oxygen time (o.time) column in m1 with formatted date/time
  m1$o.time<-as.POSIXct(paste(m1$date,m1$time,sep=' '), format="%m/%d/%Y
              %H:%M:%S")

  #loop through light steps and get slope characteristics of the oxygen data
  #using the fluoresncence times
  for(i in 1:nrow(f.times)){
    x0<-f.times$start.times[i]
    x1<-f.times[i,2]

    y<-m1$ch1_O2[which(m1$o.time==x0):which(m1$o.time==x1)]
    x<-m1$time_s[which(m1$o.time==x0):which(m1$o.time==x1)]

    my.lm<-lm(y~x)

    f.times$slope[i]<-summary(my.lm)$coefficients[2,1] #this is umol O2 per sec
    f.times$r.sq[i]<-round(summary(my.lm)$r.squared,3)
    f.times$p.val[i]<-summary(my.lm)$coefficients[2,4]
  }#end of i loop


  #make plots:
  if(data_output==T){
  #initiate pdf to catch all plots
  nam<-strsplit(firesting_filename, split='[.]')[[1]][1]
  pdf(file=paste(nam, '_oxygen_evolution.pdf',sep=''))

  #main plot of oxygen trace showing time intervals for oxygen evolution

  #regression
  par(mar=c(4,4,1,1), mgp=c(2.4,0.4,0), las=1, tck=-0.01)
  plot(m1$o.time, m1$ch1_O2, type='l', ylab=expression(paste(O[2],' (',mu,'mol ',l^-1,')')), xlab='Time')
  rect(f.times$start.times, 0, f.times[,2], 400, col='lightgrey', xpd=F,border=NA)
  points(m1$o.time, m1$ch1_O2, type='l')
  axis(1, labels=F, lwd.ticks = 0)
  axis(3, labels=F, lwd.ticks = 0)
  text(f.times$meas_time-((f.times$meas_time-f.times$start)/2), max(m1$ch1_O2), labels=1:nrow(f.times))

  if(end_respiration==T){
    arrows(f.times$start.times[nrow(f.times)], 0,f.times$start.times[nrow(f.times)],400, length=0, lty=3)
  }
  #make individual regression plots

  for(i in 1:nrow(f.times)){
    x0<-f.times$start.times[i]
    x1<-f.times[i,2]

    y<-m1$ch1_O2[which(m1$o.time==x0):which(m1$o.time==x1)]
    x<-m1$time_s[which(m1$o.time==x0):which(m1$o.time==x1)]

    my.lm<-lm(y~x)

    slope<-round(summary(my.lm)$coefficients[2,1]*60,3) #this is umol O2 per minute.
    r.sq<-round(summary(my.lm)$r.squared,3)
    p.val<-round(summary(my.lm)$coefficients[2,4],3)

    par(mar=c(3,3,1,1), xpd=NA,tck=-0.01,las=1,mgp=c(2.2,0.5,0))
    plot(x,y, ylab=expression(paste(mu,'mol ',O[2],' ',l^-1)), xlab='Time (s)', main=paste('Light step',i))
    abline(my.lm, xpd=F, col='red', lty=3)

    legend('topright', legend=c(paste('slope =',slope), paste('R2 =',round(r.sq,3)), paste('p.val =',round(p.val, 3) )), bty='n')

  }#end of i loop


  #make plot of the oxygen data (slope) ~ light step, with PAR in background as
  #blue trace
  if(end_respiration==T){
    my.xlim=c(0,no_light_steps+1)
    col<-c(rep('black',no_light_steps),'red')
  } else {
    my.xlim=c(0,no_light_steps)
    col<-rep('black',no_light_steps)
  }
  par(mar=c(3,3,1,3), mgp=c(1.6,0.4,0), tck=-0.01, las=1, xpd=F)
  plot(f.times$light_step, f.times$slope, type='o', pch=19, ylab=expression(paste('Oxygen Evolution ( ',mu,'mol ',O[2],' ',l^-1,' ',s^-1,')')), xlab='Light Step', xlim=my.xlim, col=col)
  arrows(0,0,nrow(f.times)+1,0,length=0, lty=3)
  par(new=T)
  plot(f.times$light_step+0.11, f.times$par, type='S', lwd=4, col=addTrans('blue',50),yaxt='n', ylab='',xlab='', xaxt='n', cex.axis=0.8, xlim=my.xlim)
  axis(4, cex.axis=0.8)
  text(nrow(f.times)*1.11, (max(f.times$par)-(max(f.times$par)/2)), labels=expression(paste('PAR (',mu,'mol photons ', m^-2,' ',s^-1,')')),xpd=T, srt=-90, cex=1)
  if(end_respiration==T){
    legend('topleft', pch=c(19,NA,19), col=c('black',addTrans('blue',50),'red'), lty=1, legend=c('Oxygen Evol', 'PAR','End Respiration'), lwd=c(1,3), bty='n')
  } else{
  legend('topleft', pch=c(19,NA), col=c('black',addTrans('blue',50)), lty=1, legend=c('Oxygen Evol', 'PAR'), lwd=c(1,3), bty='n')
  }


  #plot of oxygen evolution (slope) ~ PAR (NB this is a scatter plot, so if using
  #the same level of PAR in the light curve more than once, this is not linear with
  #the progression of the light curve)
  plot(f.times$par, f.times$slope, ylab=expression(paste('Oxygen Evolution ( ',mu,'mol ',O[2],' ',l^-1,' ',s^-1,')')), xlab=expression(paste('PAR (',mu,'mol photons ',m^-2,' ',s^-1,')')), pch=19)

  #turn off pdf catcher.
  dev.off()
  }

  #print the main matched oxygen plot to the screen regarless of data_output
  #could be improved with a dev.copy command to reduce the code
  #make plot of the oxygen data (slope) ~ light step, with PAR in background as
  #blue trace
  if(end_respiration==T){
    my.xlim=c(0,no_light_steps+1)
    col<-c(rep('black',no_light_steps),'red')
  } else {
    my.xlim=c(0,no_light_steps)
    col<-rep('black',no_light_steps)
  }
  par(mar=c(3,3,1,3), mgp=c(1.6,0.4,0), tck=-0.01, las=1, xpd=F)
  plot(f.times$light_step, f.times$slope, type='o', pch=19, ylab=expression(paste('Oxygen Evolution ( ',mu,'mol ',O[2],' ',l^-1,' ',s^-1,')')), xlab='Light Step', xlim=my.xlim, col=col)
  arrows(0,0,nrow(f.times)+1,0,length=0, lty=3)
  par(new=T)
  plot(f.times$light_step+0.11, f.times$par, type='S', lwd=4, col=addTrans('blue',50),yaxt='n', ylab='',xlab='', xaxt='n', cex.axis=0.8, xlim=my.xlim)
  axis(4, cex.axis=0.8)
  text(nrow(f.times)*1.11, (max(f.times$par)-(max(f.times$par)/2)), labels=expression(paste('PAR (',mu,'mol photons ', m^-2,' ',s^-1,')')),xpd=T, srt=-90, cex=1)
  if(end_respiration==T){
    legend('topleft', pch=c(19,NA,19), col=c('black',addTrans('blue',50),'red'), lty=1, legend=c('Oxygen Evol', 'PAR','End Respiration'), lwd=c(1,3), bty='n')
  } else{
    legend('topleft', pch=c(19,NA), col=c('black',addTrans('blue',50)), lty=1, legend=c('Oxygen Evol', 'PAR'), lwd=c(1,3), bty='n')
  }


  #returns list including raw oxygen (processed and tidied fire sting data) and
  #matched_oxygen (FL3500 measurements matched to oxygen slope)
  oxy_results<-list()
  oxy_results$raw_oxygen<-m1
  oxy_results$matched_oxygen<-f.times
  oxy_results<<-oxy_results
  return(oxy_results)


  if(data_output==T){
    nam<-strsplit(firesting_filename, split='[.]')[[1]][1]
    write.csv(m1, file=paste(nam, '_raw_oxygen.csv',sep=''))
    write.csv(f.times, file=paste(nam, '_matched_oxygen.csv',sep=''))
  }


}#end of oxygen_evol function












