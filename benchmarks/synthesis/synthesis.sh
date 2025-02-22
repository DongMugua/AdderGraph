#!/bin/bash
# use bash master.sh not sh master.sh
vivadoOutput="vivadoOutput"
vivadoSynResultsCSV="vivadoSynResults.csv"
flopocoLog="flopocoLog"
powerReport="./tmp/power_report.rpt"
tmp=""
vhdName=""

touch $vivadoOutput
chmod 666 $vivadoOutput
touch $vivadoSynResultsCSV
chmod 666 $vivadoSynResultsCSV
touch $flopocoLog
chmod 666 $flopocoLog
> $vivadoSynResultsCSV # clear file
> $flopocoLog
> $vivadoOutput

echo "Filter;type;dw;cw;method;LUTS;DSPs;data path delay;Total On-Chip Power (W);Device Static (W);Dynamic (W); Clocks (dyn); Logic (dyn); Signals (dyn);i DSPs; I/0 (dyn)" > $vivadoSynResultsCSV
# no need to clear vivadoOutput since it will be overwritten

# read 
filename=$1
while read line; do
	# $line # build one operator using a flopoco call # doesnt work
	# this is the flopoco call, save output of the flopoco call in flopocoLog
	echo $line >> $flopocoLog
	echo $line | bash 2>&1 | tee -a $flopocoLog

	read line
	> $vivadoOutput
	echo $line | bash 2>&1 | tee -a $vivadoOutput 

# formatted output for csv
	# get .vhd name fron python call

	# filter name
	tmp=`echo $line | awk '{print $5}'`
	echo -n $tmp >> $vivadoSynResultsCSV
        echo -n ";" >> $vivadoSynResultsCSV
	vhdName=$(awk '{print $5}' <<< $line)

	# filter type
	#echo -n $(awk -F_ '{print $5}' <<< $TMP) >> $vivadoSynResultsCSV
	echo -n $(awk -F_ '{print $1}' <<< $vhdName) >> $vivadoSynResultsCSV
	echo -n ";" >> $vivadoSynResultsCSV

	# dw, do this TMP and TMP2 thing to avoid printing a new line
	TMP=$(awk -F_ '{print $2}' <<< $vhdName)
	TMP2=`echo -n $TMP | cut -c 3-`
	#echo -n $TMP | cut -c 3- >> $vivadoSynResultsCSV
	#echo -n $(awk -F_ '{print $2}' <<< $vhdName) | cut -c 3- >> $vivadoSynResultsCSV
	echo -n $TMP2 >> $vivadoSynResultsCSV
	echo -n ";" >> $vivadoSynResultsCSV

	#cw
	TMP=$(awk -F_ '{print $3}' <<< $vhdName)
	TMP2=`echo -n $TMP | cut -c 3-`
	#echo -n $TMP | cut -c 3- >> $vivadoSynResultsCSV
	#echo -n $(awk -F_ '{print $3}' <<< $vhdName) | cut -c 3- >> $vivadoSynResultsCSV
        echo -n $TMP2 >> $vivadoSynResultsCSV
	echo -n ";" >> $vivadoSynResultsCSV

	# method
	METHOD=$(awk -F_ '{print $4}' <<< $vhdName) 
	METHOD=$(awk -F. '{print $1}' <<< $METHOD)
	echo -n $METHOD >> $vivadoSynResultsCSV
        echo -n ";" >> $vivadoSynResultsCSV

	# LUTs
	tmp=`grep "Slice LUTs" $vivadoOutput | awk '{print $5}'`
	echo -n $tmp >> $vivadoSynResultsCSV
	echo -n ";" >> $vivadoSynResultsCSV
	# DSPs
	tmp=`grep "DSPs  " $vivadoOutput | awk '{print $4}'`
	echo -n $tmp >> $vivadoSynResultsCSV
        echo -n ";" >> $vivadoSynResultsCSV
	# Delay
	# FixIIR has other output so grep for output delay instead of Data path delay
	#if [[ "$vhdName" == "fixIIR.vhd" ]]; then
	#	tmp=`grep "output delay" $vivadoOutput | awk '{print $4}'`
	#else
		tmp=`grep "Data Path Delay" $vivadoOutput | awk '{print $4}'`
#	fi
	echo -n $tmp >> $vivadoSynResultsCSV
        echo -n ";" >> $vivadoSynResultsCSV

	# total power
	tmp=`grep "Total On-Chip Power (W)" $powerReport | awk '{print $7}'`
        echo -n $tmp >> $vivadoSynResultsCSV
        echo -n ";" >> $vivadoSynResultsCSV

	# static power
	tmp=`grep "Device Static (W)" $powerReport | awk '{print $6}'`
        echo -n $tmp >> $vivadoSynResultsCSV
        echo -n ";" >> $vivadoSynResultsCSV

	#dynamic power
        tmp=`grep "Dynamic (W)" $powerReport | awk '{print $5}'`
        echo -n $tmp >> $vivadoSynResultsCSV
        echo -n ";" >> $vivadoSynResultsCSV

	#clocks
	tmp=`grep "Clocks" $powerReport | awk '{print $4}'`
        echo -n $tmp >> $vivadoSynResultsCSV
        echo -n ";" >> $vivadoSynResultsCSV

	#slice logic
	tmp=`grep "Slice Logic" $powerReport | awk '{print $5}'`
        echo -n $tmp >> $vivadoSynResultsCSV
        echo -n ";" >> $vivadoSynResultsCSV

	#signals
	tmp=`grep "Signals" $powerReport | awk '{print $4}'`
        echo -n $tmp >> $vivadoSynResultsCSV
        echo -n ";" >> $vivadoSynResultsCSV

	#DSPs
        tmp=`grep "DSPs" $powerReport | awk '{print $4}'`
        echo -n $tmp >> $vivadoSynResultsCSV
        echo -n ";" >> $vivadoSynResultsCSV


	#I/O
        tmp=`grep "I/O            |" $powerReport | awk '{print $4}'`
        #echo -n $tmp >> $vivadoSynResultsCSV
        #echo -n ";" >> $vivadoSynResultsCSV

        echo -n $tmp >> $vivadoSynResultsCSV
	# last one with line break
        echo ";" >> $vivadoSynResultsCSV
	
	# delte tmp folder
	rm -r tmp/
done < $filename
