Scanner.py is a program for numerical calculations of plasma parameters from mathematical model developed by Golubev M.S.

When run from terminal it accepts up to 3 optional positional arguments:
1. config file name. Config must be within \conf folder in the root of the program. When empty defaults to conf.json5
2. output file name (without extension). Output .xlsx file is placed in \output upon calculation completion. When empty defaults to "parameters"
3. filter file name. Filter must be within \conf folder in the root of the program. When empty defaults to False.

Config file specified in confname MUST abide by formatting:
{
	"Variable name1": Value1,
	"Variable name2": Value2
}

Passed variables can be either float or [range]. e.g.: "p": [1,2,3],
When variable passed as [range] the scanner will iteratively scan through every combination of floats in passed ranges 
	and add the result to \output\%filename%_%Y%-%m%-%d%_%H%-%M%-%S%.xlsx

Non-optional parameters of calc.json5 input are:
{
	'gas': %, # atmosphere of plasma discharge: He, Ne or Ar.
	'l': %, # length of discharge gap in meters.
	'R': %, # radius of circular cathode in meters.
	'f': %, # frequency of RF volatage in MHz, after specification this parameter will automatically be multiplied by 1e6.
	'Assy': %, # Asymmetry coefficient
	'p': %, # Pressure within discharge cahamber.
	*'Pwr': %, # Power that will be absorbed by plasma in Watts. 
		*Can be replaced by Vrf.
	*'Vrf': %, # Amplitude of RF voltage on the electrodes. 
		*Can replace Pwr and override Vrf calculation. 	
			If both Vrf & Pwr are specified, Vrf takes priority.
}

Filter file specified must abide by formatting:
['Column1','Column2','Column3']
And contain valid column keys. For valid values of keys look up Sample.xlsx in \output.
If given filtername, output file will contain ONLY column specified within filtername.
Be careful to include all nessesary columns within filter, othwerwise those values won't be saved in output.