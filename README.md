# REPrise
Download highest version number from 
```
https://github.com/hmdlab/REPrise
```

This assumes you have a C++ compiler. 
```
make
```

# Command and Options
```
Usage

REPrise [-input genome file] [-output outputname] [Options]

Options
(Required)
   -input  STR         input file name. You can input assembled genome file, or hard masked genome file
   -output STR         output file name. REPrise outputs STR.freq, STR.bed STR.masked and STR.reprof (consensus seqnences)

(Optional)
   -h                  Print help and exit
   -v                  Verbose

   -match INT          Match score of the extension alignment (default = 1)
   -match INT          Mismatch score of the extension alignment (default = -1)
   -gap   INT          Gap open score of the extension alignment (default = -5)
   -gapex  INT         Gap extension score of the extension alignment (default = -1)
   -capplenalty INT    Penalty of the imcomplete length alignment (default = -20)
   -dist INT           Number of mismatches allowed in inexact seed (default = 0)

   -maxextend INT      Uppler limit length of extension in one side direction of consensus repeat (default = 10000)
   -maxrepeat INT      Maximum Number of elements belonging to one repeat family (default = 100000)
   -maxgap INT         Band size(= maximum number of gaps allowed) of extension alignment (default = 5)
   -stopafter INT      If the maximum score of extension alignment does not change INT consecutive times, that alignment will stop (default = 100)
   -minlength INT      Minimum number of length of the consensus sequence of repeat family(default = 50)
   -minfreq INT        Minimum number of elements  belonging to one repeat family (default = 3)
   -minimprovement INT Penalty associated with the number of regions to be extended as the repeat regions (default = 3)
   -tandemdist INT     Interval to match the same seed to avoid seed matching with tandem repeats(default = 500)

   -pa INT             Number of openMP parallel cores
```
# License

# Changelogs
2023/01/XX Vesion 1.0.0 relased. 
# Reference
