# EVAL
Python scripts for evaluation of repeat prediction. 
## Usage
These scripts require `.fai` files of `seqkit faidx`.
```
seqkit faidx genome.fa
``` 
You can get `genome.fa.fai` by the above command.

### eval_repeat.py
Please prepare repeat annotation files (`.gff`, `.bed`, or `.out` format) for prediction and benchmark, respectively.
```
# python3 eval_repeat.py fai benchmark pred 
python3 eval_repeat.py genome.fa.fai benchmark.fa.out pred.bed
```
You can see the coverage rating of the prediction results. 

### eval_repeat_bulk.py
This script depends on eval_repeat.py. It performs multiple coverage comparisons and outputs them to a single csv file. 
```
echo genome.fa.fai benchmark.fa.out library.fa pred.bed > eval_def.txt
python3 eval_repeat_bulk.py eval_def.txt
```
If you don't have the `library` file, you can replace `library.fa` with any string. 

