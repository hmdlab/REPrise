import os
import subprocess
from argparse import ArgumentParser
from eval_repeat import read_fai, read_annot, evaluate

columns = [
    'FASTA index file',
    'Ground truth file',
    'Detected family file',
    'Prediction file',
    '#Family',
    'Sensitivity',
    'Specificity',
    '%Positive',
    'PNR',
    'Precision',
    'F-score'
]

def main(args):
    print(','.join(columns))

    with open(args.def_filename, 'r') as f:
        line = f.readline()
        while line:
            fai, gt_filename, seq_filename, pred_filename = line.strip().split(' ')

            seq_size, seq_names = read_fai(fai)
            gt = read_annot(gt_filename, seq_names)

            if os.path.exists(seq_filename):
                stat = subprocess.check_output('seqkit stat -T'.split() + [seq_filename])
                num_families = stat.decode('utf-8').split('\n')[1].split('\t')[3]
            else:
                num_families = '-'

            if os.path.exists(pred_filename):
                pred = read_annot(pred_filename)
                evaluations = ['{:.2f}'.format(x * 100) for x in evaluate(seq_size, seq_names, gt, pred)]
            else:
                evaluations = ['-'] * 6
            
            print(','.join([fai, gt_filename, seq_filename, pred_filename, num_families] + evaluations))

            line = f.readline()

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('def_filename')

    main(parser.parse_args())
