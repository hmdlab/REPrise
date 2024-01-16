import sys
from argparse import ArgumentParser
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from tabulate import tabulate

def read_fai(filename):
    with open(filename, 'r') as f:
        total_size = 0
        seqs = []
        
        line = f.readline()
        while line:
            s = line.split()
            name = s[0]
            size = int(s[1])
            seqs.append(name)
            total_size += size
            line = f.readline()

    return (total_size, seqs)

def read_bed(filename, seq_names=None):#20220212_武田追記
    with open(filename, 'r') as f:
        seqs = defaultdict(lambda: [])

        line = f.readline()
        while line:
            s = line.split()
            name = s[0]
            begin = int(s[1])
            end = int(s[2])
            if seq_names is not None and not name in seq_names:
                continue
            seqs[name].append((begin, end))
            line = f.readline()

    return seqs


def read_repeatmasker_out(filename, seq_names=None):
    with open(filename, 'r') as f:
        for _ in range(3):
            f.readline()

        seqs = defaultdict(lambda: [])

        line = f.readline()
        while line:
            s = line.split()
            name = s[4]
            begin = int(s[5]) - 1
            end = int(s[6]) - 1
            if seq_names is not None and not name in seq_names:
                continue
            seqs[name].append((begin, end))
            line = f.readline()

    return seqs

def read_gff(filename, seq_names=None):
    with open(filename, 'r') as f:
        seqs = defaultdict(lambda: [])

        line = f.readline()
        while line:
            if line.startswith('#'):
                line = f.readline()
                continue
            s = line.split()
            name = s[0]
            begin = int(s[3]) - 1
            end = int(s[4]) - 1
            if seq_names is not None and not name in seq_names:
                continue
            seqs[name].append((begin, end))
            line = f.readline()

    return seqs

def read_phraider_elements(filename, seq_names=None):
    with open(filename, 'r') as f:
        seqs = defaultdict(lambda: [])

        line = f.readline()
        while line:
            if line.startswith('#'):
                line = f.readline()
                continue
            s = line.split()
            name = s[3]
            begin = int(s[4])
            end = int(s[5])
            if seq_names is not None and not name in seq_names:
                continue
            seqs[name].append((begin, end))
            line = f.readline()

    return seqs

def read_annot(filename, seq_names=None):
    if filename.endswith('.gff'):
        return read_gff(filename, seq_names)
    elif filename.endswith('elements'):
        return read_phraider_elements(filename, seq_names)
    elif filename.endswith('.out'):
        return read_repeatmasker_out(filename, seq_names)
    elif filename.endswith('.bed'):
        return read_bed(filename, seq_names)
    else:
        raise Exception('unknown format')

def overlap(x, y):
    if y[1] - 1 < x[0] or x[1] - 1 < y[0]:
        return 0
    else:
        return min(x[1], y[1]) - max(x[0], y[0])

def evaluate(seq_size, seq_names, gt, pred):
    truth = 0
    positive = 0
    true_positive = 0
    false_positive = 0

    for name in seq_names:
        gt_inttree = IntervalTree(Interval(x[0], x[1] + 1) for x in gt[name])
        gt_inttree.merge_overlaps()
        truth += sum(i.length() for i in gt_inttree)

        pred_inttree = IntervalTree(Interval(x[0], x[1] + 1) for x in pred[name])
        pred_inttree.merge_overlaps()
        p = sum(i.length() for i in pred_inttree)
        positive += p
        false_positive += p

        for pred_interval in pred_inttree:
            tp = sum(overlap(o, pred_interval) for o in gt_inttree.overlap(pred_interval))
            true_positive += tp
            false_positive -= tp

    untruth = seq_size - truth
    true_negative = untruth - false_positive

    assert(false_positive >= 0)
    assert(untruth >= 0)
    assert(true_negative >= 0)

    sens = true_positive / truth
    spec = true_negative / untruth
    positive_rate = positive / seq_size
    pnr = false_positive / seq_size
    prec = true_positive / positive
    f_score = 2 * sens * prec / (sens + prec)

    return (sens, spec, positive_rate, pnr, prec, f_score)

def main(args):
    seq_size, seq_names = read_fai(args.fai)
    gt = read_annot(args.gt, seq_names)
    pred = read_annot(args.pred)

    sens, spec, positive_rate, pnr, prec, f_score = evaluate(seq_size, seq_names, gt, pred)

    tab = [
        ['Sensitivity', sens * 100],
        ['Specificity', spec * 100],
        ['Positive rate', positive_rate * 100],
        ['Potential Novel Repeat', pnr * 100],
        ['Precision', prec * 100],
        ['F-score', f_score * 100],
    ]
    print('index:\t{}'.format(args.fai))
    print('ground truth:\t{}'.format(args.gt))
    print('prediction:\t{}'.format(args.pred))
    print(tabulate(tab, headers=['', '%bp']))

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('fai')
    parser.add_argument('gt')
    parser.add_argument('pred')

    main(parser.parse_args())
