#!/usr/bin/env python
#
# collapsed reads .txt
#


import argparse

def make_and_print_table(datas):
    """
        Generated a contamination table

        columns are reference genomes, and rows are samples.
    """
    sample_names = {}
    
    header = ['sample_id']
    with open(datas[0]) as data_file:
        for line in data_file:
            header.append(line.split()[0])

    for data_f in datas:
        sample_name = data_f.split('.')[0]
        sample_names[sample_name] = []
        with open(data_f) as data_file:
            for line in data_file:
                float_f = line.split()[1]
                sample_names[sample_name].append(float_f)
    print('\t'.join(header))
    for sample, line in sample_names.items():
        print(sample+'\t'+'\t'.join(line))    


def main():
    """
        This program processes the final.merged files.
    """
    parser = argparse.ArgumentParser(description="Contamination table")
    parser.add_argument('data', nargs='*', help='contamination data')
    args = parser.parse_args()
    assert len(args.data) >=1, "Program requires at least one data file to run" 
    make_and_print_table(args.data)

if __name__ == "__main__":
    main()
