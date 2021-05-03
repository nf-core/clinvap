#!/usr/bin/env python3
import sys
import subprocess


def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    # get the rows passing the filter
    filter_process = subprocess.run(['bcftools', 'filter', '--include', 'FILTER="PASS" || FILTER="."',
                                    input_file, '--output', output_file, '--output-type', 'v'], stderr=subprocess.PIPE)
    warnings = filter_process.stderr.decode('utf-8').split('\n')
    for warn in range(len(warnings)):
        print(warnings[warn])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} <input_vcf_file> <output_file>".format(sys.argv[0]))
        exit(1)
    main()
