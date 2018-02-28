#!/usr/bin/bash

grep -f <(tail -n 6 ~/sexual_pop/filtering_vcf_from_1102/pass_fail_summary.txt) ~/sexual_pop/sequencekeyfiles/CB7GTANXX_8_barcode.txt | \
cut -f3,4 > fail_barcodes.txt

grep -f <(cut fail_barcodes.txt -f1) ~/sexual_pop/old_tassel_runs/tassel_run_0928/barcode_counts.txt > fail_counts.txt

join fail_barcodes.txt fail_counts.txt -1 1 -2 2

