# demultiplex_satay

Overview:    
Used for demultiplexing pooled Satay screen reads.

Usage:    
- Files should NOT be gzipped
- Files should end in the suffix, `.fastq`
- Single-end FASTQ file for reads should be 4 lines each, like below:
```
@NB501960:698:HMTN7BGXL:1:11101:20226:1065 1:N:0:1
CACATAGTTGGCAAATTGATCCTTGATCATTAAAATCATTNGAGAATCGTCTNCCANNNNNANNNNNNNNGTNNTN
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE#EEEEEEEEEEE#EEE#####E########EE##E#
```
- FASTQ file with read indices should be 4 lines each, like below:
```
@NB501960:698:HMTN7BGXL:1:11101:20226:1065 2:N:0:1
CGTACTAG
+
AAA6AE/E
```
- Barcodes file should be formatted as below:
```
TCGCCTTA	sample_1
CTAGTACG	sample_2
TTCTGCCT	sample_3
GCTCAGGA	sample_4
```

```
Usage: 

./demultiplex_satay --reads=<filename> --index=<filename> --barcodes=<filename> [options] ...


Parameters:

  -r, --reads              Input fastq file name (string)
  -i, --index              Index fastq file name (string)
  -b, --barcodes           File name for barcodes (.tsv or .txt) (string)
  -f, --fuzzy-threshold    Length of UMI, default is 1 (int [=1])
  -?, --help               print this message

```

Build on Linux:
```
$ make -f Makefile_Linux
```

Build on MacOS:
```
$ make -f Makefile_macOS
```

Test:
```
$ ./demultiplex_satay -r test/test_1000_R1.fastq -i test/test_1000_R2.fastq -b test/barcodes.txt
```

Note:
Built off of [`fastp_lite`](https://github.com/XPRESSyourself/XPRESSpipe/tree/main/fastp_lite).