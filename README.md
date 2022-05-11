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
- Barcodes file should a tab-delimited table formatted as below:
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
  -b, --barcodes           Barcodes table file name (tab-delimited) (string)
  -f, --fuzzy-threshold    Fuzzy index match threshold (int [=1])
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
$ ./demultiplex_satay -r test/test_reads.fastq -i test/test_index.fastq -b test/barcodes.txt
```
Expected matching rates:
```
demultiplex_satay v0.0.1

Provided sequencing reads file name:           test/test_reads.fastq
Provided index reads file name:                test/test_index.fastq
Provided barcode file name:                    test/barcodes.txt
Provided fuzzy mapping threshold:              1

--------------------------------------------------------------

Demultiplexing reads...

Reading barcodes file...
Read in 18 barcodes

Reading index file...
Read 8 total index reads
Matched 5 index reads
Sample index match rate: 62.5%

Reading sequence read file...
Read 8 total sequencing reads
Matched 5 sequencing reads
Sequencing read match rate: 62.5%


Processing complete.

Elapsed time:                                  0.01562s
```
7 reads should map to `total_library_DpnII_701`, 1 to `high_Pi_top_10_20_DpnII_710`, and 3 to `unassigned`

```
$ ./demultiplex_satay -r test/test_1000_R1.fastq -i test/test_1000_R2.fastq -b test/barcodes.txt
```

Note:    
Built off of [`fastp_lite`](https://github.com/XPRESSyourself/XPRESSpipe/tree/main/fastp_lite).