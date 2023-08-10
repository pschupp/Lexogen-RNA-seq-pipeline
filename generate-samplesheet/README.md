# Usage


Sample_Name,Sample_Plate,Sample_Well
sample1,Set_A1,F4
sample2,Set_A1,F5
Sample_Name,Sample_Well_i7,Sample_Well_i5
sample1,F4,A1
sample2,F5,A1

Attaching package: ‘renv’

The following objects are masked from ‘package:base’:

    autoload, load, remove

Bioconductor version '3.14' is out-of-date; the current release version '3.17' is available with R version '4.3'; see
  https://bioconductor.org/install

Attaching package: ‘BiocManager’

The following object is masked from ‘package:renv’:

    install

[1] 2
Using BLAS/LAPACK: ATLAS
[0;38;5;40mnormal[0m [38;2;203;75;22mx[x<=-0.01][0m [38;2;131;148;150mx[abs(x)<0.01][0m [38;2;181;137;0mx[x>=0.01][0m [38;2;211;54;130m19/01/2038 03:14:07[0m [38;2;133;153;0m"string"[0m
[38;2;108;113;196mNA/NaN/NULL[0m [38;2;220;50;47mFALSE[0m [38;2;133;153;0mTRUE[0m [38;2;211;54;130mInf[0m [38;2;131;148;150m[index][0m [38;2;211;54;130mstderror[0m [38;2;181;137;0mwarn[0m [38;2;211;54;130merror[0m
usage: ./generate-samplesheet.R [-h] [--rc] [--samplebarcs] [--workingdir]
                                [--investigator] [--experiment] [--date]
                                [--read1] [--read2] [--i7r] [--i5r]

Generate Sample Sheets for sequencing of Lexogen FWD 3' Libraries. Supports
both basic i5/i7 6nt indexing and the UDIs.

optional arguments:
  -h, --help            show this help message and exit
  --rc , -c             reverse complement the i5 barcode
  --samplebarcs , -s    location of csv containing columns: sample name,
                        sample barcode plate, sample barcode well
  --workingdir , -w     current working directory
  --investigator , -i   name of the investigator, default is: 'anonymous'
  --experiment , -e     experiment name, default is: 'UCSF experiment'
  --date , -d           date of the run, default is today: '08/10/2023'
  --read1 , -r          length of read 1, default is: '66'
  --read2 , -u          length of read 2, default is: '0'
  --i7r , -v            length of the first index read, i7
  --i5r , -x            length of the second index read, i5
