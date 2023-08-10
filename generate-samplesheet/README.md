# Usage
The function can accessed directly by calling `./generate-samplesheet.R`, with the following arguments:
```
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
  --investigator , -i   name of the investigator, default is: 'anonymous'
  --experiment , -e     experiment name, default is: 'UCSF experiment'
  --date , -d           date of the run, default is today: '08/10/2023'
  --read1 , -r          length of read 1, default is: '66'
  --read2 , -u          length of read 2, default is: '0'
  --i7r , -v            length of the first index read, i7
  --i5r , -x            length of the second index read, i5
```

# Output

The new samplesheet is automatically named in the format of `samplesheet_[experiment]_[date].csv`.

# Sample barcoding input file

The most important input is the sample barcoding file, which should formated according to the type of Lexogen barcode used.
## 12nt UDI

```
Sample_Name,Sample_Plate,Sample_Well
sample1,Set_A1,A1
sample2,Set_A3,A2
```
*n.b.* not all 12 nucleotides of the barcodes need to be used. The length of the barcodes to be sequenced can be modified using the `i7r` and `i5r` flags. 

## 6nt barcodes

```
Sample_Name,Sample_Well_i7,Sample_Well_i5
sample1,A1,B1
sample2,A1,B2
```

*n.b.* if only using the i7 barcode and not the i5 (or vice versa) barcodes can be excluded from the final Sample Sheet by setting either the `i5r` or `i7r` flags to 0.
