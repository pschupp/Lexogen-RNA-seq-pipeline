# Usage
The function can accessed directly by calling `./generate-samplesheet.R`, with the following arguments:
```
usage: ./generate-samplesheet.R [-h] [--rc] [--investigator] [--experiment]
                                [--date] [--read1] [--read2] [--i7r] [--i5r]
                                samplebarcs

Generate Sample Sheets for sequencing of Lexogen FWD 3' Libraries. Supports
both basic i5/i7 6nt indexing and the UDIs.

positional arguments:
  samplebarcs      location of the csv file matching sample names to barcode
                   plates and wells

optional arguments:
  -h, --help       show this help message and exit
  --rc             reverse complement the i5 barcode, default is 'FALSE'
  --investigator   name of the investigator, default is: 'anonymous'
  --experiment     experiment name, default is: 'UCSF experiment'
  --date           date of the run, default is today: 'MM/DD/YYYY'
  --read1          length of read 1, default is: '66'
  --read2          length of read 2, default is: '0'
  --i7r            length of the first index read, i7, default is: '6'
  --i5r            length of the second index read, i5, default is: '6'
```

# Output

The new Sample Sheet is named in the format of `samplesheet_[experiment]_[date].csv`.

# Sample barcoding input file

The only required input is the sample barcoding file, which should be formated according to the type of Lexogen barcode used.

## 12 nt UDI

```
Sample_Name,Sample_Plate,Sample_Well
sample1,Set_A1,A1
sample2,Set_A3,A2
```
*n.b.* not all 12 nucleotides of the barcodes need to be used. The length of the barcodes to be sequenced can be modified using the `i7r` and `i5r` flags. Lexogen advises shortening to either 8 or 10 nt.

## 6 nt barcodes

```
Sample_Name,Sample_Well_i7,Sample_Well_i5
sample1,A1,B1
sample2,A1,B2
```

*n.b.* if only using the i7 barcode and not the i5 barcode (or vice versa), barcodes can be excluded from the final Sample Sheet by setting either the `i5r` or `i7r` flags to 0.
