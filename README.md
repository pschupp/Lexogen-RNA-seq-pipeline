# Introduction
This repository deals with ["QuantSeq 3’ mRNA-Seq Library Prep Kit FWD"](https://www.lexogen.com/quantseq-3mrna-sequencing/).

After generating your libraries you can proceed to:

1. Generating the Sample Sheet in the `generate-samplesheet` folder. For advice on barcode type selection see the folder `lexogen-barcodes` which lists the two barcode types used by Lexogen.
2. After sequencing run the main workflow (`main-workflow` folder) which goes from raw bcl files to expression matrices.
3. (Optional) For saturation curve analysis, see the folder `saturation-curve-workflow`.
