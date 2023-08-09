# Lexogen Barcoding Schemes

|Name|SKU|Max multiplex|Single/dual|BC length|URL|
|---|---|---|---|---|---|
|Lexogen UDI v1 <sup>\*</sup>|191-203;A1-A4,B1|384 with HD=2|Dual|12, 10, or 8nt|https://www.lexogen.com/wp-content/uploads/2021/04/107UI264V0104_UDI-12-nt-Index-Sequences-for-Illumina_2021-03-26.xlsx|
|Lexogen UDI v2 <sup>\+</sup>|191-203;A1-A4,B1|384 with HD=2|Dual|12, 10, or 8nt|https://www.lexogen.com/wp-content/uploads/2021/04/107UI264V0104_UDI-12-nt-Index-Sequences-for-Illumina_2021-03-26.xlsx|
|Lexogen i5 6nt dual indexing|047.4x96|96 with HD=0|Single|6|https://www.lexogen.com/wp-content/uploads/2023/01/047UG109V0300_Lexogen-i5-6-nt-Dual-Indexing-Add-on-Kits-5001-5096_2023-01-03.pdf|
|Lexogen i7 6nt dual index|044.96(included with QuantSeq)|96 with HD=0|Single|6|https://www.lexogen.com/wp-content/uploads/2023/01/047UG109V0300_Lexogen-i5-6-nt-Dual-Indexing-Add-on-Kits-5001-5096_2023-01-03.pdf|
|Lexogon i7/i5 6nt dual index|combo of above|9216 with HD=0|dual|6|https://www.lexogen.com/wp-content/uploads/2021/04/107UI264V0104_UDI-12-nt-Index-Sequences-for-Illumina_2021-03-26.xlsx|

<sup>\*</sup> Set B1 corresponds to Set A1, but with i5 sequences in the reverse complement direction (to accommodate for different Illumina chemistries); see not below on sequencing order.\
<sup>\+</sup> The Add-on UDI V2 kits include an improved amplification module that ensures better resistance to overcycling effects.\
HD refers to Hamming distance, or the number of errors that can be corrected after sequencing to recover more "Undetermined Reads".


# Important notes

## RC Barcodes

Barcodes must be reverse-complemented (RC) if using the following instrument/kit combination because the `i5 Read` is read **after** `Read 2 Resynthesis`:
- iSeq 100™ (PE flow cell)
- MiniSeq™ (all (PE) flow cells)
- NextSeq® 500/550/2000 (all (PE) flow cells)
- HiSeq® 3000/4000 (PE flow cells only)
- NovaSeq™ 6000 with v1.5 reagent kits (all (PE) flow cells)

## Sequencing order

The order of sequencing is always: 
```
1. Read 1
2. index read preparation
3. Index 1 Read (i7)
4. Index 2 Read (i5)
5. Read 2 Resynthesis
6. Read 2 for paired-end flow cells / runs (except for sequencing above where i5 read comes after resynthesis)
```
- Stylized as: R1 x i7 x i5 x R2
