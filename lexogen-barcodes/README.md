# Lexogen Barcoding Schemes

|Name|SKU|Max multiplex|Single/dual|BC length|URL|
|---|---|---|---|---|---|
|Lexogen UDI v2|Dual|12nt||
|Lexogen i5 6nt dual indexing|047.4x96|96|Single|6|https://www.lexogen.com/wp-content/uploads/2023/01/047UG109V0300_Lexogen-i5-6-nt-Dual-Indexing-Add-on-Kits-5001-5096_2023-01-03.pdf|
|Lexogen i7 6nt dual index|044.96(included with QuantSeq)|96|Single|6|https://www.lexogen.com/wp-content/uploads/2023/01/047UG109V0300_Lexogen-i5-6-nt-Dual-Indexing-Add-on-Kits-5001-5096_2023-01-03.pdf|
|Lexogon i7/i5 6nt dual index|combo of above|96x96|dual|6|https://www.lexogen.com/wp-content/uploads/2021/04/107UI264V0104_UDI-12-nt-Index-Sequences-for-Illumina_2021-03-26.xlsx|


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
1. Read 1
2. index read preparation
3. Index 1 Read (i7)
4. Index 2 Read (i5)
5. plus Read 2 Resynthesis
6. and Read 2 for paired-end flow cells / runs (except for sequencing above where i5 read comes after resynthesis)
- Stylized as: R1 x i7 x i5 x R2