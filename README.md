# ahpC

Code used to generate table summary tables of ahpC promoter mutations from a tb-profiler variant dump csv file and a csv file containing the isoniazid dst.

### Input file formats

#### Variant dump file:

| sample_id  | genome_pos | gene | change      | freq | type                  | sublin           | drtype    | drugs     |
|------------|------------|------|-------------|------|-----------------------|------------------|-----------|-----------|
| SRR5073738 | 2154724    | katG | p.Arg463Leu | 1    | missense_variant      | lineage2.2.1     | Sensitive |           |
| ERR181815  | 2155168    | katG | p.Ser315Thr | 0.7  | missense_variant      | lineage3.1.1     | HR-TB     | isoniazid |
| ERR2514579 | 2726105    | ahpC | c.-88G>A    | 1    | upstream_gene_variant | lineage3         | Sensitive |           |
| ERR3276067 | 2726051    | ahpC | c.-142G>A   | 1    | upstream_gene_variant | lineage1.2.1.2.1 | RR-TB     |           |

#### DST file:


| wgs_id     | dst |
|------------|-----|
| SRR5073738 | 0   |
| ERR181815  | 1   |
| ERR2514579 | NA  |