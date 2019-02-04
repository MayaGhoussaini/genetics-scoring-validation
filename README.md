Open Target Genetics Validation
===============================

Scripts to validate Open Targets Genetics scores.

### Input requirements

Columns

```
set (str, required): set name
source (str, required): source of the association (e.g. pumed ID)
rsid (str, required if varid null)
varid (str, required if rsid null): chrom_pos_ref_alt (GRCh37)
gene_id (str, required): Ensembl gene ID of validated target (GRCh37)
confidence (High/Medium/Low): confidence of the association
study_id (str, nullable): we need to be able to trace the publication to a study ID in Open Targets in order to validate L2G scoring (not necessary for V2G)
```

### Requirements

- gql

### Todo

- Normalise within groups (`varid`)
