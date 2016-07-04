modules::import_package('dplyr', attach = TRUE)
io = modules::import('ebi-predocs/ebits/io')
cds_ = modules::import('./cds')

load_cds = function () {
    path = 'data/cds.rds'
    cds = try(readRDS(path), silent = TRUE)

    if (inherits(cds, 'try-error')) {
        bm = modules::import_package('biomaRt')
        ensembl = bm$useMart('ensembl', 'celegans_gene_ensembl')
        attributes = c('ensembl_gene_id', 'external_gene_name', 'coding')
        cds = bm$getBM(attributes = attributes, mart = ensembl, bmHeader = TRUE) %>% tbl_df()
        saveRDS(cds, path)
    }

    cds
}

cds = load_cds()

coding_sequences = cds %>%
    select(`Ensembl Gene ID`,
           Gene = `Associated Gene Name`,
           Sequence = `Coding sequence`) %>%
    filter(cds_$valid_cds(Sequence)) %>%
    group_by(Gene) %>%
    top_n(nchar(Sequence)) %>%
    slice(1)
    ungroup()

codon_usage = cds_$cu$cu(coding_sequences)
relative_codon_usage = cds_$cu$norm(codon_usage)

codon_pattern = '^.A[CU]$'
filtered_codon_usage = relative_codon_usage %>%
    filter(grepl(codon_pattern, Codon)) %>%
    group_by(Gene) %>%
    summarize(SumUsage = sum(CU)) %>%
    mutate(SSE = (SumUsage - mean(SumUsage)) ** 2,
           Score = SSE / max(SSE)) %>%
    arrange(desc(Score))

io$write_table(select(filtered_codon_usage, Gene, Score), 'data/genelist.csv', col.names = FALSE)
