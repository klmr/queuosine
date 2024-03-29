#+ echo=FALSE
knitr::opts_chunk$set(dev = c('png', 'pdf'))

#+ echo=TRUE
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
    top_n(1, nchar(Sequence)) %>%
    slice(1) %>% # Break ties
    ungroup()

codon_usage = cds_$cu$cu(coding_sequences)
relative_codon_usage = codon_usage %>%
    group_by(Gene) %>%
    mutate(CU = CU / sum(CU)) %>%
    ungroup()

codon_pattern = '^.A[CU]$'
filtered_codon_usage = relative_codon_usage %>%
    filter(grepl(codon_pattern, Codon)) %>%
    group_by(Gene) %>%
    summarize(SumUsage = sum(CU)) %>%
    mutate(SSE = (SumUsage - mean(SumUsage)) ** 2,
           Score = SSE / max(SSE)) %>%
    arrange(desc(Score))

#+ codon-usage
modules::import_package('ggplot2', attach = TRUE)
theme_set(theme_bw() + theme(panel.border = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             axis.line = element_line(colour = "black")))
ggrepel = modules::import_package('ggrepel')

ggplot(filtered_codon_usage) +
    aes(x = '', SumUsage) +
    geom_boxplot() +
    labs(x = '', y = 'Usage of codons •A[CU]') +
    ggrepel$geom_text_repel(aes(label = Gene),
                            data = head(filtered_codon_usage, 14),
                            nudge_x = 0.1) +
        theme(axis.ticks.x = element_blank())

#+ genelist
head(filtered_codon_usage, 14)

io$write_table(select(filtered_codon_usage, Gene, Score), 'data/genelist.csv', col.names = FALSE)

piano = modules::import_package('piano')

load_gene_set = function () {
    path = 'data/gene_association.wb.gz'

    read = function () {
        on.exit(close(f))
        f = gzfile(path, 'r')
        suppressWarnings(io$read_table(f, sep = '\t', comment.char = '!')) %>%
            select(Gene = 3, GO = 5) %>%
            piano$loadGSC(type = 'data.frame')
    }

    go = try(read(), silent = TRUE)

    if (! inherits(go, 'try-error'))
        return(go)

    download.file('http://geneontology.org/gene-associations/gene_association.wb.gz',
                  path)
    read()
}

go_genes = load_gene_set()

untidy = function (tidy_data, rownames = 1)
    `rownames<-`(as.data.frame(tidy_data[-rownames]), tidy_data[[rownames]])

gsa = function (.data, go_genes, genes = Gene, stat = Stat)
    gsa_(.data, go_genes, substitute(genes), substitute(stat))

gsa_ = function (.data, go_genes, genes = 'Gene', stat = 'Stat') {
    stopifnot(inherits(go_genes, 'GSC'))
    data = .data %>%
        # In case `stat` contains a calculation.
        mutate_(` stat ` = stat) %>%
        select_(genes, '` stat `') %>%
        untidy()
    piano$runGSA(data, gsc = go_genes, verbose = FALSE)
}

go_label = function (terms) {
    get_label = function (term) {
        json = modules::import_package('jsonlite')
        uri = 'http://www.ebi.ac.uk/ols/api/ontologies/go/terms?short_form=%s'
        result = json$fromJSON(sprintf(uri, sub(':', '_', term)))
        result$`_embedded`$terms$label[[1]]
    }

    vapply(terms, get_label, character(1))
}

gsa_result = gsa(filtered_codon_usage, go_genes, Gene, 1 - Score)
enriched_go_terms = piano$GSAsummaryTable(gsa_result) %>%
    select(Name, padj = `p adj (non-dir.)`) %>%
    filter(padj < 0.99) %>%
    mutate(Label = go_label(Name)) %>%
    arrange(padj)

#+ enriched-go-terms
enriched_go_terms
io$write_table(enriched_go_terms, 'data/enriched-go-terms.txt')
