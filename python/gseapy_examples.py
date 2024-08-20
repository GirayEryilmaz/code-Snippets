import gseapy as gp

# Get avaliable libraries:
gp.get_library_name()

# Get the gene lists of a library (returns dict)
gp.get_library('GO_Biological_Process_2023')

# Gene set enrichment analysis with a gene list
enr = gp.enrichr(gene_list=['ISG15', 'MX1', 'STAT2'],
                 gene_sets=['MSigDB_Hallmark_2020','KEGG_2021_Human', 'GO_Biological_Process_2023', 'WikiPathway_2023_Human'],
                 organism='human',
                 outdir=None, # don't write to disk
                )
enr.results.sort_values(by = 'Combined Score', ascending=False)
