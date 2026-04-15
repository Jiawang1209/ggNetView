# Build a graph object from GO enrichment analysis

Build a graph object from GO enrichment analysis

## Usage

``` r
build_graph_from_enrichGO(df)
```

## Arguments

- df:

  Data frame. The output from GO enrichmen analysis in clusterProfiler

## Value

An graph object with GO enrichment analysis

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires clusterProfiler enrichment result
ego <- clusterProfiler::enrichGO(gene, OrgDb = org.Hs.eg.db, ont = "ALL")
obj <- build_graph_from_enrichGO(df = as.data.frame(ego))
} # }
```
