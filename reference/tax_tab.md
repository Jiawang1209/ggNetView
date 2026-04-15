# Taxonomy annotation table

Taxonomic classification (Kingdom to Species) for the OTUs in
\[otu_tab\]. The first column \`OTUID\` is used as the join key.

## Usage

``` r
tax_tab
```

## Format

A \`tbl_df\` with 2,859 rows and 8 columns (\`OTUID\`, \`Kingdom\`,
\`Phylum\`, \`Class\`, \`Order\`, \`Family\`, \`Genus\`, \`Species\`).

## Source

Internal example data shipped with ggNetView.

## Examples

``` r
data(tax_tab)
head(tax_tab)
#> # A tibble: 6 × 8
#>   OTUID  Kingdom  Phylum          Class              Order  Family Genus Species
#>   <chr>  <chr>    <chr>           <chr>              <chr>  <chr>  <chr> <chr>  
#> 1 ASV_2  Archaea  Thaumarchaeota  Unassigned         Nitro… Nitro… Nitr… Unassi…
#> 2 ASV_3  Bacteria Verrucomicrobia Spartobacteria     Unass… Unass… Spar… Unassi…
#> 3 ASV_31 Bacteria Actinobacteria  Actinobacteria     Actin… Mycob… Myco… Unassi…
#> 4 ASV_27 Archaea  Thaumarchaeota  Unassigned         Nitro… Nitro… Nitr… Unassi…
#> 5 ASV_9  Bacteria Unassigned      Unassigned         Unass… Unass… Unas… Unassi…
#> 6 ASV_30 Bacteria Acidobacteria   Acidobacteria_Gp16 Unass… Unass… Gp16  Unassi…
```
