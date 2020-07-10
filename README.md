This is repository contains the scripts to reproduce all figures and analysis from the publication:
<b>Time-Restricted Feeding Alters Lipid and Amino Acid Metabolite Rhythmicity without Perturbing Clock Gene Expression</b> published in
<b>Nature Communications</b>

Periodicly rhythmic features are detected using the rhythmic.detection.R, and the figures and integrated analysis is performed in the integrated.analysis.R script. Both depend on the functions found in rhythmic.functions.R

SessionInfo for all analysis

R version 3.6.3 (2020-02-29)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 17763)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
 [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] org.Hs.eg.db_3.10.0    AnnotationDbi_1.48.0   IRanges_2.20.2         S4Vectors_0.24.4       Biobase_2.46.0        
 [6] BiocGenerics_0.32.0    openxlsx_4.1.5         ggupset_0.3.0          waffle_1.0.1           Rtsne_0.15            
[11] impute_1.60.0          VennDiagram_1.6.20     futile.logger_1.4.3    UpSetR_1.4.0           scales_1.1.1          
[16] patchwork_1.0.1        clusterProfiler_3.14.3 data.table_1.12.8      ggplot2_3.3.2          reshape2_1.4.4        

loaded via a namespace (and not attached):
 [1] enrichplot_1.6.1     bit64_0.9-7          RColorBrewer_1.1-2   progress_1.2.2       httr_1.4.1          
 [6] tools_3.6.3          DT_0.14              R6_2.4.1             DBI_1.1.0            colorspace_1.4-1    
[11] withr_2.2.0          tidyselect_1.1.0     gridExtra_2.3        prettyunits_1.1.1    extrafontdb_1.0     
[16] curl_4.3             bit_1.1-15.2         compiler_3.6.3       formatR_1.7          xml2_1.3.2          
[21] labeling_0.3         triebeard_0.3.0      ggridges_0.5.2       stringr_1.4.0        digest_0.6.25       
[26] DOSE_3.12.0          htmltools_0.5.0      pkgconfig_2.0.3      extrafont_0.17       htmlwidgets_1.5.1   
[31] rlang_0.4.6          rstudioapi_0.11      RSQLite_2.2.0        gridGraphics_0.5-0   farver_2.0.3        
[36] generics_0.0.2       jsonlite_1.6.1       gtools_3.8.2         BiocParallel_1.20.1  zip_2.0.4           
[41] GOSemSim_2.12.1      dplyr_1.0.0          magrittr_1.5         ggplotify_0.0.5      GO.db_3.10.0        
[46] Matrix_1.2-18        Rcpp_1.0.4.6         munsell_0.5.0        viridis_0.5.1        lifecycle_0.2.0     
[51] stringi_1.4.6        ggraph_2.0.3         MASS_7.3-51.6        plyr_1.8.6           qvalue_2.18.0       
[56] blob_1.2.1           gdata_2.18.0         ggrepel_0.8.2        DO.db_2.9            crayon_1.3.4        
[61] lattice_0.20-41      graphlayouts_0.7.0   cowplot_1.0.0        splines_3.6.3        hms_0.5.3           
[66] pillar_1.4.4         fgsea_1.12.0         igraph_1.2.5         futile.options_1.0.1 fastmatch_1.1-0     
[71] glue_1.4.1           packrat_0.5.0        lambda.r_1.2.4       BiocManager_1.30.10  vctrs_0.3.0         
[76] tweenr_1.0.1         urltools_1.7.3       Rttf2pt1_1.3.8       gtable_0.3.0         purrr_0.3.4         
[81] polyclip_1.10-0      tidyr_1.1.0          ggforce_0.3.2        europepmc_0.4        tidygraph_1.2.0     
[86] viridisLite_0.3.0    tibble_3.0.1         rvcheck_0.1.8        memoise_1.1.0        ellipsis_0.3.1 
