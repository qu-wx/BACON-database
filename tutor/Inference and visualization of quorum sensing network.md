

# Inference and visualization of quorum sensing network


---

# 1.Load BACON library

```r
library(BACON)
#>载入需要的程辑包：dplyr

#>载入程辑包：‘dplyr’

#>The following objects are masked from ‘package:stats’:

#>    filter, lag

#>The following objects are masked from ‘package:base’:

#>    intersect, setdiff, setequal, union

#>载入需要的程辑包：igraph

#>载入程辑包：‘igraph’

#>The following objects are masked from ‘package:dplyr’:

#>    as_data_frame, groups, union

#>The following objects are masked from ‘package:stats’:

#>    decompose, spectrum

#>The following object is masked from ‘package:base’:

#>    union

#>载入需要的程辑包：ggplot2
#>Warning messages:
#>1: 程辑包‘igraph’是用R版本4.2.3 来建造的 
#>2: 程辑包‘ggplot2’是用R版本4.2.3 来建造的 
```

# 2.Create BACON object

```r
## load sample data, here we use single-microbe RNA sequencing data of Escherichia coli as the example. please download "example_data.rda" at https://github.com/qu-wx/BACON-database/tree/main/tutor; please download E.coli database "Escherichia_coli_db.rda" at https://github.com/qu-wx/BACON-database/tree/main/Intra-species
load("example_data.rda")
load("Escherichia_coli_db.rda")
## creat BACON object; choose the database Escherichia_coli_db
## set column leiden in sample_metadata as group annotation
ecoli_sample <- createBACON(as.matrix(sample_matrix),group.by = sample_metadata$leiden,database = Escherichia_coli_db)
#> Create a BACON object from a data matrix
```
# 3.RUN BACON and visualization
```r
## run BACON 
ecoli_sample <- run_BACON(ecoli_sample,M=100)
#>Time difference of 0.620218 secs
##see strength after p-value filt
ecoli_sample@net
#>$LuxS_LsrB
#>            cluster_0 cluster_7
#>cluster_0 0.000000000         0
#>cluster_7 0.001243413         0
#>
#>$Tdh_QseC
#>          cluster_0 cluster_7
#>cluster_0         0         0
#>cluster_7         0         0
## circlr plot visual
netVisual_circle_BACON(ecoli_sample@net$LuxS_LsrB,vertex.label.cex = 1)
```
![sample_output](https://i-blog.csdnimg.cn/direct/a77a88a113784257914ed35be5c7dd5d.jpeg#pic_center)

```r
## chord diagram visual 
netVisual_chord_BACON(ecoli_sample,interaction_use = 'LuxS_LsrB')
```
![在这里插入图片描述](https://i-blog.csdnimg.cn/direct/369173b9e4f2487f91fb05764c31d59b.jpeg#pic_center)
