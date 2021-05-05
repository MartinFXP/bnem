# B-NEM

Boolean Nested Effects Models (B-NEM) are used to infer signalling
pathways. In different experiments (conditions) members of a pathway
(S-genes) are stimulated or inhibited, alone and in combination. In
each experiment transcriptional targets (E-genes) of the pathway react
differently and are higher or lower expressed depending on the
condition. From these differential expression profiles B-NEM infers Boolean
functions presented as hyper-edges of a hyper-graph connecting parents
and children in the pathway. For example if the signal is transducted
by two parents A and B to a child C and the signal can be blocked with
a knock-down of either one, they are connected by a typical
AND-gate. If the signal is still transduced during a single
knock-down, but blocked by the double knock-down of A and B, they
activate C by an OR-gate. In general the state of child C is defined
by a Boolean function

f: {0,1}^n -> {0,1}, C = f(A_1, ... , A_n)

with its parents A_i, i ∈ {1,...,n}.


Install:
--------

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("bnem")
```

Most recent (devel) version:

```r
install.packages("devtools")

library(devtools)

install_github("MartinFXP/bnem")

library(bnem)
```

Then check out the vignette for working examples.

```r
vignette("bnem")
```

Publication result:
-------------------

Use the function ```?processDataBCR``` to reproduce the data analysed in 
the publication (Pirkl et. al., 2016).
 
References:
-----------
 
Pirkl, Martin, Hand, Elisabeth, Kube, Dieter, & Spang,
Rainer. 2016. Analyzing synergistic and non-synergistic interactions
in signalling pathways using Boolean Nested Effect
Models. \textit{Bioinformatics}, 32(6), 893–900.

Pirkl, Martin. 2016. Indirect inference of synergistic and
alternative signalling of intracellular pathways. University of
Regensburg.
