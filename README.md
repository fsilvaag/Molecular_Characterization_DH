# Molecular Characterization of DH Lines of Maize

Esta carpeta contiene la siguiente informacion:

## Carpetas:
01. Filtrado --> Contiene todos los arhivos intermedios resultado del filtrado de la base d datos.

                 El filtrado corresponde a los siguientes:
                 
                             - Remover SNPs duplicados (142 SNPs en total)
                             
                             - Remover INDEL (0 SNPs)
                             
                             - Remover SNPs con mas de 2 alelos (0 SNPs)
                             
                             - Remover SNPs cpn > 5% informacion perdida (22589 SNPs)
                             
                             - Imputacion usando Beagle 
                             
                             - Remover Monomorphic markers (2204 SNPs)
                             
                 Igualmente se incluye el resultado del analisis de LD usando TASSEL ("LD_Decay.txt")
                 

02. IBDLD    --> Contiene los resultados de correr el programa IBDLD

                             - Archivos terminados en "ibdtxt.gz"          = Probabilidad de IBD de cada SNP
                             
                             - Archivos terminados en ".kinship"           = Provee el kinship de cada par de individuos Par = [(502 x 501)/2 + 502]
                             
                             - Archivos terminados en ".gtype" y ".mthd"   = Archivo de uso interno por IBDLD program
                             
                             - Archivos terminados en ".hbdtxt"            = Probabilidad de "Homozygous by descent". Esto seria la probabilidad de hbd de cada SNP
                             
                             - Archivos terminados en ".segment"           = Los IBD segments identificados despues de usar una longitud minima de 500 Kb, minimo 50 SNPs, y que cada SNP tenga una probabilidad de IBD >= 0.7
                             
                             - Archivos terminados en ".segment.transosed.txt" = el mismo .segment pero transpuesto
                             
                             - Archivos terminados en ".primal"            = Es el "Chromosome-wide Identity Coefficient" que puede ser usado para imputacion y phasing.
                             
                             - Archivos terminados en ".reg"               = Archivo de uso interno por IBDLD program
                             
                             
03. Resultados  --> Contiene 2 archivos:

                             - Figure_Structure.xlsx = La tabla de contribucion genetica de cada parental a las 487 DH lines. La contribucion fue calculada como
                                                          100*(# IBD SNPs entre progenitor y DH Line)/(7993 SNPs usados)
                                                          
                             - ld_decay_table.csv    = Resultados del analisis LD Decay ("LD_Decay.R")
                             



## Archivos

1. transposed.awk                  --> codigo para trasponer un archivo en Unix

2. Dendrogram.R                    --> Codigo en R para hacer dendrogramas (En este caso es muy dificil usar el circular dendrogram por la cantidad de genotipos. tocaria buscar otra forma de representarlos. Tal vez eliminar el nombre de las DH lines y colocarlos como codigo o algo asi)

3. LD_Decay.R                      --> Codigo en R para estimar el LD decay. El arcivo ha sido modificado del original dado por Leandro. Las modificaciones fueron asi: Leanddro - Arthur - Fernando

4. Remove_Monomorphic_markers.Rmd  --> Codigo en R modificado de Gustavo de los campos en Cimmyt. El codigo fue dado por Edna Mageto y modificado por Fernando Silva

5. Structure_Graph_Data.Rmd        --> Codigo en R para producir la grafica de contribucion de cada progenitor a los DH lines
                                                 Esta forma de visuaizacion debe ser mejor pensada ya que hay muchas DH lineas por cada ciclo. Podria mostrarse solo algunas por cada ciclo?
                                                 
6. Analisis.txt                    --> Archivo explicando el paso a paso del filtrado y de como correr IBDLD

7. No_Duplicates.hmp.txt           --> Base de datos genotipica despues de remover los 142 SNPs duplicados

8. Raw_geno_data_502.hmp.txt       --> Base de datos genotipica original

9. IBDLDv3.38_Manual.pdf           --> Manual del programa IBDLD

10. README.txt                     --> Este archivo
