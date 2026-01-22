


## 1. Pipeline bio-informatique 


## 1.1. Base calling, démultiplexage et contrôle qualité	

Le sequençage des librairies a été réalisé selon le parametre OverrideCycles,U7N1Y143;I10;I10;U7N1Y143. En effet, la librairie contient 2 index de 10 nucleotides, le sequençage de l’insert est fait dans les 2 sens sur 151 cycles (151 nt) et les inserts sont encadrés par des UMI de 7 + 1 nt. (voir samplesheet en annexe)
Les fichiers bcl contiennent les images sortant du séquenceur illumina (données brutes). Le base Calling (appel des bases) est la première étape et consistait à inférer la séquence des inserts dans les deux sens (matrice et anti-sens) et à démultiplexer ces séquences entre tous les patients présents en 2 fichiers fastq contenant respectivement leur reads dans le sens 5’ à 3’ et 3’ à 5’ à l’aide de l’outil bcl-convert fournit par Illumina selon les instructions contenues dans notre samplesheet  :
```
nohup /usr/local/bin/bcl-convert -bcl-input-directory --output-directory --samplesheet

```
Une fois les fichiers FASTQ générés, l’outil Multiqc a interrogé tous les fichiers générés par le Dragen et résumé l’ensemble des criteres qualités dans un rapport lisible. 
Il faut s’assurer principalement que la qualité de séquençage est supérieure à 30 (phred score)et que les reads sont d’une taille entre 75 et 200 pb comme attendu.

## 1.2. Alignement et détection des variants 

Le pipeline Dragen (version4) détecte les variants suivant ce schéma 

<img width="1017" height="566" alt="image" src="https://github.com/user-attachments/assets/eee62cc7-00e6-4ff9-9d45-2b8e6b491f41" />


Les fichiers Fastq générés par le démultiplexage des short reads issus du séquenceur contiennent toutes les séquences du patients traités. Ces sequences ont été alignés avec le génome de référence afin de générer un fichier binaire d’extension «. Bam » comprenant les séquences alignées dans le sens allant de 5’ à 3’ et de 3’ à 5’ avec le génome de référence hg38. Ce fichier est appelé Bam. Il contient également des donnees de qqualité de l’alignement pour chaque read.
Les Bams ont été ensuite utilisés par le module Dragen de détection des variants en comparant les allèles avec les allèles de référence du génome hg38. Pour chaque variation statistiquement relevante avec la reference est associé un « allele alternatif » ou « variant ». L’ensemble de ces variants est compilé dans un fichier VCF. 
La qualité de l’alignement a été vérifié à travers le rapport Html fourni par le logiciel MultiQC 

## 1.3. Annotation des variants 


 Il existe de nombreuses bases de données publiques qui recensent des informations concernant les variants connu, notamment gnomAD. Cette derniere collectionne la frequence des variants observés dans la population generale (hors pathologie). l’annotation des variants a été faite avec l’outil ANNOVAR sur les fichiers VCF selon la référence exome hg38 de la base de données des variants « GnomAD_211 ».
 
## 1.4. Filtration des variants	




L’utilisation d’outil informatique de prediction basé sur un apprentissage necessite des donnees de bonne qualitée. Il a été essentiel de filtrer les artefacts techniques. les artefacts qui sont des faux positifs dus à des erreurs de séquençage ou de préparation des librairies mais classifies en tant que variants ont été filtrés suivant deux  étapes :

•	Filtrer les indels 

 Si au moins deux allèles alternatifs  sont localises sur la même position par rapport au génome de référence et qu’au moins l’un des allèles alternatif est un indel, alors tous les allèles alternatifs dans cette même  position ont été filtrés.
Cette étape a été réalisée en python avec le script « vcf_indel_filter.py » disponible depuis GitHub.

•	Filtrer les artéfacts avec Random Forest

<img width="919" height="550" alt="image" src="https://github.com/user-attachments/assets/2c00bbfe-5c02-408a-a620-3a9f0e2c1fec" />

<img width="1014" height="563" alt="Modele Random forest 2 " src="https://github.com/user-attachments/assets/ab78dcda-3f9c-44a5-90ef-73ed2ce615b4" />


Un modèle de machine Learning RandomForest non supervisé a été développé dans le but de classer les variants selon un score de confiance. Pour cela, le modèle doit apprendre sur un set de données considéré comme vrai. L’idée originale de Brand et al a été de se baser sur les variants fréquents dans la population générale. En d’autres termes, le modèle a été entrainé a reconnaitre les variants fréquents dans la population générale en utilisant des donnees intrinseques (genotype) et technique uniquement. L’outil ne sera pas capable de prédire la fréquence du variant dans la population générale mais reconnaitra les caractéristique des vraix variants.
(AF_popmax > 0.1) classés en  en labels positifs dans le but de pouvoir prédire avec précision la probabilité qu’un variant  soit  fréquent ou non dans la population générale  en attribuant un score binaire de 1 attribué aux variants ayant un( AF_popmax >  0.1 )  et de 0 pour le reste.
Des features d’entrainement ont été prises en compte comme : 

● SOR : Score du test de biais de brin mesurant  la présence de l’allèle de référence et alterne dans le brin du sens 5’ a 3’ et de et de 3’ a 5’.
● MQ : Qualité moyenne de l’alignement sur le site de l’allèle alterne 
● MQRankSum : Score du test de biais de la qualité moyenne de l’alignement entre l’allèle de référence et l’allèle alterne .
● ReadPosRankSum : Score du test de la position moyennes des bases des allèles alterne par rapport aux allèle référence sur le brin 
●Depth : Profondeur des reads = 
●Indel :  en leur attribuant un score binaire qui qualifie leur absence par 0 ou présence par 1 

Par la suite, 800 estimateurs ont été fixés comme paramètre. 75% des données ont été réparties dans le groupe d’apprentissage et 25% dans le groupe test.
Après entrainement, le modèle a calculé la probabilité d’appartenance de chaque variant aux variants fréquents, et par conséquent le label 1.
Ensuite, Ces variants ont été tries selon leur probabilité à travers un score de probabilité Cutoff déterminé à partir de la courbe de performance du modèle, la courbe de ROC pour une sensibilité de 0,99 .
Les variants ont été afin de classer selon leur probabilité d’appartenance aux variants fréquents en attribuant un label « PASS » aux variants fréquents qui a été donc supérieure au Cutoff et un label « FAIL »  pour le reste des variants .

Les données utiles sont extraites du fichier VCF à l’aide de bcftools en utilisant l’option 
query -f '%CHROM\t%POS\t%REF\t%ALT\t%SOR\t%MQ\t%MQRankSum\t%ReadPosRankSum\ t%DP’ 

Le script permettant la filtration selon le modèle « RandomForest2.py » est disponible depuis GitHub. 

## 1.5. Estimation de la fraction fœtale et génotypage des variants 

<img width="1027" height="498" alt="image" src="https://github.com/user-attachments/assets/18c11efb-7d08-42f2-8393-556b7d0b98ba" />

•	Filtrer les  VAF et la longueur des inserts

Les variants ont été filtres selon leurs VAF et leurs Z-scores (voir annexe 2) en admettant un intervalle de [0.025 0.975] pour les VAF et un intervalle [-4 4] pour le Z-score.

•	Génotypage et modélisation des de la fraction fœtale 

Le modèle de l’estimation de la fraction fœtale est un modèle de Deep Learning probabiliste qui va simultanément prédire la fraction fœtale et assigner un génotype pour les variants répertoriés dans notre séquence d’ADN circulant.
Le modèle a été défini sur 2 dimensions qui sont respectivement la VAF et le Z-score du FragmentSizeRankSumTest.
La fraction fœtale a été définie en tant que variable latente (f) qui a permis de définir les distributions moyenne de 5 clusters de 5 génotypes différents dont les moyennes de distribution respectives sont calculées comme suit : 

• (“cluster 0”: foetal 0/1, maternel 0/0): f / 2
• (“cluster 1: foetal 0/0, maternel 0/1”): (1 - f) / 2 
• (“cluster 2: foetal 0/1, maternel 0/1”): 0.5
• (“cluster 3: foetal 1/1, maternel 0/1”): f + (1 - f) / 2
• (“cluster 4: foetal 0/1, maternel 1/1”): 1 - (f / 2)



(0 : allèle de réference,1 : allèle alterne).
Le modèle a été établi à partir de la librairie Pyro de Pytorch. Chaque cluster a été modélisé de manière indépendante selon la variation stochastique avec le guide Auto_Delta de Pyro pour l’apprentissage sur les paramètres latents.  Ceci a permis d’assigner une probabilité de chaque variant d’appartenir aux 5 différents clusters décrits ci-dessus.
Par conséquent, les génotypes maternels et fœtaux ont été déduits à partir des probabilités p de chaque variant d’appartenir à un cluster spécifique :

Fœtal homozygote référence 0/0= p (cluster 1)
Fœtal heterozygote référence 0/1= p (cluster 0) + p (cluster 2) + p (cluster 4) 
Fœtal homozygote alterne 1/1= p (cluster 3)
Maternel homozygote référence 0/0= p (cluster 0)
Maternel hétérozygote 0/1= p (cluster 1) + p (cluster 2) + p (cluster 3)
Maternel homozygote alterne 1/1= p (cluster 3)
Le script du modèle « ML_fetal.py » est disponible depuis le lien GitHub
