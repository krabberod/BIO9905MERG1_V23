# Bioinformatics for Environmental Sequencing (DNA metabarcoding)
**These pages contain course material for BIO9905MERG1 Spring 2023**  
For the course page on the University of Oslo [webpage, click here](https://www.uio.no/studier/emner/matnat/ibv/BIO9905MERG1/).

## Course content
For mapping and exploring communities of both micro- and macroorganisms, high throughput sequencing (HTS) of environmental DNA has become a powerful approach. One can either analyze the total DNA content to obtain knowledge about which genes are present (DNA metagenomics) or sequence a selected PCR-amplified marker (DNA metabarcoding) to obtain information about the taxonomic composition. We will focus on the latter approach in this course. The students will be introduced to important analytical bioinformatics approaches from the processing of raw sequence data to establishment of the OTU/sample matrix and retrieval of taxonomic identity of the sequences.

Important themes will be (1) filtering and quality assessment of high throughput sequence data, (2) error correction and/or clustering of high throughput sequence data, and (3) taxonomic annotation of high throughput sequence data. We will also touch upon some further downstream analyses, including network analyses and evolutionary placement of HTS onto backbone phylogenies. Applications of a wide suite of tools will be presented, including VSEARCH and DADA2.

The course will be a blend of presentations, guest lectures, discussion and a few hands-on sessions. All hands-on secession will be run in R on your local laptop/computer. Hence, all participants should have R and selected R packages installed – see information below.

## Schedule

The course will run from 17-21 April, 9:00-17:00 (times may vary). For a detailed overview of the program, see below.

## Report
Those of you that attend the course through the research schools or UiO and want to obtain ECTS credits, will have to hand in a report before June 1th.
For the report, you should write a 5-page text (minimum) about a fictive research project where you will use DNA-metabarcoding to explore the community composition and diversity of a certain habitat and/or ecological gradient. You are free to select the organismal group(s) and the habitat/gradient. In the text you should: (1) define the goal(s) of the study, (2) describe the sampling design, (3) the wet-lab work (briefly) and, most important, (4) the bioinformatics analyses. Not only describe how you plan to carry out the research, but also why you make your choices. On point (4), describe in detail how you plan to analyze your data and which bioinformatics approaches you will use and why. Also mentioned what you expect to be problematic and which steps that might introduce bias(es) to your results, as well as what type of bias(es). Concerning the format, you should use Times New Roman size 12, 1.5 line spacing and 2.5 cm margins.

The report should be sent to: haavarka@ibv.uio.no

## Teachers
The main teachers will be Ramiro Logares, Anders K. Krabberød, Micah Dunthorn, Torbjørn Rognes, Frederic Mahé and Håvard Kauserud (organizer), but other experts will provide guest lectures (see table).



## Program

| Day           | Time (start) | Topic                                                                                                       | Responsible                   |
| ------------- | ------------ | ----------------------------------------------------------------------------------------------------------- | ----------------------------- |
| **Monday**    | 09:00        | [Introduction to DNA metabarcoding](Lectures_and_groups/Intro_lecture_Kauserud.pdf)                         | Håvard Kauserud               |
|               | 10:00        | *Introduction continued*                                                                                    |                               |
|               | 11:00        | [Introduction to sequencing techniques](Lectures_and_groups/20210503_Lyle.pdf)                              | Håvard Kauserud               |
|               | 12:00        | _Lunch break_                                                                                               |                               |
|               | 13:00        | Introduction to [Linux](intro.to.unix) and [R](intro.to.r)/[Rstudio](intro.to.Rstudio)                      | Ramiro Logares                |
|               | 14:00        | Get to know each other (group work)                                                                         | Håvard Kauserud               |
|               | 15:00        | Cutadapt and sequence cleaning                                                                              | Torbjørn Rognes               |
|               | 16:00        | [Help with setup and installation of required packages](Setup)                                              | Anders K. Krabberød           |
|               | 17:00        | *PIZZA*                                                                                                     |                               |
|               |              |                                                                                                             |                               |
| **Tuesday**   | 09:00        | [Introduction to DADA2](Dada2_Pipeline)                                                                     | A. Krabberød / R. Logares     |
|               | 10:00        | *DADA2 continued*                                                                                           |                               |
|               | 11:00        | *DADA2 continued*                                                                                           |                               |
|               | 12:00        | _Lunch break_                                                                                               |                               |
|               | 13:00        | *DADA2 continued*                                                                                           |                               |
|               | 14:00        | *DADA2 continued*                                                                                           |                               |
|               | 15:00        | About DADA2 and its further applications                                                                    | Ben Callahan (Zoom)           |
|               | 16:00        | *About DADA2 continued*                                                                                     |                               |
|               |              |                                                                                                             |                               |
| **Wednesday** | 09:00        | Introduction to VSEARCH and Swarm                                                                           | Torbjørn Rognes               |
|               | 10:00        | from FASTQ files to OTU tables                                                                              | Frédéric Mahé                 |
|               | 11:00        | *continued*                                                                                                 |                               |
|               | 12:00        | _Lunch break_                                                                                               |                               |
|               | 13:00        | LULU and MUMU (C++ version of LULU)                                                                         | Frédéric Mahé                 |
|               | 14:00        | [Multivariate analyses of DNA-metabarcoding data](community.ecology)                                        | Douglas Yu (Zoom)             |
|               | 15:00        | [Group discussions](Lectures_and_groups/Group_work_Wednesday.pdf)                                           |                               |
|               |              |                                                                                                             |                               |
| **Thursday**  | 09:00        | Introduction to long-read DNA metabarcoding                                                                 | Mahawash Jamy                 |
|               | 10:00        | Phylogenetic placement/binning of HTS data                                                                  | Micah Dunthorn/Lucas Czech    |
|               | 11:00        | *continued*                                                                                                 |                               |
|               | 12:00        | _Lunch break_                                                                                               |                               |
|               | 13:00        | OTUs, ASVs and phylospecies                                                                                 | Micah Dunthorn                |
|               | 14:00        | Contamination, library prep, aerial metabarcoding                                                           | Kristine Bohmann              |
|               | 15:00        | *continued*                                                                                                 |                               |
|               |              |                                                                                                             |                               |
| **Friday**    | 09:00        | Metacoder                                                                                                   | Ella Thoen                    |
|               | 10:00        | *continued*                                                                                                 |                               |
|               | 11:00        | Downstream analyses                                                                                         |                               |
|               | 12:00        | _Lunch break_                                                                                               |                               |
|               | 13:00        | [Methods to retrieve intra-species diversity information](Lectures_and_groups/Wangensteen_Intraspecies.pdf) | Owen S. Wangensteen Fuentes   |
|               | 14:00        | Case studie on DNA metabarcoding                                                                            | Sundy Maurice (guest lecture) |
|               | 15:00        | Group work, Q&A                                                                                             |                               |


---
# Software
We will use R (version 4.0.5 or later) and Rstudio (version 1.4.1 or later) in this course. In addition, we will use Google Colab for programs that require a Linux/Unix environment.

**Everybody should download and install R (https://www.r-project.org/), Rstudio (https://www.rstudio.com/) and the required packages before the course starts**.

For more information about the required packages [Click](Setup/) here](Setup/).

---


# Suggested reading (reviews)
You can find the PDFs [here:](Suggested_reading/)
- Zinger et al. 2019. DNA metabarcoding—Need for robust experimental designs to draw sound ecological conclusions. Molecular Ecology, 28, 1857-1862
- Deiner et al. 2017. Environmental DNA metabarcoding: Transforming how we survey animal and plant communities. Molecular Ecology, 26, 5872-5895.
- Bohmann et al. 2014. Environmental DNA for wildlife biology and biodiversity monitoring. TREE, 29.
- Alberdi et al. 2017. Scrutinizing key steps for reliable metabarcoding of environmental samples. Methods in Ecology and Evolution, 9, 134-147.
- Ficetola et al. 2016. How to limit false positives in environmental DNA and metabarcoding? Molecular Ecology, 16, 604-607.
- Dickie et al. 2018. Towards robust and repeatable sampling methods ineDNA-based studies. Molecular Ecology Resources, 18, 940-952.
- Schnell et al. 2015. Tag jumps illuminated – reducing sequence‐to‐sample misidentifications in metabarcoding studies. Molecular Ecology Resources, 15, 1289-1303.Antich et al. 2021. To denoise or to cluster, that is not the question: optimizing pipelines for COI metabarcoding and metaphylogeography. BMC Bioinformatics, 22,177.
- Lamb et al. 2019. How quantitative is metabarcoding: A meta‐analytical approach. Molecular Ecology, 28, 420-430.
 -Jamy et al. 2019. Long-read metabarcoding of the eukaryotic rDNA operon to phylogenetically and taxonomically resolve environmental diversity. Molecular Ecology Resources, 20, 429-443.
 ----
### Supported by [Digitalt Liv Norge](https://www.digitallifenorway.org/), [ForBio](https://www.forbio.uio.no/), and [Norbis](https://norbis.w.uib.no/)
![](images/2021/04/Artboard2x.png)  
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fkrabberod%2FBIO9905MERG1_V23&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)  
All the keywords in this explanation, by the way, are totally misleading, due to the everyday quirks of language. **Don DeLillo, Ratner's Star**.

---
