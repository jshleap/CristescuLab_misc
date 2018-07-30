# Blast Command Line Tutorial

We will discuss basic usage of the blast+ tools (available in graham and cedar or at the [NCBI FTP](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)).

## Objectives / learning outcomes:
At the end of this tutorial you shoud be able to:
1. Download databases from NCBI 
2. Create your own custom local database
3. Do a basic/intermediate nucleotide blast search
4. Change the search parameters to fulfill your needs
5. Search for help regarding the available parameters

## Prerequisites:
Before we start, make sure you went over your bash notes, since we will be using several of the commands we saw in prevoious tutorials. Despite I will touch briefly on the web-based blast, I am assuming you know the basics on how to make blast searches on the [BLAST website](https://blast.ncbi.nlm.nih.gov/Blast.cgi).


I will discuss only the `blastn` (nucleotide blast) for time sake, but will comment on other types of blasts briefly.

## Outline of the tutorial
1. Introduction to BLAST: Types of databases, searches (this is brief, remember to brush up on this)
2. Downloading databases from NCBI FTP with `update_blastdb.pl`
3. Basic `blastn` search
4. Creating local databases with `makeblastdb`
5. Tunning parameters
6. Basic bash manipulation of output
