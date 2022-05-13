# stack-to-base
 A tool to search ShortStack sRNA annotations against miRBase
 
### ...A warning
This tool is in development, it is likely that many of the options and outputs may change. There are many bugs still being worked out.

### Why and how does it work
[ShortStack](https://github.com/MikeAxtell/ShortStack) is a powerful tool for sRNA-seq analysis. It produces an annotation, including high-precision predictions of miRNAs. However, precision comes with a trade-off of sensitivity and loci which contain miRNAs may not always be identified as such. It can also be useful to find family names associated with a sRNA locus. This tool is meant to bridge that gap and help identify loci that might be overlooked otherwise.

This tool uses 3 approaches to search for miRBase miRNA loci all with varying degrees of sensitivity.  
• Exact mature miRNA matches between a submitted miRNAs and MajorRNAs from ShortStack (the most precise).
• Genomic overlap of BLAT alignment of submitted hairpins vs ShortStack annotated loci.  
• BLAST alignment of ShortStack annotated loci vs submitted hairpins.  

Using these 3 approaches, it can find what family members are most similar to a locus and give a strong inference what miRNA families might be contained in a locus.  

### More warnings...
* MiRBase is not perfect. As of now this tool includes low-confidence miRNAs. These annotations should be considered miRNAs at the user's peril - many of them are likely false.  
* This will likely be only specific to the family of a miRNA. Any species information or paralog information should be disregarded - they mean nothing in the context of this tool. So, a hit that said Ath-miR166b when you are working in A. lyrata, is not evidence of horizontal gene transfer... You should only trust this as far as this is a miR166-family miRNA.  

### Installation dependencies
All of these software must be executable and in the $PATH variable:
* samtools  
* bedtools  
* blastn  
* blat  
* (python3)  

Past these, there is no installation. Simply run the stack-to-base.py script. To run this executable from anywhere, add the repo to your $PATH variable in your bash_rc or bash_profile, and confirm the tool is executable:

    chmod +x ./stack-to.base.py


### Input data

This tool uses hairpin and mature sequences from mirbase (version?). Both of these files have been supplied in the github repo and stack-to-base.py should recognize their location if the script is in the repo directory. Otherwise, you can specify their location manually.


### Syntax

Check the help document for syntax information: 

    stack-to-base.py -h
    
The tool simply takes the Results.txt output from a ShortStack run and the genome used to perform the annotation. An output name is supplied to handle intermediate files. 

    stack-to-base.py -r example/Results.txt -g ~/+Genomes/Plant_genomes/Solyc.chromosomes.4.00.fa -o test_output
    
A possibly better way to use this tool is to include non-mirbase annotations as well. These can be submitted as more arguments for the -m (Mature) and -p (hairPin) options.

    stack-to-base.py -r example/Results.txt -g ~/+Genomes/Plant_genomes/Solyc.chromosomes.4.00.fa -o test_output \
        -m miRBase21.mature.fa Lunardon2020.mature.fa \
        -p miRBase21.hairpin.fa Lunardon2020.hairpin.fa









