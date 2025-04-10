# ECCFP - EccDNA Caller based on Consecutive Full Pass
A software package for identifying eccDNAs in the ONT sequencing data of RCA-amplified eccDNA
## Introduction
ECCFP utilizes all consecutive full passes to determine accurate eccDNA positions and generate consensus sequences, based on rolling circle amplification and nanopore sequencing.

#### Reference
ECCFP is optimized from the [Flec](https://github.com/icebert/eccDNA_RCA_nanopore.git) software.
- [eccDNA_RCA_nanopore](https://github.com/icebert/eccDNA_RCA_nanopore.git)
#### Dependency (python packages)
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [pyfaidx](https://pypi.org/project/pyfaidx/)
- [pyfastx](https://pypi.org/project/pyfastx/)
- [Biopython](https://biopython.org)
## Install
```
git clone https://github.com/WSG-Lab/ECCFP.git
cd ECCFP
pip install .
```
## Usage
#### simple usage
```
eccfp --fastq input.fastq --paf mapping.paf --reference ref.fa -o output
```
#### Parameter
```
usage: eccfp.py [options]

options:
  -h, --help            show this help message and exit

  Generate candidate eccDNA from rolling circle reads:

  --fastq FASTQ         input reads in fastq format
  --paf PAF             input alignments in PAF format
  --reference REFERENCE
                        reference genome sequences in fasta format
  --output OUTPUT, -o OUTPUT
                        output folder, the default is the current folder
  --maxOffset MAXOFFSET
                        maximum offset of start/end positions between two sub-reads to be considered as mapping to the same location, default is 20
  --minMapQual MINMAPQUAL
                        minimum mapping quality of sub-reads, default is 30

  Get accurate eccDNA location from candidate eccDNA:

  --fluctuate FLUCTUATE
                        maximum offset between the candidate eccDNA and the candidate eccDNA in a group, default is 20bp
  --nf NF               minimum number of fullpass to support candidate eccDNA, default is 2
  --cov COV             minimum fullpass used for merging a set of candidate eccDNA, default is 1
  --nc NC               minimum number of candidate eccDNA required to merge overlapped candidate eccDNAs, default is 2

  Generate consensus sequences and variants:

  --minDP MINDP         minimum depth to call variants, default is 4
  --minAF MINAF         minimum alternative allele frequency to call variants, default is 0.75

```
## Output
Five files are generated following the completion of the pipeline: unit.txt and candidate_consolidated.txt(intermediate files), final_eccDNA.txt, consensus_sequence.fasta, and variant.txt (result files).   
|file|details|
|----|-------|  
|unit.txt|full pass alignment details for all candidate eccDNAs detected in reads|
|candidate_consolidated.txt|the consolidating steps used to derive accurate eccDNAs from the candidate eccDNAs|
|final_eccDNA.txt|accurate eccDNA information|
|consensus_sequence.fasta|the consensus sequences of accurate eccDNAs|
|variants.txt|variant profiles specific to these accurate eccDNAs|

###### Example final_eccDNA.txt file
|eccDNApos|Nfullpass|Nfragments|Nreads|refLength|seqLength|
|---------|---------|----------|------|---------|---------|
|chr16:757924-758274(+)|4|1|1|351|350|
|chr1:21828265-21828439(+)\|chr7:132167689-132167880(+)|24|2|2|367|367|

||description|
|---------|---------------|
|eccDNApos|eccDNA position|
|Nfullpass|Number of consecutive full pass for this eccDNA covered by all reads|
|Nfragments|Number of fragment that form this eccDNA|
|Nreads|Number of reads identified for the eccDNA|
|refLength|The length of reference genome that this eccDNA |
|seqLength|The length of consensus sequence that this eccDNA |

###### Example variants.txt file
|col1|col2|col3|col4|col5|col6|col7|col8|
|----|-------|-|-|--|--|------------|-----------------------|
|chr2|4898178|TGTG|-|11|13|InDel|chr2:4897293-4898597(-)|
|chr3|4745441|G|A|23|23|transition|chr3:4744817-4746914(+)|

||description|
|-------|-------|
|col1|chromsome|
|col2|position in the reference genome|
|col3|reference base|
|col4|variant|
|col5|supportive coverage depth|
|col6|total coverage depth|
|col7|type|
|col8|eccDNApos|

###### Example unit.txt file
|reads|Nfullpass|fragments|read_start|read_end|chr|start|end|strand|candidate_eccDNA|cigar|
|-|-|-|-|-|-|-|-|-|-|-|
|07ead43b-dc39-4ba2-a38a-bd8310e16b45|2|1|1117|2444|chr1|22649376|22650447|-|chr1:22649376-22650451(-)\|chr21:15109952-15110325(-)|51M1D284M3D465M272I129M2D46M2D38M8D43M|
|07ead43b-dc39-4ba2-a38a-bd8310e16b45|2|1|2821|3848|chr1|22649376|22650442|-|chr1:22649376-22650451(-)\|chr21:15109952-15110325(-)|343M2D422M28D186M2D38M8D34M1I4M|
|07ead43b-dc39-4ba2-a38a-bd8310e16b45|2|1|4215|4785|chr1|22649951|22650451|-|chr1:22649376-22650451(-)\|chr21:15109952-15110325(-)|220M80I186M2D38M8D47M|
|07ead43b-dc39-4ba2-a38a-bd8310e16b45|2|2|2445|2820|chr21|15109952|15110325|-|chr1:22649376-22650451(-)\|chr21:15109952-15110325(-)|4M1D114M1I87M1D4M1I148M2I15M|
|07ead43b-dc39-4ba2-a38a-bd8310e16b45|2|2|3849|4214|chr21|15109962|15110325|-|chr1:22649376-22650451(-)\|chr21:15109952-15110325(-)|109M2I87M1D4M1I163M|

###### Example candidate_consolidated.txt
|eccDNApos|eccDNA_len|consolidated|cand_eccDNA|cand_len|cand_Nfullpass|
|-|-|-|-|-|-|
|chr1:4496303-4496663(-)|361|False|chr1:4496303-4496663(-)|361|5|
|chr1:36762100-36762408(-)|309|True|chr1:36762083-36762406(-)|324|2|
|chr1:36762100-36762408(-)|309|True|chr1:36762083-36762408(-)|326|2|
|chr1:36762100-36762408(-)|309|True|chr1:36762083-36762409(-)|327|2|
|chr1:36762100-36762408(-)|309|True|chr1:36762084-36762408(-)|325|4|
|chr1:36762100-36762408(-)|309|True|chr1:36762084-36762414(-)|331|2|
|chr1:36762100-36762408(-)|309|True|chr1:36762097-36762403(-)|307|2|
|chr1:36762100-36762408(-)|309|True|chr1:36762100-36762408(-)|309|2|
|chr10:24967157-24968326(+)\|chr12:26660679-26660872(+)|1364|True|chr10:24967149-24968326(+)\|chr12:26660679-26660872(+)|1372|2|
|chr10:24967157-24968326(+)\|chr12:26660679-26660872(+)|1364|True|chr10:24967154-24968326(+)\|chr12:26660679-26660872(+)|1367|2|
|chr10:24967157-24968326(+)\|chr12:26660679-26660872(+)|1364|True|chr10:24967157-24968325(+)\|chr12:26660679-26660870(+)|1361|2|
|chr10:24967157-24968326(+)\|chr12:26660679-26660872(+)|1364|True|chr10:24967157-24968326(+)\|chr12:26660679-26660872(+)|1364|2|
|chr10:24967157-24968326(+)\|chr12:26660679-26660872(+)|1364|True|chr10:24967158-24968326(+)\|chr12:26660679-26660872(+)|1363|2|
