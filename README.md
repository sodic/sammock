[![Build Status](https://travis-ci.org/sodic/sammock.svg?branch=master)](https://travis-ci.org/sodic/sammock)

# sammock
A simple program for generating [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) alignment files together with their corresponding reference [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files using a visual (and hopefully intuitive) alignment description format. The goal is to simplify the creation of arbitrary SAM files and make testing easier. 

## 1. Dependencies
- Python 3.6+

## 2. Installation
The program comes as a single _Python3_ module. You can either <a href="sammock.py" download>download only the required file</a>  or clone the entire repository:
```bash
git clone https://github.com/sodic/sammock.git
```
Position yourself in the same directory as the _sammock.py_ file and use _Python3.6+_ to run the program:
```bash
python3.6 sammock.py --help
```
## 3. Usage
The program reads data form an input text file describing alignments between multiple reads and a single reference sequence in a symbolic, human-readable way. It then produces a `FASTA` file containg the reference sequence and a (properly sorted) `SAM` file containing the alignment information. It can be used to quickly create aligment toy samples and simplify testing. Run `python3.6 sammock.py --help` for information about the command line parameters and their default values.


### 3.1. Input format specification
This is the specification of the symbolic alignment format which the program expects (the format is very simple so it might help to look at [the example](#32-example) prior to reading the specification):
- All empty rows are ignored (any row without a non-whitespace character is considered empty).
- Each non empty row represents one sequence.
- The file must contain at least two sequences (one read and a reference).
- The reference is always considered to be the last sequence in the file.
- Each sequence consists of symbols denoting the bases, allowed symbols are `A`, `C`, `G`, `T` and `-` which denotes a 'missing base'.
- Spaces between bases are ignored and can be used to enhance readability (e.g. `A C - C T` instead of `AC-CT`).
- An index of a base is determined by the number of non-whitespace characters preceding it.
- The reference sequence must not start with a 'missing base' symbol
- Preceding 'missing bases' on a read determine its starting position wrt. the reference (e.g. read `---ACG` starts at position `4` of the reference, one-indexed as in the SAM format)
- Trailing 'missing bases' and are removed from the sequence.

### 3.2. Example
Specify the input file, e.g. `sample.txt` with the following content:
```
A  A  A  A  -  G  C  C  T  T  A  C  T  A  A 
-  -  -  -  -  -  C  -  -  T  A  -  -  A  A  G
A  A  C  A  C  G  C  C  T  T  A  -  -  -  A  G  T
```
Run the program:
```bash
python3.6 sammock.py --alignments alignments_file.sam --reference ref_file.fa
```
This produces two files: 
- `ref_file.fa` containing the full reference:
```
>ref
AACACGCCTTAAGT
```
- and the `alignments_file.sam` containing the alignment data:
```
@HD	VN:1.4	SO:coordinate
@SQ	SN:ref	LN:14
1	0	ref	1	60	4M1D6M3I1M	*	0	0	AAAAGCCTTACTAA	@@@@@@@@@@@@@@
2	0	ref	7	60	1M2D2M1I2M	*	0	0	CTAAAG	@@@@@@
```

**NOTE:** Following from [the specification](#31-input-format-specification), the same ouput would also be produced with the following input files:
```
AAAA-GCCTTACTAA
------C--TA--AAG
AACACGCCTTA---AGT
```
or
```
 A   A   A   A   -   G   C   C   T   T   A   C   T   A   A   -   -  -
 
 -   -   -   -   -   -   C   -   -   T   A   -   -   A   A   G   -  -
 
 A   A   C   A   C   G   C   C   T   T   A   -   -   -   A   G   T  -
```
etc.

## 4. Additional
You can use the `prepareSample` bash script to fully automate the process of creating `FASTA` and `SAM` files and creating the indices (a `BAM` file):
```bash
./prepareSample <input_file_name> <sam_file_name> <ref_file_name>
```
