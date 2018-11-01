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
- Each read consists of symbols denoting the bases, allowed symbols are `A`, `C`, `G`, `T` and `-` which denotes a 'missing base'. The bases must be delimited by one or more spaces.
- Base qualities are specified after a colon following the base value, like this: `A:45 C:23 G:54 T:12`. One read may or may not contain base qualities. However, if it does, a quality must be specified for every non missing base of the read.
- The reference sequence must not start with a 'missing base' symbol.
- The reference sequence must not contain base qualities.
- Preceding 'missing bases' on a read determine its starting position wrt. the reference (e.g. read `- - - A C G` starts at position `4` of the reference, one-indexed as in the SAM format).
- Trailing 'missing bases' are truncated from both the reads and the reference.

### 3.2. Example
Specify the input file, e.g. `sample.txt` with the following content:
```
A:42  A:42  A:42  A:42  -     G:43  C:45  C:42  T:44  T:44  A:43  -     -     A:42  A:43
-     -     -     -     -     -     C:24  -     -     T:25  A:26  -     -     A:27  A:28  G:29
A     A     A     A     -     G     C     C     T     T     A     C     T     A     A

A     A     C     A     C     G     C     C     T     T     A     -     -     -     A     G      T
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
1	0	ref	1	60	4M1D6M1I1M	*	0	0	AAAAGCCTTAAA	KKKKLNKMMLKL
3	0	ref	1	60	4M1D6M3I1M	*	0	0	AAAAGCCTTACTAA	*
2	0	ref	7	60	1M2D2M1I2M	*	0	0	CTAAAG	9:;<=>
```

**NOTE:** Spaces can be added to increase readibility. Sequence bases don't need to be aligned for the program to work. Following from [the specification](#31-input-format-specification), the same ouput would also be produced with this (much less readable) input file:
```
A:42  A:42  A:42  A:42  -           G:43  C:45  C:42  T:44  T:44  A:43  -     -     A:42  A:43
- - - - - - C:24 - - T:25 A:26 - - A:27 A:28 G:29
A A A A - G C C T T A C T A A

A A C A C G C C T T A - - - A G T
```

## 4. Additional
You can use the `prepareSample` bash script to fully automate the process of creating `FASTA` and `SAM` files and creating the indices (a `BAM` file). However, this requires `samtools` to be installed/in scope:
```bash
./prepareSample <input_file_name> <sam_file_name> <ref_file_name>
```
