# FindRear :smile:
## A tool to construct somatic rearrangements :bowtie:

### Requirement
* python (>3.7)
* numpy
* pandas
* matplotlib
* argparse

### A. Input file 
#### Input_SV_file (6 columns about one SV are required)
1. SN: sample name
2. chr_1: the chromosome of low-end breakpoint
3. pos_1: genomic position of low-end breakpoint
4. flag_1: the orientation of low-end breakpoint
5. chr_2: the chromosome of low-end breakpoint
6. pos_2: genomic position of low-end breakpoint
7. flag_2: the orientation of low-end breakpoint

#### Input_bp/kb list (two columns are required)
1. sn: sample name that is consistent with sample name in input_SV_file
2. link: the bp/kb path 

#### bp/kb file (normalized coverage of windows (200-bp or 10kbp))
1. bp_file: the normalized coverage of non-overlapping windows (200bp)
2. kp_file: the normalized coverage of non-overlapping windows (10kbp)<br>

    *note: 

    a. normalzied coverage is obtained from PATCHWORK results;
    
    b. the shape of bp/kp matrix is N and M (N = the number of SV breakpoints; M = 201, 
    100 windows on the either left or right side of breakpoint) 

#### optional: --type
1. -type: choose the complex SV type you want to identify, 
consisting of complex TD/DEL, fold and unbalanced inversions.


### B. identify rearrangements (example)<br>
a. identify fold-back inversion<br>
```
python  ../construct_rearrangements.py --type fold-back test.sv test.out bp.list kp.list
```

b. identify unbalnced inversions<br>

```
python  ../construct_rearrangements.py --type inv test.sv test.out bp.list kp.list
```

c. identify TD-del rearrangements<br>
```
python ../construct_rearrangements.py --type td-del test.sv test.out bp.list kp.list
```

### C. compute the genomic path based on two anchored breakpoints (under development)

### D. visualize the SV with coverage (under construction)<br>
a. visualize simple SV<br>
b. visualize complex deletions/ TDs<br>
c. visualize fold-back rearrangements<br>
d. visualize chromothripsis<br>
