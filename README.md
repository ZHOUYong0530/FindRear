# findRear :smile:
## a tool to construct somatic rearrangements :bowtie:

### 1. data-preprocessing
a. extract coverage of breakpoints<br>

### 2. identify rearrangements<br>
a. identify fold-back inversion<br>

python  ../construct_rearrangements.py --type fold-back test.sv test.out bp.list kp.list
  
b. identify unbalnced inversions<br>

python  ../construct_rearrangements.py --type inv test.sv test.out bp.list kp.list
	
c. identify TD-del rearrangements<br>

python ../construct_rearrangements.py --type td-del test.sv test.out bp.list kp.list

### 3. visualize the coverage<br>
a. visualize simple SV

b. visualize complex deletions/ TDs

c. visualize fold-back rearrangements

d. visualize chromothripsis
