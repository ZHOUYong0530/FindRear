# consRear3
a tool to construct rearrangements


1. data-preprocessing
a. extract coverage of breakpoints<br>

2. identify rearrangements<br>
  a. identify fold-back inversion<br>
  python  ../construct_rearrangements.py --type fold-back test.sv test.out bp.list kp.list
  
  b. identify unbalnced inversions<br>
  python  ../construct_rearrangements.py --type inv test.sv test.out bp.list kp.list
