# consRear3
a tool to construct rearrangements


1. data-preprocessing
a. extract coverage of breakpoints

2. identify rearrangements
  a. identify fold-back inversion
  python  ../construct_rearrangements.py --type fold-back test.sv test.out bp.list kp.list
  
  b. identify unbalnced inversions
  python  ../construct_rearrangements.py --type inv test.sv test.out bp.list kp.list
