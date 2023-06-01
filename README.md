# Towards Finding S-box Circuits with Optimal Multiplicative Complexity

This repository includes a tool for finding MC-optimal S-box circuits.  
If you need to generate circuits quickly, you can choose an option that enables finding a circuit that is not MC-optimal (e.g., -b and -c).  
Please refer to https://github.com/usnistgov/Circuits/tree/master/data/mc_dim for all files in **mc_dim** directory.
For more details, please refer to the paper.

## Options:
-S <value> (or <file>) : S-box (or S-boxes) to implement  
-N <value> : the name of the S-box  
-n <value> : the size of the s-box  
-a : Use affine equivalence option (Section 4.1)  
-b <value> : Use the optional paprameter kappa (Section 4.2)  
-c <value> : Use the optional paprameter kappa without searching with distance of 1 (Section 4.2)  
-M <value> : Number of threads to use for multithreading
  
## Usage: 
Run the command below.
>python Astar.py -S 1,a,4,c,6,f,3,9,2,d,b,7,5,0,8,e -N GIFT -n 4  
>python Astar.py -S 0x1,0xa,0x4,0xc,0x6,0xf,0x3,0x9,0x2,0xd,0xb,0x7,0x5,0x0,0x8,0xe -N GIFT -n 4  
>python Astar.py -S 1,10,4,12,6,15,3,9,2,13,11,7,5,0,8,14 -N GIFT -n 4  
>python Astar.py -S 5bitSboxes.txt -n 5 -a -b 4 -M 8  
>python Astar.py -S 6bitSboxes.txt -n 5 -a -c 4 -M 48  
  
Please note that the format inside {name}.txt should follow the format of the files inside the **example** directory.
