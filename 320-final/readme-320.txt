To find a set of strings of length 320 with rate of stars 8, use the following code.

C-code: ced-320.cpp

To compile: g++ ced-320.cpp -o ced-320 -O3

To run: ced-320

Generates initial set of strings satisfying property 3a.

To generate more strings faster in parallel, you can run "ced-320-paral.cpp" once you have some initial set (100+) of strings and certain initial thresholds are calculated.
To do so, from the output of ced-320, save the three lines starting by "int Pthr[]=", "int Sthr[]=" and "int Mthr[]=" into a file "320-thr.h". 
This contains calculated thresholds and it is necessary to compile "ced-320-paral.cpp".
Also in the output, locate the last set of lines starting by "col[   0]=\"*0010110*1010000..." and ending by e.g "col[ 224]=\"*0000101*1111110*...",
copy them into a file "ced-extend-parallel-320-8-0001.out". 
From each line in "ced-extend-parallel-320-8-0001.out" remove the initial "col[ nnn]=\"" and final "\";", so each line should contain just the string, i.e. it should look like "*0111100*1001100*1001...10*1011110\n". The length of the file should be divisible by 321=320+1.
You are ready to proceed with ced-320-paral.cpp.

Relevant parameters in the C-code:

#define LEN  320        // string length
#define RATE 8          // rate of stars in the strings
#define COLSIZE 800     // maximal collection size
#define RANKFRAC 0.0001 // ced-320 first samples about 10^5 pairs of strings and measures their pairwise edit distance, and the distance of their various fragments. Based on that it prepares some set of thresholds. The collection is then built by trying new random strings. New strings that are too close to strings already in the collection will be rejected. So only strings whose distance and the distance of their fragments is above the thresholds will be retained for the collection. RANKFRAC is used to determine the thresholds - it specifies what fraction of randomly chosen pairs of strings violates the threshold. So the lower the number, the fewer random strings will violate it but the worse the guarantee on the final edit distance. The target edit distance for the chosen thresholds is reported on the 3rd line of the output from ced-320 as the "min". So one would like to set RANKFRAC as high as possible to maximize min. But then random strings will be rejected at higher rate. To build the collection in reasomable time one needs to set RANKFRAC to a suitable number. RANKFRAC in [0.0001,0.001] seems to work fairly well.


----

C-code: ced-320-paral.cpp

To compile: g++ ced-320-paral.cpp -o ced-320-paral -O3

To run: ced-320-paral [<seed>]

Parallel version of ced-320.
It looks for additional strings, and adds them to the file "ced-extend-parallel-320-8-0001.out".
Multiple instances of ced-320-paral can be run in parallel accessing the same file "ced-extend-parallel-320-8-0001.out".
From time to time each instance re-reads the current file "ced-extend-parallel-320-8-0001.out" and tries to extend it by a next possible string.
Each instance should use a different <seed> for its pseudoraodnom generator.
Each seed is a non-negative integer.

Due to the parallelism, several instances might contribute by incompatible strings at the very same time.
Hence, the final set of strings in "ced-extend-parallel-320-8-0001.out" might not satisfy condition 3a.
Running "ted-320 3" will re-test the property 3a on the set and suggest which strings to remove to fix the property.
See further ted-320.cpp

Relevant parameters in the C-code:

#define LEN  320        // string length
#define RATE 8          // rate of stars in the strings
#define COLSIZE 800     // maximal collection size
#define ALPHA  0.1625   // This is the distance parameter alpha. It should be set to min/LEN, where "min" is described in the comment for RANKFRAC of ced-320. 
It is used to cap all the thresholds at LEN*ALPHA.

----

C-code: ted-320.cpp  

To compile: g++ ted-320.cpp -o ted-320 -O3

Tests properties 2-4 for a candidate collection of strings specified in the array "colstr[]" in the C-code.

Before compilation one ought to include the tested collection of strings from "ced-extend-parallel-320-8-0001.out" into the C-code.
Variable "const char *colstr[]={}" should contain the list of the strings in quotation marks and comma-separated.

To run: ted-320 [<property> [<range>]]

where <property> could be:

  2 ... Test property 2.
  3 ... Test property 3a.
  4 ... Test property 4.
  5 ... Test property 3b.
  6 ... Test a simplified property 2 quickly.

<range> is an integer i, where the test will be performed only for strings in the range [20*i, 20*(i+1)] for properties 3 and 4, and in the range [10*i, 10*(i+1)] for property 2. This allows to parallelize the test. 
For example, running "ted-320 4" on a collection of 67 strings is equivalent to running "ted-320 4 0 & ted-320 4 1 & ted-320 4 2 & ted-320 4 3"

Property 3a can be tested quite quickly.
If it is violated it suggests which strings to remove from the collection.
You can remove all of them or some, and re-test the property. (See useful tricks to remove all of them.)

Property 2 is the slowest to test. A quick pre-test can be done by running "ted-320 6". This tests the property on pairs of strings instead of triples.
It identifies violating strings so that the full test can be run only once on the "final" candidate set.

Relevant parameters in the C-code:

#define LEN  320        // string length
#define RATE 8          // rate of stars in the strings
#define COLSIZE 800     // maximal collection size
#define ALPHA  0.1625   // This is the distance parameter alpha for which all properties are tested. Distances that are exactly at the threshold ALPHA*LEN are accepted. (So we use non-strict inequality to pass the strings.)

----

Useful tricks

1) To remove strings with double monochromatic blocks and starting by a monochromatic block (which will eventually violate Property 2) run:

  egrep -v "(^.(0000000|1111111)|(0000000|1111111).(0000000|1111111))"  ced-extend-parallel-320-8-0001.out >cleaned.out
  mv cleaned.out ced-extend-parallel-320-8-0001.out

(One can do it from time to time even while instances of ced-320-paral are running.)

2) To remove strings that violate property 3a after running "ted-320 3 >ted-320-violated-sample.out" you can use:

  egrep "^[ 0-9]* .[01]{7}" ted-320-violated-sample.out|sort|uniq -c|sort -n

This will list all strings from the tested collection. Strings to remove are on lines starting by "   2 ". Strings to keep are on lines starting by "   1 ".

----

Shell scripts to run the full set of tests on 32-core machines:

ted-kamenik3      
ted-kamenik4      
ted-kamenozrout2
ted-lomikamen2

-----

Output files:

ced-320.out: output of ced-320

ced-320-paral-<i>.out: output of "i"-th instance of ced-320-paral

ted-320-0-p3a.out: Protocol of testing property 3a on the collection.

ted-320-<i>-p<prop>.out: Protocol of testing property "prop" for range "i".

