## Test cases for core-accessory classification.

In this case, Test1 is the reference genome, it's genome file is Test1.genome.tab:

`$ cat Test1.genome_file.tab`

    Test1_contig1   15000  
    Test1_contig2   20000  
    Test1_contig3   10000  
    Test1_contig4   5000  
    Test1_contig5   22000  

The folder `coords_files` contains two files of coordinates of alignments of genome assemblies Test2 and Test3 to the assembly of Test1. E.g.:

`$ cat coords_files/Test1.vs.Test2.coords`  

    /path/to/fastafile/genome_Test1.fa /path/to/fastafile/genome_Test2.fa
    NUCMER

    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
    =====================================================================================
    2001     3100  |        0        0  |     1100        0  |    89.47  | Test1_contig1    Test2_Contig1
    3001     9000  |        0        0  |     6000        0  |    97.69  | Test1_contig1    Test2_Contig1
    8526    10000  |        0        0  |     1475        0  |    79.19  | Test1_contig1    Test2_Contig1
    10501    11000 |        0        0  |      500        0  |    95.28  | Test1_contig1    Test2_Contig1
      11      300  |        0        0  |      290        0  |    97.98  | Test1_contig2    Test2_Contig1
     501      550  |        0        0  |       50        0  |    96.61  | Test1_contig2    Test2_Contig1
    19001    19999 |        0        0  |      999        0  |    86.55  | Test1_contig2    Test2_Contig1
    1001     9000  |        0        0  |     8000        0  |    85.71  | Test1_contig3    Test2_Contig2
       1     4000  |        0        0  |     4000        0  |    95.73  | Test1_contig4    Test2_Contig2
    3501     4500  |        0        0  |     1000        0  |    96.58  | Test1_contig4    Test2_Contig2
    2001     3000  |        0        0  |     1000        0  |    95.29  | Test1_contig4    Test2_Contig2
    2501     3100  |        0        0  |      600        0  |    85.06  | Test1_contig4    Test2_Contig2
    8001     9000  |        0        0  |     1000        0  |    86.45  | Test1_contig5    Test2_Contig1
    8101     8600  |        0        0  |      500        0  |    95.10  | Test1_contig5    Test2_Contig3
    8501     8700  |        0        0  |      200        0  |    99.75  | Test1_contig5    Test2_Contig3
    8601     9100  |        0        0  |      500        0  |    99.68  | Test1_contig5    Test2_Contig3
    20001    20100 |        0        0  |      100        0  |    99.68  | Test1_contig5    Test2_Contig3
    20051    20150 |        0        0  |      100        0  |    96.80  | Test1_contig5    Test2_Contig3
    20011    20220 |        0        0  |      200        0  |    83.33  | Test1_contig5    Test2_Contig3

`$ cat coords_files/Test1.vs.Test3.coords `  

    /path/to/fastafile/genome_Test1.fa /path/to/fastafile/genome_Test3.fa
    NUCMER

    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
    =====================================================================================
    1001     3100  |        0        0  |     2100        0  |    89.47  | Test1_contig1    Test3_Contig1
    3001     5000  |        0        0  |     2000        0  |    97.69  | Test1_contig1    Test3_Contig1
    8501    10000  |        0        0  |     1500        0  |    79.19  | Test1_contig1    Test3_Contig1
    10001    11000 |        0        0  |     1000        0  |    95.28  | Test1_contig1    Test3_Contig1
     501      550  |        0        0  |       50        0  |    96.61  | Test1_contig2    Test3_Contig1
    19001    20000 |        0        0  |     1000        0  |    86.55  | Test1_contig2    Test3_Contig1
    1001     9000  |        0        0  |     8000        0  |    85.71  | Test1_contig3    Test3_Contig4
       1     4000  |        0        0  |     4000        0  |    95.73  | Test1_contig4    Test3_Contig4
    3501     4500  |        0        0  |     1000        0  |    96.58  | Test1_contig4    Test3_Contig4
    2001     3000  |        0        0  |     1000        0  |    95.29  | Test1_contig4    Test3_Contig4
    2501     2700  |        0        0  |      200        0  |    85.06  | Test1_contig4    Test3_Contig2
    8001     9000  |        0        0  |     1000        0  |    86.45  | Test1_contig5    Test3_Contig1
    8101     8600  |        0        0  |      500        0  |    95.10  | Test1_contig5    Test3_Contig3
    8501     8700  |        0        0  |      200        0  |    99.75  | Test1_contig5    Test3_Contig3
    8601     9100  |        0        0  |      500        0  |    99.68  | Test1_contig5    Test3_Contig3
    20001    20100 |        0        0  |      100        0  |    99.68  | Test1_contig5    Test3_Contig3
    20051    20150 |        0        0  |      100        0  |    96.80  | Test1_contig5    Test3_Contig2
    20001    22000 |        0        0  |     2000        0  |    93.33  | Test1_contig5    Test3_Contig3

These are fake files and only contain information that is used in the scripts that classify contigs as core or accessory.  

If we run the classification with default parameters, that is we only include alignments that are longer than 1 kb and in which aligned segments are more than 90% identical, consider a contig as 'core' if at least 75% of the contig is present in at least 90% of the other genomes. In this case, with only two other genomes, it has to be present in all.

If we filter the alignments we are left with the following regions:  
`$ awk '{if ($10 >= 90 && $7 >= 1000) print }' coords_files/Test1.vs.Test2.coords`   

    3001     9000  |        0        0  |     6000        0  |    97.69  | Test1_contig1    Test2_Contig1
       1     4000  |        0        0  |     4000        0  |    95.73  | Test1_contig4    Test2_Contig2
    3501     4500  |        0        0  |     1000        0  |    96.58  | Test1_contig4    Test2_Contig2
    2001     3000  |        0        0  |     1000        0  |    95.29  | Test1_contig4    Test2_Contig2

`$ awk '{if ($10 >= 90 && $7 >= 1000) print }' coords_files/Test1.vs.Test3.coords`   

    3001     5000  |        0        0  |     2000        0  |    97.69  | Test1_contig1    Test3_Contig1
    10001    11000 |        0        0  |     1000        0  |    95.28  | Test1_contig1    Test3_Contig1
       1     4000  |        0        0  |     4000        0  |    95.73  | Test1_contig4    Test3_Contig4
    3501     4500  |        0        0  |     1000        0  |    96.58  | Test1_contig4    Test3_Contig4
    2001     3000  |        0        0  |     1000        0  |    95.29  | Test1_contig4    Test3_Contig4
    20001    22000 |        0        0  |     2000        0  |    93.33  | Test1_contig5    Test3_Contig3

This means only contig4 would be considered a core contig and
contig1, 3000-5000  
contig4, 1-4500
are considered core regions. Hence Test1 has 6500 core bp and 72000 - 6500 = 65500 accessory bp. Note that the offset on the MUMmer output is 1, hence -like e.g. BLAST- MUMmer starts counting at 1, we correct for this because bedtools starts counting at 0, and would interpret output like `1     4000  |        0        0  |     4000        0  |    95.73  | Test1_contig4    Test3_Contig4` as though the first base does not align.  

Run it like this:
`python ../classify_contigs.py -genome_file Test1.genome_file.tab -ref_name Test1 -in_dir coords_files/ -out_dir OUTPUT`

    Test1 6500 65500

If we now include more alignments, say any alignment longer than 100 bp in which the aligned segments are more than 80% identical, we include all except four alignments and the cumulative size of core regions is much bigger.
For example for contig1 we have:

    2001     3100  |        0        0  |     1100        0  |    89.47  | Test1_contig1    Test2_Contig1
    3001     9000  |        0        0  |     6000        0  |    97.69  | Test1_contig1    Test2_Contig1
    10501    11000 |        0        0  |      500        0  |    95.28  | Test1_contig1    Test2_Contig1

and  

    1001     3100  |        0        0  |     2100        0  |    89.47  | Test1_contig1    Test3_Contig1
    3001     5000  |        0        0  |     2000        0  |    97.69  | Test1_contig1    Test3_Contig1
    10001    11000 |        0        0  |     1000        0  |    95.28  | Test1_contig1    Test3_Contig1

hence 2001-5000 and 10501-11000 are core on contig1 (3500 bp). Contig2 has 999 core bp, contig3 has 8000 core bp and contig4 4500 and contig5 3000. This adds up to
Run it like this:  
`python ../classify_contigs.py -genome_file Test1.genome_file.tab -ref_name Test1 -in_dir coords_files/ -out_dir OUTPUT -min_identity 80 -min_length 100`  

    Test1 18319 53681
