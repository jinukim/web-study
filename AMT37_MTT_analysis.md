

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import subprocess
import seaborn as sns
%matplotlib inline
```

# Trimming using trim_galore


```python
filelist = !ls *l.fastq.gz
filelist
```




    ['v3l.fastq.gz', 'v4l.fastq.gz', 'v5l.fastq.gz', 'v8l.fastq.gz']




```python
#forfile = filelist[0:][::2]
#revfile = filelist[1:][::2]
```


```python
forfile = filelist
forfile
```




    ['v3l.fastq.gz', 'v4l.fastq.gz', 'v5l.fastq.gz', 'v8l.fastq.gz']




```python
#trimming fastq using trim_galore

for i in range(len(forfile)):
    file1 = forfile[i]
#    file2 = revfile[i]
#    subprocess.Popen(['trim_galore','--paired',file1,file2], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#    subprocess.Popen(['trim_galore',file1], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print '!trim_galore %s' %file1
```

    !trim_galore v3l.fastq.gz
    !trim_galore v4l.fastq.gz
    !trim_galore v5l.fastq.gz
    !trim_galore v8l.fastq.gz



```python
!trim_galore v3l.fastq.gz
!trim_galore v4l.fastq.gz
!trim_galore v5l.fastq.gz
!trim_galore v8l.fastq.gz
```

    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)
    
    Path to Cutadapt set as: 'cutadapt' (default)
    1.11
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    
    
    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> v3l.fastq.gz <<)
    
    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	529362	AGATCGGAAGAGC	1000000	52.94
    Nextera	2	CTGTCTCTTATA	1000000	0.00
    smallRNA	0	TGGAATTCTCGG	1000000	0.00
    Using Illumina adapter for trimming (count: 529362). Second best hit was Nextera (count: 2)
    
    
    gzip: stdout: Broken pipe
    Writing report to 'v3l.fastq.gz_trimming_report.txt'
    
    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: v3l.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.4.2
    Cutadapt version: 1.11
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Output file(s) will be GZIP compressed
    
    Writing final adapter and quality trimmed output to v3l_trimmed.fq.gz
    
    
      >>> Now performing quality (cutoff 20) and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file v3l.fastq.gz <<< 
    10000000 sequences processed
    20000000 sequences processed
    30000000 sequences processed
    40000000 sequences processed
    50000000 sequences processed
    This is cutadapt 1.11 with Python 2.7.12
    Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC v3l.fastq.gz
    Trimming 1 adapter with at most 10.0% errors in single-end mode ...
    Finished in 704.42 s (12 us/read; 5.08 M reads/minute).
    
    === Summary ===
    
    Total reads processed:              59,699,262
    Reads with adapters:                33,546,898 (56.2%)
    Reads written (passing filters):    59,699,262 (100.0%)
    
    Total basepairs processed: 8,954,889,300 bp
    Quality-trimmed:           1,177,398,820 bp (13.1%)
    Total written (filtered):  6,814,104,950 bp (76.1%)
    
    === Adapter 1 ===
    
    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 33546898 times.
    
    No. of allowed errors:
    0-9 bp: 0; 10-13 bp: 1
    
    Bases preceding removed adapters:
      A: 15.0%
      C: 36.5%
      G: 33.0%
      T: 14.9%
      none/other: 0.6%
    
    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	5959070	14924815.5	0	5959070
    2	1430612	3731203.9	0	1430612
    3	631367	932801.0	0	631367
    4	162864	233200.2	0	162864
    5	318705	58300.1	0	318705
    6	133597	14575.0	0	133597
    7	549439	3643.8	0	549439
    8	232336	910.9	0	232336
    9	1045620	227.7	0	1024879 20741
    10	558344	56.9	1	478735 79609
    11	350384	14.2	1	289683 60701
    12	292455	3.6	1	246170 46285
    13	607440	0.9	1	512188 95252
    14	300567	0.9	1	253308 47259
    15	319121	0.9	1	274164 44957
    16	357450	0.9	1	299349 58101
    17	465305	0.9	1	399545 65760
    18	339566	0.9	1	290180 49386
    19	242613	0.9	1	212086 30527
    20	346972	0.9	1	308984 37988
    21	194226	0.9	1	173116 21110
    22	190732	0.9	1	170567 20165
    23	179719	0.9	1	152512 27207
    24	601623	0.9	1	527017 74606
    25	428200	0.9	1	374075 54125
    26	276114	0.9	1	230419 45695
    27	687776	0.9	1	584680 103096
    28	999630	0.9	1	877690 121940
    29	461827	0.9	1	406372 55455
    30	492769	0.9	1	438694 54075
    31	481795	0.9	1	436617 45178
    32	634666	0.9	1	581632 53034
    33	269105	0.9	1	242673 26432
    34	637758	0.9	1	586107 51651
    35	262421	0.9	1	237983 24438
    36	88850	0.9	1	74738 14112
    37	610150	0.9	1	561670 48480
    38	354471	0.9	1	321807 32664
    39	579592	0.9	1	530936 48656
    40	249442	0.9	1	217399 32043
    41	537681	0.9	1	496793 40888
    42	370615	0.9	1	338115 32500
    43	558617	0.9	1	527395 31222
    44	216502	0.9	1	196812 19690
    45	158746	0.9	1	144011 14735
    46	171104	0.9	1	157142 13962
    47	379943	0.9	1	355440 24503
    48	206753	0.9	1	190013 16740
    49	374958	0.9	1	353570 21388
    50	230273	0.9	1	212771 17502
    51	298201	0.9	1	279451 18750
    52	155735	0.9	1	143266 12469
    53	135833	0.9	1	120697 15136
    54	601375	0.9	1	568175 33200
    55	314737	0.9	1	295154 19583
    56	267146	0.9	1	251241 15905
    57	294458	0.9	1	279518 14940
    58	117237	0.9	1	106446 10791
    59	303266	0.9	1	288170 15096
    60	244466	0.9	1	232195 12271
    61	161642	0.9	1	151638 10004
    62	135428	0.9	1	126651 8777
    63	188855	0.9	1	179326 9529
    64	74449	0.9	1	67328 7121
    65	252895	0.9	1	242059 10836
    66	125452	0.9	1	118132 7320
    67	45324	0.9	1	38835 6489
    68	210589	0.9	1	198363 12226
    69	231112	0.9	1	217049 14063
    70	234016	0.9	1	219886 14130
    71	225704	0.9	1	213584 12120
    72	236137	0.9	1	223690 12447
    73	245514	0.9	1	233096 12418
    74	256352	0.9	1	243565 12787
    75	290344	0.9	1	274927 15417
    76	415663	0.9	1	401125 14538
    77	231133	0.9	1	220579 10554
    78	126750	0.9	1	119597 7153
    79	86769	0.9	1	80782 5987
    80	67800	0.9	1	62993 4807
    81	52604	0.9	1	49124 3480
    82	45396	0.9	1	42442 2954
    83	38457	0.9	1	35979 2478
    84	34445	0.9	1	32337 2108
    85	30064	0.9	1	28160 1904
    86	26528	0.9	1	24899 1629
    87	23916	0.9	1	22450 1466
    88	21211	0.9	1	19922 1289
    89	19347	0.9	1	18170 1177
    90	16419	0.9	1	15296 1123
    91	15079	0.9	1	14125 954
    92	13705	0.9	1	12832 873
    93	12482	0.9	1	11665 817
    94	10707	0.9	1	10005 702
    95	9712	0.9	1	9091 621
    96	9007	0.9	1	8457 550
    97	7942	0.9	1	7500 442
    98	6868	0.9	1	6442 426
    99	6029	0.9	1	5690 339
    100	5397	0.9	1	5105 292
    101	4408	0.9	1	4183 225
    102	4319	0.9	1	4100 219
    103	3770	0.9	1	3595 175
    104	3123	0.9	1	2971 152
    105	2959	0.9	1	2813 146
    106	2455	0.9	1	2331 124
    107	2052	0.9	1	1950 102
    108	1886	0.9	1	1816 70
    109	1668	0.9	1	1599 69
    110	1440	0.9	1	1376 64
    111	1200	0.9	1	1144 56
    112	1073	0.9	1	1014 59
    113	848	0.9	1	818 30
    114	852	0.9	1	824 28
    115	654	0.9	1	631 23
    116	468	0.9	1	450 18
    117	411	0.9	1	399 12
    118	316	0.9	1	301 15
    119	342	0.9	1	330 12
    120	265	0.9	1	248 17
    121	131	0.9	1	111 20
    122	131	0.9	1	124 7
    123	106	0.9	1	98 8
    124	56	0.9	1	55 1
    125	45	0.9	1	41 4
    126	36	0.9	1	33 3
    127	35	0.9	1	33 2
    128	29	0.9	1	26 3
    129	29	0.9	1	28 1
    130	24	0.9	1	22 2
    131	14	0.9	1	11 3
    132	24	0.9	1	15 9
    133	41	0.9	1	12 29
    134	9	0.9	1	8 1
    135	15	0.9	1	14 1
    136	30	0.9	1	20 10
    137	31	0.9	1	7 24
    138	4	0.9	1	2 2
    139	11	0.9	1	5 6
    140	7	0.9	1	3 4
    141	10	0.9	1	7 3
    142	5	0.9	1	3 2
    143	2	0.9	1	0 2
    144	9	0.9	1	4 5
    145	8	0.9	1	1 7
    146	10	0.9	1	1 9
    147	17	0.9	1	6 11
    148	34	0.9	1	1 33
    149	118	0.9	1	1 117
    150	1786	0.9	1	5 1781
    
    
    RUN STATISTICS FOR INPUT FILE: v3l.fastq.gz
    =============================================
    59699262 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	180200 (0.3%)
    
    
    gzip: stdout: Broken pipe
    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)
    
    Path to Cutadapt set as: 'cutadapt' (default)
    1.11
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    
    
    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> v4l.fastq.gz <<)
    
    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	409297	AGATCGGAAGAGC	1000000	40.93
    smallRNA	0	TGGAATTCTCGG	1000000	0.00
    Nextera	0	CTGTCTCTTATA	1000000	0.00
    Using Illumina adapter for trimming (count: 409297). Second best hit was smallRNA (count: 0)
    
    
    gzip: stdout: Broken pipe
    Writing report to 'v4l.fastq.gz_trimming_report.txt'
    
    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: v4l.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.4.2
    Cutadapt version: 1.11
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Output file(s) will be GZIP compressed
    
    Writing final adapter and quality trimmed output to v4l_trimmed.fq.gz
    
    
      >>> Now performing quality (cutoff 20) and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file v4l.fastq.gz <<< 
    10000000 sequences processed
    20000000 sequences processed
    30000000 sequences processed
    40000000 sequences processed
    50000000 sequences processed
    This is cutadapt 1.11 with Python 2.7.12
    Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC v4l.fastq.gz
    Trimming 1 adapter with at most 10.0% errors in single-end mode ...
    Finished in 615.01 s (12 us/read; 5.03 M reads/minute).
    
    === Summary ===
    
    Total reads processed:              51,528,832
    Reads with adapters:                26,212,830 (50.9%)
    Reads written (passing filters):    51,528,832 (100.0%)
    
    Total basepairs processed: 7,729,324,800 bp
    Quality-trimmed:             899,996,763 bp (11.6%)
    Total written (filtered):  6,189,804,487 bp (80.1%)
    
    === Adapter 1 ===
    
    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 26212830 times.
    
    No. of allowed errors:
    0-9 bp: 0; 10-13 bp: 1
    
    Bases preceding removed adapters:
      A: 20.1%
      C: 32.0%
      G: 29.4%
      T: 17.9%
      none/other: 0.7%
    
    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	6722844	12882208.0	0	6722844
    2	1442698	3220552.0	0	1442698
    3	628156	805138.0	0	628156
    4	167988	201284.5	0	167988
    5	251197	50321.1	0	251197
    6	113359	12580.3	0	113359
    7	412687	3145.1	0	412687
    8	157348	786.3	0	157348
    9	754475	196.6	0	738323 16152
    10	368461	49.1	1	314727 53734
    11	281200	12.3	1	236302 44898
    12	225910	3.1	1	193081 32829
    13	423552	0.8	1	361510 62042
    14	216347	0.8	1	185321 31026
    15	257883	0.8	1	225397 32486
    16	235332	0.8	1	200716 34616
    17	342188	0.8	1	299475 42713
    18	236867	0.8	1	206436 30431
    19	192771	0.8	1	171988 20783
    20	271347	0.8	1	245697 25650
    21	172812	0.8	1	155662 17150
    22	152511	0.8	1	138018 14493
    23	139726	0.8	1	120300 19426
    24	440872	0.8	1	392974 47898
    25	304366	0.8	1	269986 34380
    26	155999	0.8	1	130396 25603
    27	380492	0.8	1	325940 54552
    28	589295	0.8	1	523862 65433
    29	277013	0.8	1	247605 29408
    30	194814	0.8	1	172950 21864
    31	296577	0.8	1	272461 24116
    32	357854	0.8	1	332562 25292
    33	187427	0.8	1	172330 15097
    34	281811	0.8	1	260146 21665
    35	233139	0.8	1	214133 19006
    36	250663	0.8	1	222188 28475
    37	492776	0.8	1	458711 34065
    38	286201	0.8	1	262605 23596
    39	382729	0.8	1	352654 30075
    40	160480	0.8	1	137388 23092
    41	463037	0.8	1	433539 29498
    42	287841	0.8	1	266711 21130
    43	453284	0.8	1	431164 22120
    44	158481	0.8	1	146166 12315
    45	130430	0.8	1	120009 10421
    46	165680	0.8	1	155299 10381
    47	254818	0.8	1	240592 14226
    48	139050	0.8	1	128751 10299
    49	269209	0.8	1	255931 13278
    50	151340	0.8	1	141172 10168
    51	202930	0.8	1	191825 11105
    52	100226	0.8	1	93236 6990
    53	90593	0.8	1	81834 8759
    54	373469	0.8	1	355876 17593
    55	193274	0.8	1	182772 10502
    56	174503	0.8	1	165703 8800
    57	183387	0.8	1	175029 8358
    58	74688	0.8	1	68367 6321
    59	201559	0.8	1	193110 8449
    60	156712	0.8	1	149710 7002
    61	106689	0.8	1	101247 5442
    62	86370	0.8	1	81627 4743
    63	121067	0.8	1	115670 5397
    64	46139	0.8	1	42237 3902
    65	165807	0.8	1	159654 6153
    66	78720	0.8	1	74760 3960
    67	26765	0.8	1	23773 2992
    68	123328	0.8	1	117382 5946
    69	129430	0.8	1	122552 6878
    70	131858	0.8	1	125249 6609
    71	129737	0.8	1	123764 5973
    72	132920	0.8	1	126649 6271
    73	136281	0.8	1	130213 6068
    74	141164	0.8	1	135302 5862
    75	168399	0.8	1	160941 7458
    76	248631	0.8	1	241110 7521
    77	140575	0.8	1	135328 5247
    78	78322	0.8	1	74663 3659
    79	53268	0.8	1	50080 3188
    80	41655	0.8	1	39104 2551
    81	32445	0.8	1	30566 1879
    82	26843	0.8	1	25338 1505
    83	23126	0.8	1	21889 1237
    84	20240	0.8	1	19170 1070
    85	17627	0.8	1	16686 941
    86	15777	0.8	1	14965 812
    87	14226	0.8	1	13454 772
    88	12464	0.8	1	11811 653
    89	11299	0.8	1	10731 568
    90	9795	0.8	1	9288 507
    91	8634	0.8	1	8136 498
    92	7699	0.8	1	7286 413
    93	6882	0.8	1	6550 332
    94	6285	0.8	1	5939 346
    95	5657	0.8	1	5356 301
    96	5312	0.8	1	4995 317
    97	4514	0.8	1	4279 235
    98	4153	0.8	1	3949 204
    99	3700	0.8	1	3524 176
    100	3063	0.8	1	2923 140
    101	2930	0.8	1	2810 120
    102	2283	0.8	1	2185 98
    103	2164	0.8	1	2078 86
    104	1875	0.8	1	1798 77
    105	1724	0.8	1	1650 74
    106	1339	0.8	1	1285 54
    107	1093	0.8	1	1045 48
    108	1046	0.8	1	1003 43
    109	959	0.8	1	916 43
    110	780	0.8	1	760 20
    111	627	0.8	1	601 26
    112	557	0.8	1	538 19
    113	448	0.8	1	427 21
    114	400	0.8	1	383 17
    115	393	0.8	1	384 9
    116	322	0.8	1	313 9
    117	260	0.8	1	236 24
    118	166	0.8	1	160 6
    119	175	0.8	1	163 12
    120	113	0.8	1	109 4
    121	118	0.8	1	112 6
    122	70	0.8	1	56 14
    123	83	0.8	1	80 3
    124	70	0.8	1	50 20
    125	37	0.8	1	33 4
    126	36	0.8	1	27 9
    127	26	0.8	1	22 4
    128	20	0.8	1	18 2
    129	9	0.8	1	8 1
    130	15	0.8	1	14 1
    131	12	0.8	1	7 5
    132	6	0.8	1	4 2
    133	13	0.8	1	11 2
    134	23	0.8	1	11 12
    135	20	0.8	1	16 4
    136	13	0.8	1	6 7
    137	19	0.8	1	9 10
    138	5	0.8	1	3 2
    139	10	0.8	1	7 3
    140	11	0.8	1	8 3
    141	4	0.8	1	2 2
    142	9	0.8	1	7 2
    143	1	0.8	1	0 1
    144	10	0.8	1	3 7
    145	9	0.8	1	6 3
    146	4	0.8	1	3 1
    147	2	0.8	1	1 1
    148	14	0.8	1	1 13
    149	61	0.8	1	0 61
    150	902	0.8	1	1 901
    
    
    RUN STATISTICS FOR INPUT FILE: v4l.fastq.gz
    =============================================
    51528832 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	148861 (0.3%)
    
    
    gzip: stdout: Broken pipe
    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)
    
    Path to Cutadapt set as: 'cutadapt' (default)
    1.11
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    
    
    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> v5l.fastq.gz <<)
    
    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	342735	AGATCGGAAGAGC	1000000	34.27
    Nextera	1	CTGTCTCTTATA	1000000	0.00
    smallRNA	1	TGGAATTCTCGG	1000000	0.00
    Using Illumina adapter for trimming (count: 342735). Second best hit was Nextera (count: 1)
    
    
    gzip: stdout: Broken pipe
    Writing report to 'v5l.fastq.gz_trimming_report.txt'
    
    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: v5l.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.4.2
    Cutadapt version: 1.11
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Output file(s) will be GZIP compressed
    
    Writing final adapter and quality trimmed output to v5l_trimmed.fq.gz
    
    
      >>> Now performing quality (cutoff 20) and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file v5l.fastq.gz <<< 
    10000000 sequences processed
    20000000 sequences processed
    30000000 sequences processed
    40000000 sequences processed
    This is cutadapt 1.11 with Python 2.7.12
    Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC v5l.fastq.gz
    Trimming 1 adapter with at most 10.0% errors in single-end mode ...
    Finished in 554.81 s (12 us/read; 5.13 M reads/minute).
    
    === Summary ===
    
    Total reads processed:              47,395,082
    Reads with adapters:                21,794,862 (46.0%)
    Reads written (passing filters):    47,395,082 (100.0%)
    
    Total basepairs processed: 7,109,262,300 bp
    Quality-trimmed:             838,588,542 bp (11.8%)
    Total written (filtered):  5,786,999,629 bp (81.4%)
    
    === Adapter 1 ===
    
    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 21794862 times.
    
    No. of allowed errors:
    0-9 bp: 0; 10-13 bp: 1
    
    Bases preceding removed adapters:
      A: 20.7%
      C: 31.7%
      G: 29.8%
      T: 17.0%
      none/other: 0.7%
    
    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	6449934	11848770.5	0	6449934
    2	1438721	2962192.6	0	1438721
    3	572879	740548.2	0	572879
    4	136477	185137.0	0	136477
    5	210848	46284.3	0	210848
    6	88531	11571.1	0	88531
    7	340796	2892.8	0	340796
    8	131795	723.2	0	131795
    9	580717	180.8	0	568576 12141
    10	291097	45.2	1	250518 40579
    11	222901	11.3	1	187616 35285
    12	182375	2.8	1	156024 26351
    13	327951	0.7	1	280570 47381
    14	167616	0.7	1	143752 23864
    15	197135	0.7	1	172317 24818
    16	185434	0.7	1	158011 27423
    17	265550	0.7	1	232004 33546
    18	183666	0.7	1	159844 23822
    19	147920	0.7	1	131395 16525
    20	212821	0.7	1	192543 20278
    21	128916	0.7	1	116175 12741
    22	117349	0.7	1	106383 10966
    23	107597	0.7	1	92555 15042
    24	339667	0.7	1	301998 37669
    25	239134	0.7	1	211205 27929
    26	134753	0.7	1	112185 22568
    27	336726	0.7	1	288024 48702
    28	474958	0.7	1	421244 53714
    29	227049	0.7	1	202417 24632
    30	181714	0.7	1	160941 20773
    31	258694	0.7	1	236084 22610
    32	295117	0.7	1	273238 21879
    33	168026	0.7	1	152946 15080
    34	242020	0.7	1	222767 19253
    35	188822	0.7	1	172786 16036
    36	132232	0.7	1	119967 12265
    37	311558	0.7	1	290010 21548
    38	199937	0.7	1	183586 16351
    39	330173	0.7	1	304775 25398
    40	83076	0.7	1	68954 14122
    41	260028	0.7	1	242491 17537
    42	208449	0.7	1	192713 15736
    43	330124	0.7	1	313586 16538
    44	118332	0.7	1	108875 9457
    45	97952	0.7	1	90171 7781
    46	103757	0.7	1	96421 7336
    47	199869	0.7	1	188431 11438
    48	108563	0.7	1	100950 7613
    49	192767	0.7	1	182745 10022
    50	116180	0.7	1	108402 7778
    51	153432	0.7	1	144989 8443
    52	76823	0.7	1	71107 5716
    53	66571	0.7	1	59898 6673
    54	285246	0.7	1	271515 13731
    55	148715	0.7	1	140348 8367
    56	130194	0.7	1	123053 7141
    57	143440	0.7	1	136648 6792
    58	56173	0.7	1	51326 4847
    59	150585	0.7	1	143911 6674
    60	120328	0.7	1	114883 5445
    61	79410	0.7	1	74927 4483
    62	65467	0.7	1	61576 3891
    63	89912	0.7	1	85761 4151
    64	34736	0.7	1	31642 3094
    65	123168	0.7	1	118454 4714
    66	58632	0.7	1	55497 3135
    67	20813	0.7	1	18134 2679
    68	94322	0.7	1	89228 5094
    69	98678	0.7	1	93091 5587
    70	99524	0.7	1	94368 5156
    71	97938	0.7	1	93224 4714
    72	98817	0.7	1	93853 4964
    73	103160	0.7	1	98321 4839
    74	106426	0.7	1	101608 4818
    75	124166	0.7	1	118164 6002
    76	181370	0.7	1	175419 5951
    77	103083	0.7	1	98884 4199
    78	55975	0.7	1	53008 2967
    79	39274	0.7	1	36642 2632
    80	29472	0.7	1	27395 2077
    81	23392	0.7	1	21921 1471
    82	19517	0.7	1	18422 1095
    83	16742	0.7	1	15783 959
    84	14607	0.7	1	13752 855
    85	13066	0.7	1	12367 699
    86	11532	0.7	1	10918 614
    87	10386	0.7	1	9806 580
    88	9132	0.7	1	8602 530
    89	8138	0.7	1	7693 445
    90	7105	0.7	1	6662 443
    91	6498	0.7	1	6121 377
    92	5765	0.7	1	5453 312
    93	5034	0.7	1	4733 301
    94	4499	0.7	1	4261 238
    95	4113	0.7	1	3860 253
    96	3572	0.7	1	3362 210
    97	3133	0.7	1	2981 152
    98	2962	0.7	1	2775 187
    99	2611	0.7	1	2472 139
    100	2149	0.7	1	2050 99
    101	1996	0.7	1	1899 97
    102	1845	0.7	1	1763 82
    103	1393	0.7	1	1332 61
    104	1355	0.7	1	1298 57
    105	1209	0.7	1	1163 46
    106	1060	0.7	1	1010 50
    107	872	0.7	1	843 29
    108	760	0.7	1	730 30
    109	685	0.7	1	664 21
    110	625	0.7	1	598 27
    111	525	0.7	1	498 27
    112	384	0.7	1	364 20
    113	417	0.7	1	391 26
    114	337	0.7	1	321 16
    115	297	0.7	1	284 13
    116	240	0.7	1	223 17
    117	199	0.7	1	175 24
    118	190	0.7	1	181 9
    119	116	0.7	1	111 5
    120	111	0.7	1	100 11
    121	79	0.7	1	74 5
    122	82	0.7	1	76 6
    123	76	0.7	1	71 5
    124	57	0.7	1	54 3
    125	43	0.7	1	32 11
    126	35	0.7	1	33 2
    127	25	0.7	1	25
    128	16	0.7	1	15 1
    129	31	0.7	1	24 7
    130	11	0.7	1	7 4
    131	9	0.7	1	9
    132	12	0.7	1	11 1
    133	9	0.7	1	6 3
    134	12	0.7	1	11 1
    135	15	0.7	1	5 10
    136	14	0.7	1	9 5
    137	8	0.7	1	7 1
    138	10	0.7	1	6 4
    139	9	0.7	1	6 3
    140	14	0.7	1	3 11
    141	6	0.7	1	5 1
    142	3	0.7	1	1 2
    143	5	0.7	1	3 2
    144	7	0.7	1	5 2
    145	12	0.7	1	4 8
    146	17	0.7	1	2 15
    147	9	0.7	1	0 9
    148	17	0.7	1	1 16
    149	72	0.7	1	0 72
    150	1025	0.7	1	2 1023
    
    
    RUN STATISTICS FOR INPUT FILE: v5l.fastq.gz
    =============================================
    47395082 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	142789 (0.3%)
    
    
    gzip: stdout: Broken pipe
    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)
    
    Path to Cutadapt set as: 'cutadapt' (default)
    1.11
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    
    
    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> v8l.fastq.gz <<)
    
    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	331279	AGATCGGAAGAGC	1000000	33.13
    Nextera	0	CTGTCTCTTATA	1000000	0.00
    smallRNA	0	TGGAATTCTCGG	1000000	0.00
    Using Illumina adapter for trimming (count: 331279). Second best hit was Nextera (count: 0)
    
    
    gzip: stdout: Broken pipe
    Writing report to 'v8l.fastq.gz_trimming_report.txt'
    
    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: v8l.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.4.2
    Cutadapt version: 1.11
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Output file(s) will be GZIP compressed
    
    Writing final adapter and quality trimmed output to v8l_trimmed.fq.gz
    
    
      >>> Now performing quality (cutoff 20) and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file v8l.fastq.gz <<< 
    10000000 sequences processed
    20000000 sequences processed
    30000000 sequences processed
    40000000 sequences processed
    This is cutadapt 1.11 with Python 2.7.12
    Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC v8l.fastq.gz
    Trimming 1 adapter with at most 10.0% errors in single-end mode ...
    Finished in 521.50 s (12 us/read; 5.04 M reads/minute).
    
    === Summary ===
    
    Total reads processed:              43,774,358
    Reads with adapters:                19,772,755 (45.2%)
    Reads written (passing filters):    43,774,358 (100.0%)
    
    Total basepairs processed: 6,566,153,700 bp
    Quality-trimmed:             775,961,737 bp (11.8%)
    Total written (filtered):  5,361,120,824 bp (81.6%)
    
    === Adapter 1 ===
    
    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 19772755 times.
    
    No. of allowed errors:
    0-9 bp: 0; 10-13 bp: 1
    
    Bases preceding removed adapters:
      A: 20.7%
      C: 31.8%
      G: 29.9%
      T: 16.8%
      none/other: 0.8%
    
    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	5979618	10943589.5	0	5979618
    2	1329701	2735897.4	0	1329701
    3	512504	683974.3	0	512504
    4	123795	170993.6	0	123795
    5	188931	42748.4	0	188931
    6	78828	10687.1	0	78828
    7	303130	2671.8	0	303130
    8	114210	667.9	0	114210
    9	545687	167.0	0	534623 11064
    10	265442	41.7	1	228244 37198
    11	196822	10.4	1	165584 31238
    12	157149	2.6	1	134353 22796
    13	302938	0.7	1	258648 44290
    14	148317	0.7	1	126992 21325
    15	179663	0.7	1	156571 23092
    16	166389	0.7	1	141364 25025
    17	239431	0.7	1	208624 30807
    18	165811	0.7	1	144052 21759
    19	130014	0.7	1	115335 14679
    20	192033	0.7	1	173519 18514
    21	116525	0.7	1	104559 11966
    22	105768	0.7	1	95391 10377
    23	97263	0.7	1	83385 13878
    24	311774	0.7	1	276692 35082
    25	216975	0.7	1	191749 25226
    26	125476	0.7	1	105114 20362
    27	302047	0.7	1	258813 43234
    28	460716	0.7	1	409052 51664
    29	207036	0.7	1	184702 22334
    30	178641	0.7	1	159270 19371
    31	236665	0.7	1	216559 20106
    32	280720	0.7	1	259780 20940
    33	168745	0.7	1	154392 14353
    34	226534	0.7	1	209211 17323
    35	125459	0.7	1	115086 10373
    36	92148	0.7	1	83043 9105
    37	293930	0.7	1	273443 20487
    38	196513	0.7	1	180917 15596
    39	265109	0.7	1	243916 21193
    40	81497	0.7	1	68578 12919
    41	259878	0.7	1	242416 17462
    42	198488	0.7	1	183607 14881
    43	356773	0.7	1	340307 16466
    44	120524	0.7	1	111105 9419
    45	98743	0.7	1	91628 7115
    46	63149	0.7	1	57909 5240
    47	170419	0.7	1	160603 9816
    48	86029	0.7	1	79299 6730
    49	170736	0.7	1	161825 8911
    50	86504	0.7	1	79883 6621
    51	128774	0.7	1	121259 7515
    52	63697	0.7	1	58841 4856
    53	54828	0.7	1	48766 6062
    54	241172	0.7	1	229196 11976
    55	122525	0.7	1	115231 7294
    56	109048	0.7	1	103121 5927
    57	118979	0.7	1	113206 5773
    58	43767	0.7	1	39555 4212
    59	122776	0.7	1	117181 5595
    60	99964	0.7	1	95292 4672
    61	66060	0.7	1	62293 3767
    62	52590	0.7	1	49368 3222
    63	73837	0.7	1	70256 3581
    64	26538	0.7	1	24026 2512
    65	104662	0.7	1	100478 4184
    66	50075	0.7	1	47318 2757
    67	15931	0.7	1	13809 2122
    68	82844	0.7	1	78639 4205
    69	86874	0.7	1	82090 4784
    70	88305	0.7	1	83573 4732
    71	88189	0.7	1	83942 4247
    72	88196	0.7	1	83689 4507
    73	90610	0.7	1	86403 4207
    74	95093	0.7	1	90801 4292
    75	109334	0.7	1	103857 5477
    76	160908	0.7	1	155644 5264
    77	89701	0.7	1	86071 3630
    78	48227	0.7	1	45805 2422
    79	32624	0.7	1	30715 1909
    80	25073	0.7	1	23422 1651
    81	21036	0.7	1	19662 1374
    82	17639	0.7	1	16538 1101
    83	14712	0.7	1	13759 953
    84	12914	0.7	1	12109 805
    85	11830	0.7	1	11153 677
    86	10245	0.7	1	9681 564
    87	9181	0.7	1	8691 490
    88	8053	0.7	1	7606 447
    89	7540	0.7	1	7131 409
    90	6225	0.7	1	5869 356
    91	5586	0.7	1	5275 311
    92	5072	0.7	1	4802 270
    93	4357	0.7	1	4076 281
    94	3925	0.7	1	3678 247
    95	3567	0.7	1	3352 215
    96	3222	0.7	1	3018 204
    97	2850	0.7	1	2687 163
    98	2482	0.7	1	2312 170
    99	2397	0.7	1	2279 118
    100	2004	0.7	1	1922 82
    101	1706	0.7	1	1629 77
    102	1581	0.7	1	1513 68
    103	1440	0.7	1	1376 64
    104	1326	0.7	1	1253 73
    105	1154	0.7	1	1118 36
    106	905	0.7	1	872 33
    107	850	0.7	1	812 38
    108	759	0.7	1	730 29
    109	610	0.7	1	590 20
    110	638	0.7	1	605 33
    111	493	0.7	1	475 18
    112	370	0.7	1	357 13
    113	380	0.7	1	340 40
    114	320	0.7	1	300 20
    115	319	0.7	1	303 16
    116	224	0.7	1	213 11
    117	226	0.7	1	213 13
    118	145	0.7	1	134 11
    119	147	0.7	1	136 11
    120	155	0.7	1	143 12
    121	103	0.7	1	92 11
    122	62	0.7	1	57 5
    123	77	0.7	1	65 12
    124	51	0.7	1	49 2
    125	33	0.7	1	30 3
    126	26	0.7	1	20 6
    127	29	0.7	1	21 8
    128	12	0.7	1	11 1
    129	13	0.7	1	11 2
    130	10	0.7	1	9 1
    131	17	0.7	1	8 9
    132	13	0.7	1	8 5
    133	23	0.7	1	10 13
    134	13	0.7	1	9 4
    135	9	0.7	1	7 2
    136	13	0.7	1	4 9
    137	15	0.7	1	8 7
    138	15	0.7	1	6 9
    139	8	0.7	1	3 5
    140	8	0.7	1	2 6
    141	9	0.7	1	4 5
    142	7	0.7	1	4 3
    143	5	0.7	1	4 1
    144	6	0.7	1	3 3
    145	4	0.7	1	2 2
    146	5	0.7	1	1 4
    147	9	0.7	1	1 8
    148	16	0.7	1	0 16
    149	75	0.7	1	2 73
    150	1086	0.7	1	2 1084
    
    
    RUN STATISTICS FOR INPUT FILE: v8l.fastq.gz
    =============================================
    43774358 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	131127 (0.3%)
    
    
    gzip: stdout: Broken pipe



```python
!ls -lh *fq.gz
```

    -rw-rw-r-- 1 jk jk 3.7G Jan  4 12:19 v3l_trimmed.fq.gz
    -rw-rw-r-- 1 jk jk 3.5G Jan  4 12:36 v4l_trimmed.fq.gz
    -rw-rw-r-- 1 jk jk 3.2G Jan  4 12:52 v5l_trimmed.fq.gz
    -rw-rw-r-- 1 jk jk 3.0G Jan  4 13:06 v8l_trimmed.fq.gz


# Briefly check composition using kraken


```python
#testfiles = !ls *1.fq.gz
testfiles = !ls *.fq.gz
testfiles
```




    ['v3l_trimmed.fq.gz',
     'v4l_trimmed.fq.gz',
     'v5l_trimmed.fq.gz',
     'v8l_trimmed.fq.gz']




```python
dbloc = '/home/jk/Downloads/kraken_db'
```


```python
#import shlex
```


```python
#Run kraken db

for i in range(len(testfiles)):
    targetfile = testfiles[i]
    targetname = targetfile + '.kraken'
    logfile = targetname + '.log'
    command_line = 'kraken --db %s %s > %s 2> %s' %(dbloc,targetfile,targetname,logfile)
#    args = shlex.split(command_line)
    print command_line
    subprocess.Popen(command_line,shell=True)
```

    kraken --db /home/jk/Downloads/kraken_db amt-control-covaris_S10_L001_R1_001_val_1.fq.gz > amt-control-covaris_S10_L001_R1_001_val_1.fq.gz.kraken 2> amt-control-covaris_S10_L001_R1_001_val_1.fq.gz.kraken.log
    kraken --db /home/jk/Downloads/kraken_db amt-control-ex_S9_L001_R1_001_val_1.fq.gz > amt-control-ex_S9_L001_R1_001_val_1.fq.gz.kraken 2> amt-control-ex_S9_L001_R1_001_val_1.fq.gz.kraken.log
    kraken --db /home/jk/Downloads/kraken_db niams37-v3-ll_S5_L001_R1_001_val_1.fq.gz > niams37-v3-ll_S5_L001_R1_001_val_1.fq.gz.kraken 2> niams37-v3-ll_S5_L001_R1_001_val_1.fq.gz.kraken.log
    kraken --db /home/jk/Downloads/kraken_db niams37-v4-ll_S6_L001_R1_001_val_1.fq.gz > niams37-v4-ll_S6_L001_R1_001_val_1.fq.gz.kraken 2> niams37-v4-ll_S6_L001_R1_001_val_1.fq.gz.kraken.log
    kraken --db /home/jk/Downloads/kraken_db niams37-v5-ll_S7_L001_R1_001_val_1.fq.gz > niams37-v5-ll_S7_L001_R1_001_val_1.fq.gz.kraken 2> niams37-v5-ll_S7_L001_R1_001_val_1.fq.gz.kraken.log
    kraken --db /home/jk/Downloads/kraken_db niams37-v8-ll_S8_L001_R1_001_val_1.fq.gz > niams37-v8-ll_S8_L001_R1_001_val_1.fq.gz.kraken 2> niams37-v8-ll_S8_L001_R1_001_val_1.fq.gz.kraken.log



```python
!head *kraken.log
```

    ==> amt-control-covaris_S10_L001_R1_001_val_1.fq.gz.kraken.log <==
    
    gzip: stdout: Broken pipe
    773 sequences (0.06 Mbp) processed in 4.811s (9.6 Kseq/m, 0.75 Mbp/m).
      162 sequences classified (20.96%)
      611 sequences unclassified (79.04%)
    
    ==> amt-control-ex_S9_L001_R1_001_val_1.fq.gz.kraken.log <==
    
    gzip: stdout: Broken pipe
    44784 sequences (5.98 Mbp) processed in 47.925s (56.1 Kseq/m, 7.49 Mbp/m).
      5906 sequences classified (13.19%)
      38878 sequences unclassified (86.81%)
    
    ==> niams37-v3-ll_S5_L001_R1_001_val_1.fq.gz.kraken.log <==
    
    gzip: stdout: Broken pipe
    35557 sequences (4.29 Mbp) processed in 60.319s (35.4 Kseq/m, 4.27 Mbp/m).
      3041 sequences classified (8.55%)
      32516 sequences unclassified (91.45%)
    
    ==> niams37-v4-ll_S6_L001_R1_001_val_1.fq.gz.kraken.log <==
    
    gzip: stdout: Broken pipe
    102766 sequences (13.19 Mbp) processed in 102.593s (60.1 Kseq/m, 7.71 Mbp/m).
      22312 sequences classified (21.71%)
      80454 sequences unclassified (78.29%)
    
    ==> niams37-v5-ll_S7_L001_R1_001_val_1.fq.gz.kraken.log <==
    
    gzip: stdout: Broken pipe
    50554 sequences (6.65 Mbp) processed in 73.320s (41.4 Kseq/m, 5.44 Mbp/m).
      8239 sequences classified (16.30%)
      42315 sequences unclassified (83.70%)
    
    ==> niams37-v8-ll_S8_L001_R1_001_val_1.fq.gz.kraken.log <==
    
    gzip: stdout: Broken pipe
    77424 sequences (10.23 Mbp) processed in 86.748s (53.6 Kseq/m, 7.08 Mbp/m).
      13351 sequences classified (17.24%)
      64073 sequences unclassified (82.76%)



```python
#Run kraken translate

for i in range(len(testfiles)):
    targetfile = testfiles[i] + '.kraken'
    labelfile = testfiles[i] + '.labels'
    command_line = 'kraken-translate --db %s %s > %s' %(dbloc,targetfile,labelfile)
    print command_line
    subprocess.Popen(command_line,shell=True)
```

    kraken-translate --db /home/jk/Downloads/kraken_db amt-control-covaris_S10_L001_R1_001_val_1.fq.gz.kraken > amt-control-covaris_S10_L001_R1_001_val_1.fq.gz.labels
    kraken-translate --db /home/jk/Downloads/kraken_db amt-control-ex_S9_L001_R1_001_val_1.fq.gz.kraken > amt-control-ex_S9_L001_R1_001_val_1.fq.gz.labels
    kraken-translate --db /home/jk/Downloads/kraken_db niams37-v3-ll_S5_L001_R1_001_val_1.fq.gz.kraken > niams37-v3-ll_S5_L001_R1_001_val_1.fq.gz.labels
    kraken-translate --db /home/jk/Downloads/kraken_db niams37-v4-ll_S6_L001_R1_001_val_1.fq.gz.kraken > niams37-v4-ll_S6_L001_R1_001_val_1.fq.gz.labels
    kraken-translate --db /home/jk/Downloads/kraken_db niams37-v5-ll_S7_L001_R1_001_val_1.fq.gz.kraken > niams37-v5-ll_S7_L001_R1_001_val_1.fq.gz.labels
    kraken-translate --db /home/jk/Downloads/kraken_db niams37-v8-ll_S8_L001_R1_001_val_1.fq.gz.kraken > niams37-v8-ll_S8_L001_R1_001_val_1.fq.gz.labels



```python
labelfiles = !ls *labels
labelfiles
```




    ['amt-control-covaris_S10_L001_R1_001_val_1.fq.gz.labels',
     'amt-control-ex_S9_L001_R1_001_val_1.fq.gz.labels',
     'niams37-v3-ll_S5_L001_R1_001_val_1.fq.gz.labels',
     'niams37-v4-ll_S6_L001_R1_001_val_1.fq.gz.labels',
     'niams37-v5-ll_S7_L001_R1_001_val_1.fq.gz.labels',
     'niams37-v8-ll_S8_L001_R1_001_val_1.fq.gz.labels']




```python
!head *.labels
```

    ==> amt-control-covaris_S10_L001_R1_001_val_1.fq.gz.labels <==
    M01048:154:000000000-D2Y86:1:1101:18404:2309	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas stutzeri group;Pseudomonas stutzeri subgroup;Pseudomonas stutzeri
    M01048:154:000000000-D2Y86:1:1101:15969:2823	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:21956:3197	root;cellular organisms;Bacteria
    M01048:154:000000000-D2Y86:1:1101:22701:4167	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:14497:4547	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas pertucinogena group;Pseudomonas denitrificans;Pseudomonas denitrificans ATCC 13867
    M01048:154:000000000-D2Y86:1:1101:15468:4841	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:21078:5642	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:15733:6586	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:17056:6855	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:6642:6955	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    
    ==> amt-control-ex_S9_L001_R1_001_val_1.fq.gz.labels <==
    M01048:154:000000000-D2Y86:1:1101:14975:1399	root;cellular organisms;Bacteria
    M01048:154:000000000-D2Y86:1:1101:13743:1479	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia;Escherichia coli
    M01048:154:000000000-D2Y86:1:1101:13045:1607	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas
    M01048:154:000000000-D2Y86:1:1101:15496:1607	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:13337:1619	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:12970:1620	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:13902:1628	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas fluorescens group;Pseudomonas fluorescens;Pseudomonas fluorescens A506
    M01048:154:000000000-D2Y86:1:1101:12799:1641	root;cellular organisms;Archaea
    M01048:154:000000000-D2Y86:1:1101:15710:1659	root;cellular organisms;Bacteria
    M01048:154:000000000-D2Y86:1:1101:14297:1716	root;cellular organisms;Archaea
    
    ==> niams37-v3-ll_S5_L001_R1_001_val_1.fq.gz.labels <==
    M01048:154:000000000-D2Y86:1:1101:15697:1599	root;cellular organisms;Bacteria;Firmicutes
    M01048:154:000000000-D2Y86:1:1101:12479:1712	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas
    M01048:154:000000000-D2Y86:1:1101:17488:1764	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria
    M01048:154:000000000-D2Y86:1:1101:16418:1806	root;Viruses;dsRNA viruses;Endornaviridae;Endornavirus;Oryza sativa endornavirus
    M01048:154:000000000-D2Y86:1:1101:13305:1838	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:19249:1885	root
    M01048:154:000000000-D2Y86:1:1101:17178:1969	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria
    M01048:154:000000000-D2Y86:1:1101:15391:1969	root;cellular organisms;Bacteria
    M01048:154:000000000-D2Y86:1:1101:17322:1988	root;cellular organisms;Archaea
    M01048:154:000000000-D2Y86:1:1101:12927:2009	root;cellular organisms;Archaea
    
    ==> niams37-v4-ll_S6_L001_R1_001_val_1.fq.gz.labels <==
    M01048:154:000000000-D2Y86:1:1101:17079:1400	root
    M01048:154:000000000-D2Y86:1:1101:14631:1415	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus
    M01048:154:000000000-D2Y86:1:1101:15100:1415	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:15852:1415	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas
    M01048:154:000000000-D2Y86:1:1101:16759:1468	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:17071:1478	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:15144:1508	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus
    M01048:154:000000000-D2Y86:1:1101:15642:1530	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus
    M01048:154:000000000-D2Y86:1:1101:15710:1534	root;cellular organisms;Archaea
    M01048:154:000000000-D2Y86:1:1101:17467:1542	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus
    
    ==> niams37-v5-ll_S7_L001_R1_001_val_1.fq.gz.labels <==
    M01048:154:000000000-D2Y86:1:1101:14826:1412	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus
    M01048:154:000000000-D2Y86:1:1101:14099:1512	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:16267:1597	root;cellular organisms;Archaea
    M01048:154:000000000-D2Y86:1:1101:13266:1673	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas
    M01048:154:000000000-D2Y86:1:1101:17099:1673	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:12497:1709	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas sp. TKP
    M01048:154:000000000-D2Y86:1:1101:12985:1724	root
    M01048:154:000000000-D2Y86:1:1101:18808:1784	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:17085:1784	root;cellular organisms;Bacteria
    M01048:154:000000000-D2Y86:1:1101:13019:1788	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae
    
    ==> niams37-v8-ll_S8_L001_R1_001_val_1.fq.gz.labels <==
    M01048:154:000000000-D2Y86:1:1101:15557:1366	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:15917:1399	root;cellular organisms;Archaea
    M01048:154:000000000-D2Y86:1:1101:17058:1407	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:16000:1408	root;cellular organisms;Bacteria;Firmicutes
    M01048:154:000000000-D2Y86:1:1101:16006:1434	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae
    M01048:154:000000000-D2Y86:1:1101:14613:1445	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus
    M01048:154:000000000-D2Y86:1:1101:15488:1466	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    M01048:154:000000000-D2Y86:1:1101:15862:1469	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus
    M01048:154:000000000-D2Y86:1:1101:14618:1487	root;cellular organisms;Archaea
    M01048:154:000000000-D2Y86:1:1101:17732:1601	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae


# Remove human and rRNA  

## remove human


```python
#Run bowtie2 against human genome

for i in range(len(testfiles)):
    humandbloc = '/home/jk/Dropbox/UCSD/Research/Data/genome_reference/human/humanh19'
    targetfile = testfiles[i]
    samfile = testfiles[i] + '.human.sam'
    logfile = testfiles[i] + '.human.bowtie2log'
    
    command_line = 'bowtie2 -p 4 -x %s -U %s -S %s 2> %s' %(humandbloc,targetfile,samfile,logfile)
    print "!" + command_line
#    subprocess.Popen(command_line,shell=True)
```

    !bowtie2 -p 4 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/human/humanh19 -U v3l_trimmed.fq.gz -S v3l_trimmed.fq.gz.human.sam 2> v3l_trimmed.fq.gz.human.bowtie2log
    !bowtie2 -p 4 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/human/humanh19 -U v4l_trimmed.fq.gz -S v4l_trimmed.fq.gz.human.sam 2> v4l_trimmed.fq.gz.human.bowtie2log
    !bowtie2 -p 4 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/human/humanh19 -U v5l_trimmed.fq.gz -S v5l_trimmed.fq.gz.human.sam 2> v5l_trimmed.fq.gz.human.bowtie2log
    !bowtie2 -p 4 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/human/humanh19 -U v8l_trimmed.fq.gz -S v8l_trimmed.fq.gz.human.sam 2> v8l_trimmed.fq.gz.human.bowtie2log



```python
!bowtie2 -p 4 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/human/humanh19 -U v3l_trimmed.fq.gz -S v3l_trimmed.fq.gz.human.sam 2> v3l_trimmed.fq.gz.human.bowtie2log
!bowtie2 -p 4 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/human/humanh19 -U v4l_trimmed.fq.gz -S v4l_trimmed.fq.gz.human.sam 2> v4l_trimmed.fq.gz.human.bowtie2log
!bowtie2 -p 4 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/human/humanh19 -U v5l_trimmed.fq.gz -S v5l_trimmed.fq.gz.human.sam 2> v5l_trimmed.fq.gz.human.bowtie2log
!bowtie2 -p 4 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/human/humanh19 -U v8l_trimmed.fq.gz -S v8l_trimmed.fq.gz.human.sam 2> v8l_trimmed.fq.gz.human.bowtie2log
```


```python
!ls -lh *sam
```

    -rw-rw-r-- 1 jk jk 20G Jan  4 14:06 v3l_trimmed.fq.gz.human.sam
    -rw-rw-r-- 1 jk jk 18G Jan  4 14:40 v4l_trimmed.fq.gz.human.sam
    -rw-rw-r-- 1 jk jk 17G Jan  4 15:16 v5l_trimmed.fq.gz.human.sam
    -rw-rw-r-- 1 jk jk 16G Jan  4 15:48 v8l_trimmed.fq.gz.human.sam



```python
#Check contamination
!head *bowtie2log
```

    ==> v3l_trimmed.fq.gz.human.bowtie2log <==
    59519062 reads; of these:
      59519062 (100.00%) were unpaired; of these:
        7775925 (13.06%) aligned 0 times
        1624874 (2.73%) aligned exactly 1 time
        50118263 (84.21%) aligned >1 times
    86.94% overall alignment rate
    
    ==> v4l_trimmed.fq.gz.human.bowtie2log <==
    51379971 reads; of these:
      51379971 (100.00%) were unpaired; of these:
        12647229 (24.62%) aligned 0 times
        2097681 (4.08%) aligned exactly 1 time
        36635061 (71.30%) aligned >1 times
    75.38% overall alignment rate
    
    ==> v5l_trimmed.fq.gz.human.bowtie2log <==
    47252293 reads; of these:
      47252293 (100.00%) were unpaired; of these:
        9104185 (19.27%) aligned 0 times
        2806174 (5.94%) aligned exactly 1 time
        35341934 (74.79%) aligned >1 times
    80.73% overall alignment rate
    
    ==> v8l_trimmed.fq.gz.human.bowtie2log <==
    43643231 reads; of these:
      43643231 (100.00%) were unpaired; of these:
        8782422 (20.12%) aligned 0 times
        1368550 (3.14%) aligned exactly 1 time
        33492259 (76.74%) aligned >1 times
    79.88% overall alignment rate



```python
#convert human mapped files to bam files

for i in range(len(testfiles)):
    targetfile = testfiles[i]
    samfile = testfiles[i] + '.human.sam'
    bamfile = testfiles[i] + '.human.sam.bam'
    
    command_line = 'samtools view -bS %s >%s' %(samfile, bamfile)
    print "!" + command_line
#    subprocess.Popen(command_line,shell=True)
```

    `samtools view -bS v3l_trimmed.fq.gz.human.sam >v3l_trimmed.fq.gz.human.sam.bam`;
    `samtools view -bS v4l_trimmed.fq.gz.human.sam >v4l_trimmed.fq.gz.human.sam.bam`;
    `samtools view -bS v5l_trimmed.fq.gz.human.sam >v5l_trimmed.fq.gz.human.sam.bam`;
    `samtools view -bS v8l_trimmed.fq.gz.human.sam >v8l_trimmed.fq.gz.human.sam.bam`;



```python
!samtools view -bS v3l_trimmed.fq.gz.human.sam >v3l_trimmed.fq.gz.human.sam.bam
!samtools view -bS v4l_trimmed.fq.gz.human.sam >v4l_trimmed.fq.gz.human.sam.bam
!samtools view -bS v5l_trimmed.fq.gz.human.sam >v5l_trimmed.fq.gz.human.sam.bam
!samtools view -bS v8l_trimmed.fq.gz.human.sam >v8l_trimmed.fq.gz.human.sam.bam
```

    [samopen] SAM header is present: 25 sequences.
    [samopen] SAM header is present: 25 sequences.
    [samopen] SAM header is present: 25 sequences.
    [samopen] SAM header is present: 25 sequences.



```python
!ls -lh *.human.sam.bam
```

    -rw-rw-r-- 1 jk jk 5.0G Jan  4 16:03 v3l_trimmed.fq.gz.human.sam.bam
    -rw-rw-r-- 1 jk jk 4.6G Jan  4 16:17 v4l_trimmed.fq.gz.human.sam.bam
    -rw-rw-r-- 1 jk jk 4.3G Jan  4 16:30 v5l_trimmed.fq.gz.human.sam.bam
    -rw-rw-r-- 1 jk jk 3.9G Jan  4 16:41 v8l_trimmed.fq.gz.human.sam.bam



```python
#convert unaligned reads from bam to fastq

for i in range(len(testfiles)):
    fastqout = testfiles[i].split('.')[0] + '.humanremoved.fastq'
    bamfile = testfiles[i] + '.human.sam.bam'
    
    command_line = 'bam2fastq --unaligned --no-aligned --force -o %s %s' %(fastqout, bamfile)
    print "`" + command_line + "`;"
#    subprocess.Popen(command_line,shell=True)
```

    `bam2fastq --unaligned --no-aligned --force -o v3l_trimmed.humanremoved.fastq v3l_trimmed.fq.gz.human.sam.bam`;
    `bam2fastq --unaligned --no-aligned --force -o v4l_trimmed.humanremoved.fastq v4l_trimmed.fq.gz.human.sam.bam`;
    `bam2fastq --unaligned --no-aligned --force -o v5l_trimmed.humanremoved.fastq v5l_trimmed.fq.gz.human.sam.bam`;
    `bam2fastq --unaligned --no-aligned --force -o v8l_trimmed.humanremoved.fastq v8l_trimmed.fq.gz.human.sam.bam`;



```python
!bam2fastq --unaligned --no-aligned --force -o v3l_trimmed.humanremoved.fastq v3l_trimmed.fq.gz.human.sam.bam
!bam2fastq --unaligned --no-aligned --force -o v4l_trimmed.humanremoved.fastq v4l_trimmed.fq.gz.human.sam.bam
!bam2fastq --unaligned --no-aligned --force -o v5l_trimmed.humanremoved.fastq v5l_trimmed.fq.gz.human.sam.bam
!bam2fastq --unaligned --no-aligned --force -o v8l_trimmed.humanremoved.fastq v8l_trimmed.fq.gz.human.sam.bam
```

    This looks like single-end data from lane 676.
    Output will be in v3l_trimmed.humanremoved.fastq
    59519062 sequences in the BAM file
    7775925 sequences exported
    This looks like single-end data from lane 676.
    Output will be in v4l_trimmed.humanremoved.fastq
    51379971 sequences in the BAM file
    12647229 sequences exported
    This looks like single-end data from lane 676.
    Output will be in v5l_trimmed.humanremoved.fastq
    47252293 sequences in the BAM file
    9104185 sequences exported
    This looks like single-end data from lane 676.
    Output will be in v8l_trimmed.humanremoved.fastq
    43643231 sequences in the BAM file
    8782422 sequences exported



```python
nohuman = !ls -1 *fastq
nohuman
```




    ['v3l_trimmed.humanremoved.fastq',
     'v4l_trimmed.humanremoved.fastq',
     'v5l_trimmed.humanremoved.fastq',
     'v8l_trimmed.humanremoved.fastq']



## remove rRNA


```python
#Check rRNA contamination of human-unaligned reads
for i in range(len(nohuman)):
    humandbloc = '/home/jk/Dropbox/UCSD/Research/Data/genome_reference/rRNA/silva'
    targetfile = nohuman[i]
    samfile = nohuman[i] + '.norrna.sam'
    logfile = nohuman[i] + '.norrna.bowtie2log'
    
    command_line = 'bowtie2 -x %s -U %s -S %s 2> %s' %(humandbloc,targetfile,samfile,logfile)
    print "!" + command_line
```

    `bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/rRNA/silva -U v3l_trimmed.humanremoved.fastq -S v3l_trimmed.humanremoved.fastq.norrna.sam 2> v3l_trimmed.humanremoved.fastq.norrna.bowtie2log`;
    `bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/rRNA/silva -U v4l_trimmed.humanremoved.fastq -S v4l_trimmed.humanremoved.fastq.norrna.sam 2> v4l_trimmed.humanremoved.fastq.norrna.bowtie2log`;
    `bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/rRNA/silva -U v5l_trimmed.humanremoved.fastq -S v5l_trimmed.humanremoved.fastq.norrna.sam 2> v5l_trimmed.humanremoved.fastq.norrna.bowtie2log`;
    `bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/rRNA/silva -U v8l_trimmed.humanremoved.fastq -S v8l_trimmed.humanremoved.fastq.norrna.sam 2> v8l_trimmed.humanremoved.fastq.norrna.bowtie2log`;



```python
!bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/rRNA/silva -U v3l_trimmed.humanremoved.fastq -S v3l_trimmed.humanremoved.fastq.norrna.sam 2> v3l_trimmed.humanremoved.fastq.norrna.bowtie2log
!bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/rRNA/silva -U v4l_trimmed.humanremoved.fastq -S v4l_trimmed.humanremoved.fastq.norrna.sam 2> v4l_trimmed.humanremoved.fastq.norrna.bowtie2log
!bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/rRNA/silva -U v5l_trimmed.humanremoved.fastq -S v5l_trimmed.humanremoved.fastq.norrna.sam 2> v5l_trimmed.humanremoved.fastq.norrna.bowtie2log
!bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/rRNA/silva -U v8l_trimmed.humanremoved.fastq -S v8l_trimmed.humanremoved.fastq.norrna.sam 2> v8l_trimmed.humanremoved.fastq.norrna.bowtie2log
```


```python
!ls -la *.norrna.sam
```


```python
!head *norrna.bowtie2log
```


```python
#convert human mapped files to bam files

for i in range(len(nohuman)):
    targetfile = nohuman[i]
    samfile = nohuman[i] + '.norrna.sam'
    bamfile = nohuman[i] + '.norrna.sam.bam'
    
    command_line = 'samtools view -bS %s >%s' %(samfile, bamfile)
    print "`" + command_line + "`;"
    subprocess.Popen(command_line,shell=True)
```

    `samtools view -bS v3l_trimmed.humanremoved.fastq.norrna.sam >v3l_trimmed.humanremoved.fastq.norrna.sam.bam`;
    `samtools view -bS v4l_trimmed.humanremoved.fastq.norrna.sam >v4l_trimmed.humanremoved.fastq.norrna.sam.bam`;
    `samtools view -bS v5l_trimmed.humanremoved.fastq.norrna.sam >v5l_trimmed.humanremoved.fastq.norrna.sam.bam`;
    `samtools view -bS v8l_trimmed.humanremoved.fastq.norrna.sam >v8l_trimmed.humanremoved.fastq.norrna.sam.bam`;



```python
!samtools view -bS v3l_trimmed.humanremoved.fastq.norrna.sam >v3l_trimmed.humanremoved.fastq.norrna.sam.bam
!samtools view -bS v4l_trimmed.humanremoved.fastq.norrna.sam >v4l_trimmed.humanremoved.fastq.norrna.sam.bam
!samtools view -bS v5l_trimmed.humanremoved.fastq.norrna.sam >v5l_trimmed.humanremoved.fastq.norrna.sam.bam
!samtools view -bS v8l_trimmed.humanremoved.fastq.norrna.sam >v8l_trimmed.humanremoved.fastq.norrna.sam.bam
```


```python
!ls -lh *sam.bam
```


```python
#convert human-unalinged reads and norRNA to fastq files

for i in range(len(nohuman)):
    fastqout = nohuman[i].split('.')[0] + '.nohuman.norrna.fastq'
    bamfile = nohuman[i] + '.norrna.sam.bam'
    
    command_line = 'bam2fastq --unaligned --no-aligned --force -o %s %s' %(fastqout, bamfile)
    print "`" + command_line + "`;"
#    subprocess.Popen(command_line,shell=True)

```

    `bam2fastq --unaligned --no-aligned --force -o v3l_trimmed.nohuman.norrna.fastq v3l_trimmed.humanremoved.fastq.norrna.sam.bam`;
    `bam2fastq --unaligned --no-aligned --force -o v4l_trimmed.nohuman.norrna.fastq v4l_trimmed.humanremoved.fastq.norrna.sam.bam`;
    `bam2fastq --unaligned --no-aligned --force -o v5l_trimmed.nohuman.norrna.fastq v5l_trimmed.humanremoved.fastq.norrna.sam.bam`;
    `bam2fastq --unaligned --no-aligned --force -o v8l_trimmed.nohuman.norrna.fastq v8l_trimmed.humanremoved.fastq.norrna.sam.bam`;



```python
!bam2fastq --unaligned --no-aligned --force -o v3l_trimmed.nohuman.norrna.fastq v3l_trimmed.humanremoved.fastq.norrna.sam.bam
!bam2fastq --unaligned --no-aligned --force -o v4l_trimmed.nohuman.norrna.fastq v4l_trimmed.humanremoved.fastq.norrna.sam.bam
!bam2fastq --unaligned --no-aligned --force -o v5l_trimmed.nohuman.norrna.fastq v5l_trimmed.humanremoved.fastq.norrna.sam.bam
!bam2fastq --unaligned --no-aligned --force -o v8l_trimmed.nohuman.norrna.fastq v8l_trimmed.humanremoved.fastq.norrna.sam.bam
```

# Run kraken using human-rRNA removed fastq


```python
#Run kraken db

for i in range(len(testfiles)):
    dbloc = '/home/jk/Downloads/kraken_db'
    targetfile = testfiles[i].split('.')[0] + '.nohuman.norrna.fastq'
    targetname = targetfile + '.kraken'
    logfile = targetname + '.log'
    command_line = 'kraken --db %s %s > %s 2> %s' %(dbloc,targetfile,targetname,logfile)
#    args = shlex.split(command_line)
    print "`" + command_line + "`;"
#    subprocess.Popen(command_line,shell=True)
```

    `kraken --db /home/jk/Downloads/kraken_db v3l_trimmed.nohuman.norrna.fastq > v3l_trimmed.nohuman.norrna.fastq.kraken 2> v3l_trimmed.nohuman.norrna.fastq.kraken.log`;
    `kraken --db /home/jk/Downloads/kraken_db v4l_trimmed.nohuman.norrna.fastq > v4l_trimmed.nohuman.norrna.fastq.kraken 2> v4l_trimmed.nohuman.norrna.fastq.kraken.log`;
    `kraken --db /home/jk/Downloads/kraken_db v5l_trimmed.nohuman.norrna.fastq > v5l_trimmed.nohuman.norrna.fastq.kraken 2> v5l_trimmed.nohuman.norrna.fastq.kraken.log`;
    `kraken --db /home/jk/Downloads/kraken_db v8l_trimmed.nohuman.norrna.fastq > v8l_trimmed.nohuman.norrna.fastq.kraken 2> v8l_trimmed.nohuman.norrna.fastq.kraken.log`;



```python
!head *nohuman.norrna.fastq.kraken.log
```

    ==> v3l_trimmed.nohuman.norrna.fastq.kraken.log <==
    11677547 sequences (1487.69 Mbp) processed in 689.191s (1016.6 Kseq/m, 129.52 Mbp/m).
      9821248 sequences classified (84.10%)
      1856299 sequences unclassified (15.90%)
    
    ==> v4l_trimmed.nohuman.norrna.fastq.kraken.log <==
    9182260 sequences (1173.80 Mbp) processed in 300.803s (1831.5 Kseq/m, 234.13 Mbp/m).
      7615342 sequences classified (82.94%)
      1566918 sequences unclassified (17.06%)
    
    ==> v5l_trimmed.nohuman.norrna.fastq.kraken.log <==
    3235361 sequences (412.21 Mbp) processed in 124.350s (1561.1 Kseq/m, 198.89 Mbp/m).
      1831586 sequences classified (56.61%)
      1403775 sequences unclassified (43.39%)
    
    ==> v8l_trimmed.nohuman.norrna.fastq.kraken.log <==
    4028653 sequences (516.41 Mbp) processed in 126.526s (1910.4 Kseq/m, 244.89 Mbp/m).
      2193612 sequences classified (54.45%)
      1835041 sequences unclassified (45.55%)



```python
#Run kraken translate

for i in range(len(testfiles)):
    targetfile = testfiles[i].split('.')[0] + '.nohuman.norrna.fastq.kraken'
    labelfile = testfiles[i].split('.')[0] + '.nohuman.norrna.labels'
    command_line = 'kraken-translate --db %s %s > %s' %(dbloc,targetfile,labelfile)
    print "`" + command_line + "`;"
    subprocess.Popen(command_line,shell=True)
```

    `kraken-translate --db /home/jk/Downloads/kraken_db v3l_trimmed.nohuman.norrna.fastq.kraken > v3l_trimmed.nohuman.norrna.labels`;
    `kraken-translate --db /home/jk/Downloads/kraken_db v4l_trimmed.nohuman.norrna.fastq.kraken > v4l_trimmed.nohuman.norrna.labels`;
    `kraken-translate --db /home/jk/Downloads/kraken_db v5l_trimmed.nohuman.norrna.fastq.kraken > v5l_trimmed.nohuman.norrna.labels`;
    `kraken-translate --db /home/jk/Downloads/kraken_db v8l_trimmed.nohuman.norrna.fastq.kraken > v8l_trimmed.nohuman.norrna.labels`;



```python
!head *nohuman.norrna.labels
```

    ==> v3l_trimmed.nohuman.norrna.labels <==
    D00361:676:H3T5YBCX2:1:1101:2546:3473	root;cellular organisms;Bacteria
    D00361:676:H3T5YBCX2:1:1101:2942:3487	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:2882:3470	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:3164:3488	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:3574:3399	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus
    D00361:676:H3T5YBCX2:1:1101:3659:3441	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:3583:3478	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:4485:3427	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:5202:3422	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:5304:3470	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    
    ==> v4l_trimmed.nohuman.norrna.labels <==
    D00361:676:H3T5YBCX2:1:1101:1723:3393	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:1999:3489	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:2456:3426	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:2612:3475	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:2460:3485	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:3132:3475	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:3198:3446	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:5785:3420	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:6100:3487	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus
    D00361:676:H3T5YBCX2:1:1101:8038:3435	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus
    
    ==> v5l_trimmed.nohuman.norrna.labels <==
    D00361:676:H3T5YBCX2:1:1101:3804:3452	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:6463:3496	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus
    D00361:676:H3T5YBCX2:1:1101:7801:3402	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:9611:3439	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:10757:3472	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae
    D00361:676:H3T5YBCX2:1:1101:11557:3420	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:12129:3407	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus
    D00361:676:H3T5YBCX2:1:1101:15974:3395	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:19799:3445	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:1533:3592	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas sp. TKP
    
    ==> v8l_trimmed.nohuman.norrna.labels <==
    D00361:676:H3T5YBCX2:1:1101:1630:3422	root;cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas
    D00361:676:H3T5YBCX2:1:1101:2369:3413	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:6883:3428	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus
    D00361:676:H3T5YBCX2:1:1101:8320:3493	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:11181:3500	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:11550:3474	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:13741:3427	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:14292:3457	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:15087:3491	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132
    D00361:676:H3T5YBCX2:1:1101:16784:3434	root;cellular organisms;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus;Staphylococcus aureus subsp. aureus;Staphylococcus aureus subsp. aureus MSHR1132



```python
!tail -n 1 *human.bowtie2log
```

    ==> v3l_trimmed.fq.gz.human.bowtie2log <==
    58.71% overall alignment rate
    
    ==> v4l_trimmed.fq.gz.human.bowtie2log <==
    63.95% overall alignment rate
    
    ==> v5l_trimmed.fq.gz.human.bowtie2log <==
    68.35% overall alignment rate
    
    ==> v8l_trimmed.fq.gz.human.bowtie2log <==
    76.05% overall alignment rate



```python
#human RNA contamination
humanRNA = pd.DataFrame({'v3ll':87.89, 'v4ll':76.53, 'v5ll':81.62, 'v8ll':80.50, 
                         'extractioncontrol':83.98, 'covariscontrol':6.34},
                        index=['human'])
humanRNA = humanRNA.transpose()
humanRNA['nonhuman']=100-humanRNA['human']
humanRNA
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>human</th>
      <th>nonhuman</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>covariscontrol</th>
      <td>6.34</td>
      <td>93.66</td>
    </tr>
    <tr>
      <th>extractioncontrol</th>
      <td>83.98</td>
      <td>16.02</td>
    </tr>
    <tr>
      <th>v3ll</th>
      <td>87.89</td>
      <td>12.11</td>
    </tr>
    <tr>
      <th>v4ll</th>
      <td>76.53</td>
      <td>23.47</td>
    </tr>
    <tr>
      <th>v5ll</th>
      <td>81.62</td>
      <td>18.38</td>
    </tr>
    <tr>
      <th>v8ll</th>
      <td>80.50</td>
      <td>19.50</td>
    </tr>
  </tbody>
</table>
</div>




```python
humanRNA.plot(kind='bar',stacked='True',color=['Grey','#4169e1'])
plt.xticks(rotation=0)
plt.legend(fontsize=10, bbox_to_anchor=(1.3, 1.0))
plt.savefig("./human_contam.png", format='png',bbox_inches='tight')
#plt.savefig("./human_contam.svg", format='svg',bbox_inches='tight')
```


![png](output_46_0.png)



```python

```

# Taxonomic analysis


```python
import pandas as pd
import numpy as np
import matplotlib as pl
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline
```


```python
#Import kraken files into dataframe

krakendict={}

for i in range(len(testfiles)):
    filetag = testfiles[i].split('_')[0]    
#    targetfile = testfiles[i] + '.labels'

    targetfile = testfiles[i].split('.')[0] + '.nohuman.norrna.labels'
    label = pd.read_table("./%s"%targetfile, sep='\t', header=None)
    label.drop(0,axis=1, inplace=True)
    labeldf = pd.DataFrame(label[1].str.split(';').tolist())
    labeldf['sample'] = filetag
    krakendict[filetag] = labeldf
```


```python
krakendict.keys()
```




    ['niams37-v5-ll',
     'amt-control-ex',
     'niams37-v4-ll',
     'niams37-v3-ll',
     'amt-control-covaris',
     'niams37-v8-ll']




```python
#Concatanate all df

krakendf=pd.DataFrame()

for i in krakendict.keys():
    krakendf = pd.concat([krakendf, krakendict[i]])

krakendf.drop([0,1,9,10,11,12], axis=1, inplace=True)
krakendf.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
      <th>7</th>
      <th>8</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Pseudomonadales</td>
      <td>Pseudomonadaceae</td>
      <td>Pseudomonas</td>
      <td>None</td>
      <td>niams37-v5-ll</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Pseudomonadales</td>
      <td>Pseudomonadaceae</td>
      <td>Pseudomonas</td>
      <td>Pseudomonas sp. TKP</td>
      <td>niams37-v5-ll</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Pseudomonadales</td>
      <td>Pseudomonadaceae</td>
      <td>Pseudomonas</td>
      <td>Pseudomonas fluorescens group</td>
      <td>niams37-v5-ll</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Bacteria</td>
      <td>Firmicutes</td>
      <td>Bacilli</td>
      <td>Bacillales</td>
      <td>Staphylococcaceae</td>
      <td>Staphylococcus</td>
      <td>Staphylococcus aureus</td>
      <td>niams37-v5-ll</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Pseudomonadales</td>
      <td>Pseudomonadaceae</td>
      <td>Pseudomonas</td>
      <td>Pseudomonas fluorescens group</td>
      <td>niams37-v5-ll</td>
    </tr>
  </tbody>
</table>
</div>




```python
krakendf.columns=['domain','phylum','class','order','family','genus','species','sample']
krakendf.fillna('unassigned', inplace=True)
krakendf.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>domain</th>
      <th>phylum</th>
      <th>class</th>
      <th>order</th>
      <th>family</th>
      <th>genus</th>
      <th>species</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Pseudomonadales</td>
      <td>Pseudomonadaceae</td>
      <td>Pseudomonas</td>
      <td>unassigned</td>
      <td>niams37-v5-ll</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Pseudomonadales</td>
      <td>Pseudomonadaceae</td>
      <td>Pseudomonas</td>
      <td>Pseudomonas sp. TKP</td>
      <td>niams37-v5-ll</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Pseudomonadales</td>
      <td>Pseudomonadaceae</td>
      <td>Pseudomonas</td>
      <td>Pseudomonas fluorescens group</td>
      <td>niams37-v5-ll</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Bacteria</td>
      <td>Firmicutes</td>
      <td>Bacilli</td>
      <td>Bacillales</td>
      <td>Staphylococcaceae</td>
      <td>Staphylococcus</td>
      <td>Staphylococcus aureus</td>
      <td>niams37-v5-ll</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Bacteria</td>
      <td>Proteobacteria</td>
      <td>Gammaproteobacteria</td>
      <td>Pseudomonadales</td>
      <td>Pseudomonadaceae</td>
      <td>Pseudomonas</td>
      <td>Pseudomonas fluorescens group</td>
      <td>niams37-v5-ll</td>
    </tr>
  </tbody>
</table>
</div>




```python
krakendf['count']=1
```


```python
#krakendf.fillna('None',inplace=True)
```


```python
speciesdf = krakendf.groupby(['sample','species']).sum()   #You can change 'species' to any level you want.
speciesdf.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>count</th>
    </tr>
    <tr>
      <th>sample</th>
      <th>species</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="5" valign="top">amt-control-covaris</th>
      <th>Alteromonas macleodii</th>
      <td>9</td>
    </tr>
    <tr>
      <th>Bradyrhizobium sp. BTAi1</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Corynebacteriaceae</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Dinoroseobacter shibae</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Escherichia coli</th>
      <td>13</td>
    </tr>
  </tbody>
</table>
</div>




```python
speciesdf_un = speciesdf['count'].unstack()
speciesdf_un.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>species</th>
      <th>Achromobacter xylosoxidans</th>
      <th>Acidovorax sp. JS42</th>
      <th>Acidovorax sp. KKS102</th>
      <th>Acinetobacter oleivorans</th>
      <th>Acinetobacter sp. ADP1</th>
      <th>Akkermansia</th>
      <th>Alteromonas macleodii</th>
      <th>Anaerococcus prevotii</th>
      <th>Azospirillum lipoferum</th>
      <th>Azotobacter</th>
      <th>...</th>
      <th>Streptococcus thermophilus</th>
      <th>Treponema primitia</th>
      <th>Variovorax paradoxus</th>
      <th>Xanthobacter autotrophicus</th>
      <th>Xanthomonas albilineans</th>
      <th>Xanthomonas campestris</th>
      <th>Yersinia enterocolitica</th>
      <th>pseudomallei group</th>
      <th>unassigned</th>
      <th>unclassified Enterobacteriaceae (miscellaneous)</th>
    </tr>
    <tr>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>amt-control-covaris</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>9.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>18.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>87.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>amt-control-ex</th>
      <td>14.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>3.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>3.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>355.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>niams37-v3-ll</th>
      <td>6.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>7.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>1.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>177.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>niams37-v4-ll</th>
      <td>17.0</td>
      <td>1.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>7.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>1.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>1.0</td>
      <td>NaN</td>
      <td>549.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>niams37-v5-ll</th>
      <td>7.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>1.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>218.0</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 104 columns</p>
</div>




```python
speciesdf_un.sum(axis=1)
```




    sample
    amt-control-covaris     153.0
    amt-control-ex          814.0
    niams37-v3-ll           466.0
    niams37-v4-ll          7466.0
    niams37-v5-ll           871.0
    niams37-v8-ll          1654.0
    dtype: float64




```python
speciesdf_un.fillna(0, inplace=True)
speciesdf_un['tot'] = speciesdf_un.sum(axis=1)


for i in speciesdf_un.index:
    speciesdf_un.loc[i,:] = speciesdf_un.loc[i,:]/speciesdf_un.loc[i,'tot'] *100
    
speciesratio = speciesdf_un.drop('tot',axis=1)
speciesratio.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>species</th>
      <th>Achromobacter xylosoxidans</th>
      <th>Acidovorax sp. JS42</th>
      <th>Acidovorax sp. KKS102</th>
      <th>Acinetobacter oleivorans</th>
      <th>Acinetobacter sp. ADP1</th>
      <th>Akkermansia</th>
      <th>Alteromonas macleodii</th>
      <th>Anaerococcus prevotii</th>
      <th>Azospirillum lipoferum</th>
      <th>Azotobacter</th>
      <th>...</th>
      <th>Streptococcus thermophilus</th>
      <th>Treponema primitia</th>
      <th>Variovorax paradoxus</th>
      <th>Xanthobacter autotrophicus</th>
      <th>Xanthomonas albilineans</th>
      <th>Xanthomonas campestris</th>
      <th>Yersinia enterocolitica</th>
      <th>pseudomallei group</th>
      <th>unassigned</th>
      <th>unclassified Enterobacteriaceae (miscellaneous)</th>
    </tr>
    <tr>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>amt-control-covaris</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>5.882353</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>11.764706</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>56.862745</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>amt-control-ex</th>
      <td>1.719902</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.368550</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.000000</td>
      <td>0.36855</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>43.611794</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>niams37-v3-ll</th>
      <td>1.287554</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.502146</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.214592</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>37.982833</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>niams37-v4-ll</th>
      <td>0.227699</td>
      <td>0.013394</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.093758</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.013394</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.013394</td>
      <td>0.0</td>
      <td>7.353335</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>niams37-v5-ll</th>
      <td>0.803674</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.229621</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.114811</td>
      <td>0.0</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>25.028703</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 104 columns</p>
</div>




```python
speciesratio['Staphylococcus aureus']
```




    sample
    amt-control-covaris     1.307190
    amt-control-ex          0.859951
    niams37-v3-ll           4.077253
    niams37-v4-ll          87.382802
    niams37-v5-ll          49.024110
    niams37-v8-ll          39.661427
    Name: Staphylococcus aureus, dtype: float64




```python
speciesratio.sum(axis=1)
```




    sample
    amt-control-covaris    100.0
    amt-control-ex         100.0
    niams37-v3-ll          100.0
    niams37-v4-ll          100.0
    niams37-v5-ll          100.0
    niams37-v8-ll          100.0
    dtype: float64




```python
speciesratio.rename(columns={'Propionibacteriaceae':'Probionibacterium acnes'},inplace=True)
```


```python
speciesratio.sum().sort_values(ascending=False).index
```




    Index([u'unassigned', u'Staphylococcus aureus',
           u'Pseudomonas fluorescens group', u'Pseudomonas sp. TKP',
           u'Staphylococcus epidermidis', u'Xanthobacter autotrophicus',
           u'Probionibacterium acnes', u'Escherichia coli',
           u'Alteromonas macleodii', u'Pseudomonas putida group',
           ...
           u'Klebsiella oxytoca', u'Yersinia enterocolitica',
           u'Serratia plymuthica', u'Listeria monocytogenes', u'Shigella boydii',
           u'Brevundimonas subvibrioides', u'Buchnera aphidicola',
           u'Variovorax paradoxus', u'Rhodobacter capsulatus',
           u'Lactobacillus salivarius'],
          dtype='object', name=u'species', length=104)




```python
#stacked bar graph using total species
sortindex = speciesratio.sum().sort_values(ascending=False).index
speciesratio_sorted = speciesratio.loc[:,sortindex] #Sort column order from the most abundant one to the least
indexorder = speciesratio_sorted.index
#speciesratio_sorted.drop('None',axis=1,inplace=True)

sns.set(style="white")
plt.figure(figsize=(12,16))

speciesratio_sorted.loc[indexorder].iloc[:,:10].plot(kind='bar', stacked=True,legend=False, ylim=(0,100), colormap='Paired')
plt.xticks(rotation=90, size=20)
plt.yticks(size=20)
plt.legend(fontsize=10, bbox_to_anchor=(1.6, 1.0))
#plt.savefig("./totalratio_bar_species_humanremoved_None_notconsidered_nocontrol.png", format='png',bbox_inches='tight')
#plt.savefig("./figure/totalratio_bar.svg", format='svg')
#plt.savefig("./figure/totalratio_bar.pdf", format='pdf')
```




    <matplotlib.legend.Legend at 0x7fbd9bf86cd0>




    <matplotlib.figure.Figure at 0x7fbd95eec190>



![png](output_64_2.png)



```python
level = 'species'

krakendf['count']=1
speciesdf = krakendf.groupby(['sample',level]).sum()   #You can change 'species' to any level that you want.

speciesdf_un = speciesdf['count'].unstack()

speciesdf_un.fillna(0, inplace=True)
speciesdf_un['tot'] = speciesdf_un.sum(axis=1)


for i in speciesdf_un.index:
    speciesdf_un.loc[i,:] = speciesdf_un.loc[i,:]/speciesdf_un.loc[i,'tot'] *100

speciesratio = speciesdf_un.drop('tot',axis=1)

#stacked bar graph using total species
sortindex = speciesratio.sum().sort_values(ascending=False).index
speciesratio_sorted = speciesratio.loc[:,sortindex] #Sort column order from the most abundant one to the least
indexorder = speciesratio_sorted.index
#speciesratio_sorted.drop('None',axis=1,inplace=True)

sns.set(style="white")
plt.figure(figsize=(12,16))

speciesratio_sorted.loc[indexorder].iloc[:,:10].plot(kind='bar', stacked=True,legend=False, ylim=(0,100), colormap='Paired')
plt.xticks(rotation=90, size=20)
plt.yticks(size=20)
plt.legend(fontsize=10, bbox_to_anchor=(1.6, 1.0))
#plt.savefig("./totalratio_bar_%s_humanremoved_None_notconsidered_nocontrol.png" %level, format='png',bbox_inches='tight')
#plt.savefig("./figure/totalratio_bar.svg", format='svg')
#plt.savefig("./figure/totalratio_bar.pdf", format='pdf')
```




    <matplotlib.legend.Legend at 0x7fbd96027d10>




    <matplotlib.figure.Figure at 0x7fbde0e26c50>



![png](output_65_2.png)



```python
speciesratio_sorted.index
```




    Index([u'niams34-covariscontrol', u'niams34-extractioncontrol',
           u'niams34-v3-ll', u'niams34-v3-lnl', u'niams34-v4-ll',
           u'niams34-v4-lnl', u'niams34-v5-ll', u'niams34-v5-lnl',
           u'niams34-v8-ll', u'niams34-v8-lnl'],
          dtype='object', name=u'sample')




```python
ll = speciesratio_sorted.loc[['niams34-v3-ll','niams34-v4-ll','niams34-v5-ll','niams34-v8-ll',
                              'niams34-covariscontrol','niams34-extractioncontrol'],]
ll.plot(kind='bar', stacked=True,legend=False, ylim=(0,100))
plt.xticks(rotation=90)
plt.legend(fontsize=10, bbox_to_anchor=(1.9, 1.0))
#plt.savefig("./figure/totalratio_bar_lesion.png", format='png', bbox_inches='tight')
#plt.savefig("./figure/totalratio_bar_lesion.svg", format='svg', bbox_inches='tight')
```




    <matplotlib.legend.Legend at 0x7f718aade290>




![png](output_67_1.png)



```python
nl = speciesratio_sorted.loc[['niams34-v3-lnl','niams34-v4-lnl','niams34-v5-lnl','niams34-v8-lnl',
                             'niams34-covariscontrol','niams34-extractioncontrol'],]
nl.plot(kind='bar', stacked=True,legend=False, ylim=(0,100))
plt.xticks(rotation=90)
plt.legend(fontsize=10, bbox_to_anchor=(1.9, 1.0))
#plt.savefig("./figure/totalratio_bar_nonlesion.png", format='png', bbox_inches='tight')
#plt.savefig("./figure/totalratio_bar_nonlesion.svg", format='svg', bbox_inches='tight')
```




    <matplotlib.legend.Legend at 0x7f7181a33290>




![png](output_68_1.png)



```python
#stacked bar graph using top5 species
species_top5 = speciesratio[['Staphylococcus aureus','Staphylococcus epidermidis','Propionibacteriaceae','Alteromonas macleodii']]
species_top5.columns = ['Staphylococcus aureus','Staphylococcus epidermidis','Propionibacterium acnes','Alteromonas macleodii']
species_top5.plot(kind='bar', stacked=True, ylim=(0,100))
plt.legend(fontsize=10, bbox_to_anchor=(1.5, 1.0))
plt.xticks(rotation=90)
#plt.savefig("./figure/top5ratio_bar.png", format='png')
#plt.savefig("./figure/top5ratio_bar.svg", format='svg')
```




    (array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]), <a list of 10 Text xticklabel objects>)




![png](output_69_1.png)



```python
ll5 = species_top5.iloc[[2,4,6,8],]
ll5.plot(kind='bar', stacked=True,legend=True, ylim=(0,100))
plt.xticks(rotation=0)
plt.legend(fontsize=10, bbox_to_anchor=(1.5, 1.0))
#plt.savefig("./figure/totalratio_bar_lesion5.png", format='png', bbox_inches='tight')
#plt.savefig("./figure/totalratio_bar_lesion5.svg", format='svg', bbox_inches='tight')
```




    <matplotlib.legend.Legend at 0x7f7174ca45d0>




![png](output_70_1.png)



```python
nl5 = species_top5.iloc[[3,5,7,9],]
nl5.plot(kind='bar', stacked=True,legend=True, ylim=(0,100))
plt.xticks(rotation=0)
plt.legend(fontsize=10, bbox_to_anchor=(1.5, 1.0))
#plt.savefig("./figure/totalratio_bar_nonlesion5.png", format='png', bbox_inches='tight')
#plt.savefig("./figure/totalratio_bar_nonlesion5.svg", format='svg', bbox_inches='tight')
```




    <matplotlib.legend.Legend at 0x7f7175142690>




![png](output_71_1.png)



```python
humanRNA['nonhuman']
```




    v3l     72.32
    v3nl    61.67
    v4l     95.25
    v4nl    20.36
    v5l     60.78
    v5nl    14.42
    v8l     15.59
    v8nl    18.19
    Name: nonhuman, dtype: float64




```python
#Normalize using nonhuman reads ratio to see absolute level (This can be done using total number of cells. If Teru has this data, we have to use it.)

species_ratio_normalized_by_number = speciesratio
for i in humanRNA.index:
    species_ratio_normalized_by_number.loc[i,:] *= humanRNA.loc[i,'nonhuman']

    
sortindex = speciesratio.sum().sort_values(ascending=False).index
speciesrationorm_sorted = species_ratio_normalized_by_number.loc[:,sortindex]

speciesrationorm_sorted.plot(kind='bar', stacked=True,legend=False)
plt.xticks(rotation=0)
```




    (array([0, 1, 2, 3, 4, 5, 6, 7]), <a list of 8 Text xticklabel objects>)




![png](output_73_1.png)


# Diversity analysis

## alpha diversity


```python
from skbio.diversity import alpha_diversity
```


```python
ratio_list = speciesdf['count'].unstack().fillna(0).astype('int')

ratio_list = ratio_list.values.tolist()
```


```python
otucounts = alpha_diversity('observed_otus',ratio_list,ids=speciesratio.index )
otucounts
```




    sample
    v3l     112
    v3nl    117
    v4l      39
    v4nl    191
    v5l      95
    v5nl    149
    v8l     103
    v8nl    141
    dtype: int64




```python
otucounts.plot(kind='bar')
plt.ylabel("OTU counts")
```




    <matplotlib.text.Text at 0x7f87428167d0>




![png](output_79_1.png)


## beta diversity


```python
def pw_distances(counts, ids=None, metric="braycurtis"):
    """Compute distances between all pairs of columns in a counts matrix
    Parameters
    ----------
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of observations in a given sample.
    ids : iterable of strs, optional
        Identifiers for each sample in ``counts``.
    metric : str, optional
        The name of the pairwise distance function to use when generating
        pairwise distances. See the scipy ``pdist`` docs, linked under *See
        Also*, for available metrics.
    Returns
    -------
    skbio.DistanceMatrix
        Distances between all pairs of samples (i.e., rows). The number of
        row and columns will be equal to the number of rows in ``counts``.
    Raises
    ------
    ValueError
        If ``len(ids) != len(counts)``.
    See Also
    --------
    scipy.spatial.distance.pdist
    pw_distances_from_table
    """
    num_samples = len(counts)
    if ids is not None and num_samples != len(ids):
        raise ValueError(
            "Number of rows in counts must be equal to number of provided "
            "ids.")

    distances = pdist(counts, metric)
    return DistanceMatrix(
        squareform(distances, force='tomatrix', checks=False), ids)
```


```python
import numpy as np

from skbio import DistanceMatrix
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

bc_dm = pw_distances(ratio_list,ids=speciesratio.index, metric="braycurtis")

bc_dm
```




![svg](output_82_0.svg)




```python
from skbio.stats.ordination import pcoa
bc_pc = pcoa(bc_dm)
```

    /home/jk/anaconda2/lib/python2.7/site-packages/skbio/stats/ordination/_principal_coordinate_analysis.py:102: RuntimeWarning: The result contains negative eigenvalues. Please compare their magnitude with the magnitude of some of the largest positive eigenvalues. If the negative ones are smaller, it's probably safe to ignore them, but if they are large in magnitude, the results won't be useful. See the Notes section for more details. The smallest eigenvalue is -0.000222798712441 and the largest is 1.30826505826.
      RuntimeWarning



```python
sample_md={
    'v3l':{'name': 'v3l', 'eczema': 'lesional'},
    'v3nl':{'name': 'v3nl', 'eczema': 'nonlesional'},
    'v4l':{'name': 'v4l', 'eczema': 'lesional'},
    'v4nl':{'name': 'v4nl', 'eczema': 'nonlesional'},
    'v5l':{'name': 'v5l', 'eczema': 'lesional'},
    'v5nl':{'name': 'v5nl', 'eczema': 'nonlesional'},
    'v8l':{'name': 'v8l', 'eczema': 'lesional'},
    'v8nl':{'name': 'v8nl', 'eczema': 'nonlesional'},
}
sample_md = pd.DataFrame.from_dict(sample_md, orient='index')
sample_md
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>name</th>
      <th>eczema</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>v3l</th>
      <td>v3l</td>
      <td>lesional</td>
    </tr>
    <tr>
      <th>v3nl</th>
      <td>v3nl</td>
      <td>nonlesional</td>
    </tr>
    <tr>
      <th>v4l</th>
      <td>v4l</td>
      <td>lesional</td>
    </tr>
    <tr>
      <th>v4nl</th>
      <td>v4nl</td>
      <td>nonlesional</td>
    </tr>
    <tr>
      <th>v5l</th>
      <td>v5l</td>
      <td>lesional</td>
    </tr>
    <tr>
      <th>v5nl</th>
      <td>v5nl</td>
      <td>nonlesional</td>
    </tr>
    <tr>
      <th>v8l</th>
      <td>v8l</td>
      <td>lesional</td>
    </tr>
    <tr>
      <th>v8nl</th>
      <td>v8nl</td>
      <td>nonlesional</td>
    </tr>
  </tbody>
</table>
</div>




```python
#pcoa plot
fig = bc_pc.plot(sample_md, 'name',
                axis_labels=('PC 1', 'PC 2', 'PC 3'),
                title='PcoA', cmap='Set1', s=50)
```


![png](output_85_0.png)



```python

```

# Species specific analysis

## S. aureus


```python
infile = !ls -1 *nohuman.norrna.fastq
infile
```


```python
!mkdir saureus
```


```python
#bowtie2 mapping to Staphylococcus aureus subsp. aureus NCTC 8325 
for i in range(len(infile)):
    saudbloc = '/home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau'
    targetfile = infile[i]
    samfile = infile[i] + '.sau.sam'
    logfile = infile[i] + '.sau.bowtie2log'
    
    command_line = 'bowtie2 -x %s -U %s -S ./saureus/%s 2> ./saureus/%s' %(saudbloc,targetfile,samfile,logfile)
    print "!" + command_line

```

    `bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau -U v3l_trimmed.nohuman.norrna.fastq -S ./saureus/v3l_trimmed.nohuman.norrna.fastq.sau.sam 2> ./saureus/v3l_trimmed.nohuman.norrna.fastq.sau.bowtie2log`;
    `bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau -U v4l_trimmed.nohuman.norrna.fastq -S ./saureus/v4l_trimmed.nohuman.norrna.fastq.sau.sam 2> ./saureus/v4l_trimmed.nohuman.norrna.fastq.sau.bowtie2log`;
    `bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau -U v5l_trimmed.nohuman.norrna.fastq -S ./saureus/v5l_trimmed.nohuman.norrna.fastq.sau.sam 2> ./saureus/v5l_trimmed.nohuman.norrna.fastq.sau.bowtie2log`;
    `bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau -U v8l_trimmed.nohuman.norrna.fastq -S ./saureus/v8l_trimmed.nohuman.norrna.fastq.sau.sam 2> ./saureus/v8l_trimmed.nohuman.norrna.fastq.sau.bowtie2log`;



```python
!bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau -U v3l_trimmed.nohuman.norrna.fastq -S ./saureus/v3l_trimmed.nohuman.norrna.fastq.sau.sam 2> ./saureus/v3l_trimmed.nohuman.norrna.fastq.sau.bowtie2log
!bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau -U v4l_trimmed.nohuman.norrna.fastq -S ./saureus/v4l_trimmed.nohuman.norrna.fastq.sau.sam 2> ./saureus/v4l_trimmed.nohuman.norrna.fastq.sau.bowtie2log
!bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau -U v5l_trimmed.nohuman.norrna.fastq -S ./saureus/v5l_trimmed.nohuman.norrna.fastq.sau.sam 2> ./saureus/v5l_trimmed.nohuman.norrna.fastq.sau.bowtie2log
!bowtie2 -x /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau -U v8l_trimmed.nohuman.norrna.fastq -S ./saureus/v8l_trimmed.nohuman.norrna.fastq.sau.sam 2> ./saureus/v8l_trimmed.nohuman.norrna.fastq.sau.bowtie2log
```


```python
#saureus mapping ratio
!head ./saureus/*bowtie2log
```


```python
#read count using Staphylococcus aureus subsp. aureus NCTC 8325 mapping file
for i in range(len(infile)):
    sausaf = '/home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau.saf'
    targetfile = infile[i]
    samfile = infile[i] + '.sau.sam'
    outfile = infile[i][0:3] + '.sau_readcounts.txt'
    
    command_line = 'featureCounts -T 4 -t GeneID -F SAF -g locus_tag -o ./saureus/%s -a %s ./saureus/%s' %(outfile,sausaf,samfile)
    print "!" + command_line

```

    !featureCounts -T 4 -t GeneID -F SAF -g locus_tag -o ./saureus/v3l.sau_readcounts.txt -a /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau.saf ./saureus/v3l_trimmed.nohuman.norrna.fastq.sau.sam
    !featureCounts -T 4 -t GeneID -F SAF -g locus_tag -o ./saureus/v4l.sau_readcounts.txt -a /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau.saf ./saureus/v4l_trimmed.nohuman.norrna.fastq.sau.sam
    !featureCounts -T 4 -t GeneID -F SAF -g locus_tag -o ./saureus/v5l.sau_readcounts.txt -a /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau.saf ./saureus/v5l_trimmed.nohuman.norrna.fastq.sau.sam
    !featureCounts -T 4 -t GeneID -F SAF -g locus_tag -o ./saureus/v8l.sau_readcounts.txt -a /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau.saf ./saureus/v8l_trimmed.nohuman.norrna.fastq.sau.sam



```python
!featureCounts -T 4 -t GeneID -F SAF -g locus_tag -o ./saureus/v3l.sau_readcounts.txt -a /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau.saf ./saureus/v3l_trimmed.nohuman.norrna.fastq.sau.sam
!featureCounts -T 4 -t GeneID -F SAF -g locus_tag -o ./saureus/v4l.sau_readcounts.txt -a /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau.saf ./saureus/v4l_trimmed.nohuman.norrna.fastq.sau.sam
!featureCounts -T 4 -t GeneID -F SAF -g locus_tag -o ./saureus/v5l.sau_readcounts.txt -a /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau.saf ./saureus/v5l_trimmed.nohuman.norrna.fastq.sau.sam
!featureCounts -T 4 -t GeneID -F SAF -g locus_tag -o ./saureus/v8l.sau_readcounts.txt -a /home/jk/Dropbox/UCSD/Research/Data/genome_reference/saureus/sau.saf ./saureus/v8l_trimmed.nohuman.norrna.fastq.sau.sam
```


```python
for i in range(len(infile)):
    readcount = infile[i][0:3] + '.sau_readcounts.txt'
    output = readcount + '2'
    
    command_line = 'tail -n +3 ./saureus/%s | cut -f 1,7 >./saureus/%s' %(readcount,output)
    subprocess.Popen(command_line,shell=True)
```


```python
!head ./saureus/*2
```


```python
count = {}
for i in range(len(infile)):
    sample = infile[i][0:3]
    countfile = infile[i][0:3] + '.sau_readcounts.txt2'
    count[sample] = pd.read_table("./saureus/%s" %countfile, header=None, index_col=0)[1]

readtable = pd.DataFrame(count)
readtable.to_csv('./saureus/niams37_saureus.txt', sep='\t')
readtable.head()
```


```python
#Use readcount table in R for DESeq2 analysis
```


```python

```


```python
#S. aureus mapping ratio
df = pd.DataFrame({'v3l':88.27, 'v3nl':87.33, 'v4l':64.68, 'v4nl':6.41, 'v5l':89.44, 'v5nl':6.02, 'v8l':44.16, 'v8nl':1.70},
            index=['saureus'])
df = df.transpose()
df['others']=100-df['saureus']
df
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>saureus</th>
      <th>others</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>v3l</th>
      <td>88.27</td>
      <td>11.73</td>
    </tr>
    <tr>
      <th>v3nl</th>
      <td>87.33</td>
      <td>12.67</td>
    </tr>
    <tr>
      <th>v4l</th>
      <td>64.68</td>
      <td>35.32</td>
    </tr>
    <tr>
      <th>v4nl</th>
      <td>6.41</td>
      <td>93.59</td>
    </tr>
    <tr>
      <th>v5l</th>
      <td>89.44</td>
      <td>10.56</td>
    </tr>
    <tr>
      <th>v5nl</th>
      <td>6.02</td>
      <td>93.98</td>
    </tr>
    <tr>
      <th>v8l</th>
      <td>44.16</td>
      <td>55.84</td>
    </tr>
    <tr>
      <th>v8nl</th>
      <td>1.70</td>
      <td>98.30</td>
    </tr>
  </tbody>
</table>
</div>




```python
df.plot(kind='bar',stacked='True',color=['#4682b4','Grey'])
plt.xticks(rotation=0)
plt.legend(fontsize=10, bbox_to_anchor=(1.3, 1.0))
plt.savefig("./figure/saureus_mapping.png", format='png')
plt.savefig("./figure/saureus_mapping.svg", format='svg')
```


![png](output_102_0.png)


## S. epidermidis


```python
#bowtie2 mapping to Staphylococcus epidermidis ATCC1228
!mkdir sepi
!bowtie2 -x ../../../genome_reference/sepi/sepi -U ./AMT-V3-L-R1.humanremoved.norrna.fastq -S ./sepi/AMT-V3-L-R1.humanremoved.norrna.fastq.sepi.sam 2>./sepi/AMT-V3-L-R1.humanremoved.norrna.fastq.sepi.bowtie2
!bowtie2 -x ../../../genome_reference/sepi/sepi -U ./AMT-V3-NL-R1.humanremoved.norrna.fastq -S ./sepi/AMT-V3-NL-R1.humanremoved.norrna.fastq.sepi.sam 2>./sepi/AMT-V3-NL-R1.humanremoved.norrna.fastq.sepi.bowtie2
!bowtie2 -x ../../../genome_reference/sepi/sepi -U ./AMT-V4-L-R1.humanremoved.norrna.fastq -S ./sepi/AMT-V4-L-R1.humanremoved.norrna.fastq.sepi.sam 2>./sepi/AMT-V4-L-R1.humanremoved.norrna.fastq.sepi.bowtie2
!bowtie2 -x ../../../genome_reference/sepi/sepi -U ./AMT-V4-NL-R1.humanremoved.norrna.fastq -S ./sepi/AMT-V4-NL-R1.humanremoved.norrna.fastq.sepi.sam 2>./sepi/AMT-V4-NL-R1.humanremoved.norrna.fastq.sepi.bowtie2
!bowtie2 -x ../../../genome_reference/sepi/sepi -U ./AMT-V5-L-R1.humanremoved.norrna.fastq -S ./sepi/AMT-V5-L-R1.humanremoved.norrna.fastq.sepi.sam 2>./sepi/AMT-V5-L-R1.humanremoved.norrna.fastq.sepi.bowtie2
!bowtie2 -x ../../../genome_reference/sepi/sepi -U ./AMT-V5-NL-R1.humanremoved.norrna.fastq -S ./sepi/AMT-V5-NL-R1.humanremoved.norrna.fastq.sepi.sam 2>./sepi/AMT-V5-NL-R1.humanremoved.norrna.fastq.sepi.bowtie2
!bowtie2 -x ../../../genome_reference/sepi/sepi -U ./AMT-V8-L-R1.humanremoved.norrna.fastq -S ./sepi/AMT-V8-L-R1.humanremoved.norrna.fastq.sepi.sam 2>./sepi/AMT-V8-L-R1.humanremoved.norrna.fastq.sepi.bowtie2
!bowtie2 -x ../../../genome_reference/sepi/sepi -U ./AMT-V8-NL-R1.humanremoved.norrna.fastq -S ./sepi/AMT-V8-NL-R1.humanremoved.norrna.fastq.sepi.sam 2>./sepi/AMT-V8-NL-R1.humanremoved.norrna.fastq.sepi.bowtie2
```


```python
#sam to bam
!samtools view -bS ./sepi/AMT-V3-L-R1.humanremoved.norrna.fastq.sepi.sam >./sepi/AMT-V3-L-R1.humanremoved.norrna.fastq.sepi.sam.bam
!samtools view -bS ./sepi/AMT-V3-NL-R1.humanremoved.norrna.fastq.sepi.sam >./sepi/AMT-V3-NL-R1.humanremoved.norrna.fastq.sepi.sam.bam
!samtools view -bS ./sepi/AMT-V4-L-R1.humanremoved.norrna.fastq.sepi.sam >./sepi/AMT-V4-L-R1.humanremoved.norrna.fastq.sepi.sam.bam
!samtools view -bS ./sepi/AMT-V4-NL-R1.humanremoved.norrna.fastq.sepi.sam >./sepi/AMT-V4-NL-R1.humanremoved.norrna.fastq.sepi.sam.bam
!samtools view -bS ./sepi/AMT-V5-L-R1.humanremoved.norrna.fastq.sepi.sam >./sepi/AMT-V5-L-R1.humanremoved.norrna.fastq.sepi.sam.bam
!samtools view -bS ./sepi/AMT-V5-NL-R1.humanremoved.norrna.fastq.sepi.sam >./sepi/AMT-V5-NL-R1.humanremoved.norrna.fastq.sepi.sam.bam
!samtools view -bS ./sepi/AMT-V8-L-R1.humanremoved.norrna.fastq.sepi.sam >./sepi/AMT-V8-L-R1.humanremoved.norrna.fastq.sepi.sam.bam
!samtools view -bS ./sepi/AMT-V8-NL-R1.humanremoved.norrna.fastq.sepi.sam >./sepi/AMT-V8-NL-R1.humanremoved.norrna.fastq.sepi.sam.bam
```

    [samopen] SAM header is present: 1 sequences.
    [samopen] SAM header is present: 1 sequences.
    [samopen] SAM header is present: 1 sequences.
    [samopen] SAM header is present: 1 sequences.
    [samopen] SAM header is present: 1 sequences.
    [samopen] SAM header is present: 1 sequences.
    [samopen] SAM header is present: 1 sequences.
    [samopen] SAM header is present: 1 sequences.



```python
#sepi mapping ratio
!head ./sepi/*bowtie2
```

    ==> ./sepi/AMT-V3-L-R1.humanremoved.norrna.fastq.sepi.bowtie2 <==
    149486 reads; of these:
      149486 (100.00%) were unpaired; of these:
        122147 (81.71%) aligned 0 times
        26771 (17.91%) aligned exactly 1 time
        568 (0.38%) aligned >1 times
    18.29% overall alignment rate
    
    ==> ./sepi/AMT-V3-NL-R1.humanremoved.norrna.fastq.sepi.bowtie2 <==
    135710 reads; of these:
      135710 (100.00%) were unpaired; of these:
        106697 (78.62%) aligned 0 times
        28500 (21.00%) aligned exactly 1 time
        513 (0.38%) aligned >1 times
    21.38% overall alignment rate
    
    ==> ./sepi/AMT-V4-L-R1.humanremoved.norrna.fastq.sepi.bowtie2 <==
    137449 reads; of these:
      137449 (100.00%) were unpaired; of these:
        73599 (53.55%) aligned 0 times
        63296 (46.05%) aligned exactly 1 time
        554 (0.40%) aligned >1 times
    46.45% overall alignment rate
    
    ==> ./sepi/AMT-V4-NL-R1.humanremoved.norrna.fastq.sepi.bowtie2 <==
    35173 reads; of these:
      35173 (100.00%) were unpaired; of these:
        25553 (72.65%) aligned 0 times
        9505 (27.02%) aligned exactly 1 time
        115 (0.33%) aligned >1 times
    27.35% overall alignment rate
    
    ==> ./sepi/AMT-V5-L-R1.humanremoved.norrna.fastq.sepi.bowtie2 <==
    104473 reads; of these:
      104473 (100.00%) were unpaired; of these:
        90603 (86.72%) aligned 0 times
        13753 (13.16%) aligned exactly 1 time
        117 (0.11%) aligned >1 times
    13.28% overall alignment rate
    
    ==> ./sepi/AMT-V5-NL-R1.humanremoved.norrna.fastq.sepi.bowtie2 <==
    21407 reads; of these:
      21407 (100.00%) were unpaired; of these:
        16706 (78.04%) aligned 0 times
        4669 (21.81%) aligned exactly 1 time
        32 (0.15%) aligned >1 times
    21.96% overall alignment rate
    
    ==> ./sepi/AMT-V8-L-R1.humanremoved.norrna.fastq.sepi.bowtie2 <==
    22724 reads; of these:
      22724 (100.00%) were unpaired; of these:
        20583 (90.58%) aligned 0 times
        2114 (9.30%) aligned exactly 1 time
        27 (0.12%) aligned >1 times
    9.42% overall alignment rate
    
    ==> ./sepi/AMT-V8-NL-R1.humanremoved.norrna.fastq.sepi.bowtie2 <==
    12739 reads; of these:
      12739 (100.00%) were unpaired; of these:
        12504 (98.16%) aligned 0 times
        232 (1.82%) aligned exactly 1 time
        3 (0.02%) aligned >1 times
    1.84% overall alignment rate



```python
#S. epidermidis mapping ratio
df = pd.DataFrame({'v3l':18.29, 'v3nl':21.38, 'v4l':46.45, 'v4nl':21.35, 'v5l':13.28, 'v5nl':21.96, 'v8l':9.42, 'v8nl':1.84},
            index=['sepi'])
df = df.transpose()
df['others']=100-df['sepi']
df
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sepi</th>
      <th>others</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>v3l</th>
      <td>18.29</td>
      <td>81.71</td>
    </tr>
    <tr>
      <th>v3nl</th>
      <td>21.38</td>
      <td>78.62</td>
    </tr>
    <tr>
      <th>v4l</th>
      <td>46.45</td>
      <td>53.55</td>
    </tr>
    <tr>
      <th>v4nl</th>
      <td>21.35</td>
      <td>78.65</td>
    </tr>
    <tr>
      <th>v5l</th>
      <td>13.28</td>
      <td>86.72</td>
    </tr>
    <tr>
      <th>v5nl</th>
      <td>21.96</td>
      <td>78.04</td>
    </tr>
    <tr>
      <th>v8l</th>
      <td>9.42</td>
      <td>90.58</td>
    </tr>
    <tr>
      <th>v8nl</th>
      <td>1.84</td>
      <td>98.16</td>
    </tr>
  </tbody>
</table>
</div>




```python
df.plot(kind='bar',stacked='True',color=['#2e8b57','Grey'])
plt.xticks(rotation=0)
plt.legend(fontsize=10, bbox_to_anchor=(1.3, 1.0))
plt.savefig("./figure/sepi_mapping.png", format='png')
plt.savefig("./figure/sepi_mapping.svg", format='svg')
```


![png](output_108_0.png)



```python

```
