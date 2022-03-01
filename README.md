# NetworkGraphMaker
this tool is a method to generate a gml file which can be the input data of Cytoscape.

# usage
To run this tool, you should prepare the result file of GTDB-Tk, CheckM, bwa (Hi-C reads to scaffolds, in sam format), readcount.
```
python gmlFileMaker.py (1)scaffoldList.txt (2)GTDB-Tk_result.tsv (3)CheckM_result.tsv (4)samfile.sam (5)readcount.tsv

scaffoldList.txt    list of scaffold name
GTDB-Tk_result.tsv  GTDB-Tk result file
CheckM_result.tsv   CheckM result file
samfile.sam         samfile
readcount.tsv       mapped reads count for each scaffold
```

## Input file format
・scaffoldList.txt
```
scaffold1
scaffold2
scaffold3
    ・
    ・
    ・
```
・readcount.tsv
```
scaffold1 3
scaffold2 10
scaffold3 5
    ・
    ・
    ・
```
