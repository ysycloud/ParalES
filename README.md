# ParalES
Enrichment Score(GSEA) compute Tools for Lincs data

## Tools

A brief description of the tools is given below. The tools implement the function of parallel compute ES in multi-core environment. The Matlab implementation of the tools is currently the most mature. I used the 1ktools(https://github.com/cmap/l1ktools) tools to parse the .gctx file which stored gene profile data defined by Lincs(CMap) based on HDF5 file format. However,the performance is not very satisfactory. The C implementation to achieve at least 30 times faster than Matlab. But 1ktools can not parse the .gctx file by C. So I used Matlab to parse the .gctx file、extract the gene profile sets and write to .txt file. C will read the file 、complete parallel computing ES and write out the result in binary file . There also is a Matlab Script to read the binary files and can be later analysis. I will update the tools as they become available.

### Matlab Tools: matlab/

#### Requirements:

1. Matlab R2009a and above

#### Setting the MATLAB path:
Enter the "pathtool" command, click "Add with Subfolders...", and select the directory ParalES/matlab.


#### Tools:
* [**ESquick.m**]: Compute single Enrichment Scores.
* [**ESScore.m**]: parallel compute Multiple pairs of ES.
* [**getSampleforMat.m**]: parse the .gctx file and extract the gene profile sets we need.
* [**ParalES.m**]：main entrance . Setting Parameters and complete parallel Computing tasks.
* [**getSample.m**]: parse the .gctx file and extract the gene profile sets we need for Writting.
* [**PreESforC.m**] : Setting Parameters 、 extract the gene profile sets and write to .txt file.
* [**importES.m**] : Read the Enrichment Scores Result Matrix computed by C.



### C Tools: c/

#### Tools:

* [**ParalES.c**] read the .txt file 、complete parallel computing ES and write out the result in binary file

#### Demo:
* [**compileParalES.sh**]: compile the C source code to Executable files.
* [**runParalESLinux.sh**]: run Executable files in Linux.



## The LINCS Dataset

The CMAP Cloud API offers programmatic access to annotations and perturbational signatures in [the LINCS L1000 dataset](http://lincscloud.org/) via a collection of HTTP-based RESTful web services. You can get the .gctx file stored gene profile data by the Website.
