[![Build Status](https://travis-ci.org/pratas/ac.svg?branch=master)](https://travis-ci.org/pratas/ac)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
<p align="center"><img src="imgs/logo.png"
alt="AC" width="200" border="0" /></p>
<p align="center"><b>AC: a lossless compression tool for amino acid sequences.</b></p></br>

<p align="justify">
<b>AC</b> is a new lossless compressor to compress efficiently amino acid sequences (proteins). It uses a cooperation between multiple context and substitutional tolerant context models. The cooperation between models is balanced with weights that benefit the models with better performance according to a forgetting function specific for each model.

## 1. INSTALLATION ##

Downloading and installing AC:
<pre>
git clone https://github.com/pratas/ac.git
cd ac/src/
cmake .
make
</pre>

Cmake is needed for the installation (http://www.cmake.org/). You can download it directly from http://www.cmake.org/cmake/resources/software.html or use an appropriate packet manager, such as:
<pre>
sudo apt-get install cmake
</pre>
An alternative to cmake, but limited to Linux, can be set using the following instructions:
<pre>
cp Makefile.linux Makefile
make
</pre>

## 2. USAGE ##

To see the possible options of AC type
<pre>
./AC
</pre>
or
<pre>
./AC -h
</pre>
These will print the following options:
<pre>
<p>
Usage: AC [OPTION]... -r [FILE]  [FILE]:[...]                          
Compression of amino acid sequences.                                   
                                                                       
Non-mandatory arguments:                                               
                                                                       
  -h                     give this help,                               
  -s                     show AC compression levels,                   
  -v                     verbose mode (more information),              
  -V                     display version number,                       
  -f                     force overwrite of output,                    
  -l &#60level&#62             level of compression [1;7] (lazy -tm setup),  
  -t &#60threshold&#62         threshold frequency to discard from alphabet,
  -e                     it creates a file with the extension ".iae" 
                         with the respective information content.      
                                                                       
  -rm &#60c&#62:&#60d&#62:&#60g&#62/&#60m&#62:&#60e&#62:&#60a&#62  reference model (-rm 1:10:0.9/0:0:0),   
  -rm &#60c&#62:&#60d&#62:&#60g&#62/&#60m&#62:&#60e&#62:&#60a&#62  reference model (-rm 5:90:0.9/1:50:0.8),
  ...                                                                  
  -tm &#60c&#62:&#60d&#62:&#60g&#62/&#60m&#62:&#60e&#62:&#60a&#62  target model (-tm 1:1:0.8/0:0:0),       
  -tm &#60c&#62:&#60d&#62:&#60g&#62/&#60m&#62:&#60e&#62:&#60a&#62  target model (-tm 7:100:0.9/2:10:0.85),  
  ...                                                                  
                         target and reference templates use &#60c&#62 for    
                         context-order size, &#60d&#62 for alpha (1/&#60d&#62), &#60g&#62
                         for gamma (decayment forgetting factor) [0;1),
                         &#60m&#62 to the maximum sets the allowed mutations,
                         on the context without being discarded (for   
                         deep contexts), under the estimator &#60e&#62, using
                         &#60a&#62 for gamma (decayment forgetting factor)   
                         [0;1) (tolerant model),                       
                                                                       
  -r &#60FILE&#62              reference file ("-rm" are loaded here),     
                                                                       
Mandatory arguments:                                                   
                                                                       
  &#60FILE&#62:&#60...&#62:&#60...&#62     file to compress (last argument). For more    
                         files use splitting ":" characters.         
                                                                       
Example:                                                               
                                                                       
  [Compress]   ./AC -v -tm 1:1:0.8/0:0:0 -tm 5:20:0.9/3:20:0.9 seq.txt 
  [Decompress] ./AD -v seq.txt.co      

Report bugs to &#60{pratas,seyedmorteza,ap}@ua.pt&#62.                            
</pre>

## 3. EXAMPLE ##

After AC intallation, run the following:
<pre>
wget http://sweet.ua.pt/pratas/datasets/AminoAcidsCorpus.zip
unzip AminoAcidsCorpus.zip
cp AminoAcidsCorpus/HI .
./AC -v -l 2 HI
./AD -v HI.co
cmp HI HI.de
</pre>
It will download nine amino acid sequences and compress and decompress one of the smallest (HI). Finally, it compares if the uncompressed sequence is equal to the original.

## 4. CITATION ##

On using this tool/method, please, cite:

Pratas, D., Hosseini, M. and Pinho, A.J., 2018, May. Compression of Amino Acid Sequences. In International Conference on Practical Applications of Computational Biology &amp; Bioinformatics (pp. 105-113). Springer, Cham.

## 5. ISSUES ##

For any issue let us know at [issues link](https://github.com/pratas/ac/issues).

## 6. LICENSE ##

GPL v3.

For more information:
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>

