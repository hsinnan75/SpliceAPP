# SpliceAPP
Splice-Alternative Profile Predictor (SpliceAPP) 

===================

Web-service: https://bc.imb.sinica.edu.tw/WebApps/SpliceAPP/

Standalone : https://github.com/hsinnan75/SpliceAPP

# Introduction

For 3'ss (splice site) prediction, variants loacated between -78 and -4 of the 3' end of the intron can be processed.

For 5'ss prediction, the predictive models accept variants from -3 to +30 of the 5' end of the intron.

Only hg38 annotation is accepted for the prediction.

# Installation

  ```
  $ git clone https://github.com/hsinnan75/SpliceAPP.git
  ```
to download the package of SpliceAPP.

# Dependencies

To run SpliceAPP, it requires docker installed and download the prediction models (16GB) at https://bc.imb.sinica.edu.tw/SpliceAPP/data.tgz

Please decompress the file in the SpliceAPP. It will create a folder of data.

# Compiling

To compile SpliceAPP, please change to the folder of SpliceAPP/src and just type 'make' to compile SpliceAPP. If the compilation or the program fails, please contact me (arith@gate.sinica.edu.tw), Thanks.

# Usage

Please type 

  ```
  $ src/SpliceAPP [vcf_file]
  ```
to run the program

