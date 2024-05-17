# CINDI3

## Overview
The CINDI3 repository contains all necessary code for processing Slant Column Densities (SCD) and for converting native `.txt` files to `.ascii` format. This repository is designed to integrate seamlessly with the existing Blick source code.

## Installation Instructions

### Prerequisites
Ensure that you have the following requirements installed:
- Git (for cloning the repository)
- Python 2.7 (required version depends on Blick source code compatibility)

### Cloning the Repository
To ensure full compatibility and functionality, clone the CINDI3 repository directly into the `C:/Blick/src` directory on your local machine. Follow these steps:
1. Open your command prompt or terminal.
2. Navigate to the Blick source directory:
   ```bash
   cd C:/Blick/src
3. Please clone the CINDI3 repository:
   ```bash
   git clone git@github.com:LuftBlick/CINDI3.git

## Running the CINDI3 code
In order to process your L0 files from the CINDI3 campaign, please place them as usually into the C:/Blick/data/L0 folder. 
1. Please create an output folder for the CINDI ASCII files named "CINDI3" into C:/Blick/data.
2. Jump into C:/Blick/src/CINDI3. In order to process L0 to L2Fit data and convert the SCD's finally into the needed CINDI ASCII format run the main function with the following command.
    ```bash
    C:\Python27\python.exe main_makeCompData.py

## What happens
The CINDI3 processing routine creates the L2Fit files and converts them into ASCII CINDI format. Here, in the folder C:/Blick/data/CINDI3 for every reference type a separate folder is created, where only the allowed processing types for each reference type are converted (e.g. Ref:PROFILE,SKY,SUN). These informations can be found in the CONFIG file "processCompDataInput.txt".
