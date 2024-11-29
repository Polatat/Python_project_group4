# FilterFASTQ

## Description

FilterFASTQ is a Python program designed to filter and analyze FASTQ data effectively. It provides a statistical visualization to differentiate between the original and the filtered file by using a pie chart and generating  CSV files. In addition, it can identified by a group of barcodes or barcode sequences as the users desire. At the end of the program, it generates a new FASTQ file that passes the biological statistic criteria.


## Features


**Quality score Filtering:** Eliminate sequences that do not meet a minimum average quality score. The default is 20.

**Sequence Length Filtering:** This filter Eliminates sequences that do not meet the minimum sequence length, which is 50 by default.

**GC Content Filtering:** Eliminate sequences that have a GC percentage above or below the appropriate range. The default range is between 30 % and 60%.

**Group Barcode Analysis:** Visualization of barcode group distribution via pie chart.

**Writing statistical CSV files:** Generate statistical files for original and filtered data.

**Writing a new FASTQ file:** Generate a new FASTQ file that contains passed criteria sequences.

## Installation

### Prerequistes
- **Python 3.6 or higher**
- **pip** (Python package installer)
- **Git** (to clone the repository)

### Steps

1. **Clone the Repository**

    ```bash
    git clone https://github.com/Polatat/Python_project_group4
    cd Python_project_group4
    ```

2. **Set Up a Virtual Environment**

    It's recommended to use a virtual environment to manage dependencies.

    ```bash
    python3 -m venv biopython_env
    source biopython_env/bin/activate
    ```

    *On Windows:*

    ```batch
    python -m venv biopython_env
    biopython_env\Scripts\activate
    ```

3. **Install Dependencies**

    ```bash
    pip install -r requirements.txt
    ```
* if you do not have the specific dependencies on your machine please using
  
    ```bash
    pip install [specific package name]
    ```








