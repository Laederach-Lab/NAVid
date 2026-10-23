# NAVid: Identification of NAGNAG-Affecting Variants
Identify genomic variants that create, alter, or destroy a NAGNAG in the human genome and understand the proteomic consequences.

## Installation and Usage
To access the app via the web, go to https://nagnag.bio.unc.edu.

To run locally:
1. Clone the repository
2. Install the required Python packages
3. Download and parse reference genome GRCh38/hg38
4. Launch app

```
git clone https://github.com/Laederach-Lab/NAVid.git
cd NAVid
pip install -r requirements.txt
python parse-hg38.py
shiny run app.py --reload --launch-browser
```

__NOTE:__ The local version does not include any information on HGMD variants. Summary statistics pertaining to HGMD variants can be found on the webserver.

### Data Updates
To update the source data:
1. Re-run the analyses from https://github.com/hallep/nagnag-3ss
2. Run `update-data.py </path/to/nagnag-3ss/> </path/to/NAVid/>`

## Support
If you have questions or run into issues, contact us at hpearce@unc.edu or alain@unc.edu.

## Authors
This project was created by Halle Pearce of the Laederach Lab at the University of North Carolina at Chapel Hill.