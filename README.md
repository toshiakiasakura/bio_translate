# bio_translate

bio_translate is a Python library for translating articles and tagging texts into Japanese, and output as markdown or html. 

The program workflow is, 
1. Convert doi (digital object identifier) into PMCID (PubMed Central ID). 
2. Fetch full texts (sometimes abstract only) in xml format with PMCID.
3. Translate texts with python googletrans API (unofficial). 
4. Format texts with markdown notation. 
5. Convert markdown notation into html notation. 

## Sample 
`python main.py --csv sample.csv` will produce sample outputs. 
See files in markdon/html directory for outputs. 


"status temp" in sample.csv indicates results of the process. 
- "1" is success flag.
- "no journal" indicates PMC do not have articles of that jornal in database. 
- "abstract only" indicates output only include abstract only.

## Usage 
`python main.py -h` will show this help.

``` 
usage: main.py [-h] [--doi DOI] [--name NAME] [--csv CSV] [--open-html] [--access-doi] [--json2html]
               [--check CHECK]

Tranlate article from doi into markdown/html.

optional arguments:
  -h, --help            show this help message and exit  
  --doi DOI, -d DOI     One doi to markdon. Should be set "--name/-n" option if this is selected.  
  --name NAME, -n NAME  Name for each file. Name should be without extension. Should be set with "--  
                        name/-n" option.  
  --csv CSV             Process doi from csv. csv should include "doi" and "name" header.  
  --open-html, -o       Open tab when html are produced.  
  --access-doi, -a      Open target doi.org page when doi is not correctly converted.  
  --json2html, -j       Convert process starts from json file. If html file does not exist, raise  
                        error.  
  --check CHECK, -c CHECK  
                        List up journals used in pubmed central database. Journal data are in  
                        jlist.csv.  
```

### options : --csv 
--csv option takes a csv file as argument.  
csv file must include the following columns: "name", "doi", "status".
- "name" column is used for naming file for xml, json, markdown, and html 
so that file extension is not needed.
- "doi" column is doi. This will be converted into PMCID.
- "status" column is used for controlling which row will be used for converting. 
In this module, if status column has nan value, the rows are used. 
If some values exist, that row is not used. 





