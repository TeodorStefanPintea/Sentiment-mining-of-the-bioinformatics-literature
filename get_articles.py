'''
    This script is meant to download all the articles/journals that have been uploaded to pubmed in the past year.
    It uses biopython to acomplish the task.

    The program starts from the current year, which it gets from the datetime python package.
    Then, it gets all the available articles from that year.
    The year can be manually set.

    The program gets chunks of 100000 articles/journals, which is the maximum number that can be retrieved in a chunk.

    The average time for retrieving the required data is 2.67 seconds, so approximately 3 seconds per chunk.

    For every pubmed id, the script searches for the Pub Med Central id, so one can look up into the database and have access to the full free text instead of 
    just the abstract.

    Code made by: Teodor Stefan Pintea

'''

from Bio import Entrez
from Bio import Medline
import datetime
import time
import pprint

Entrez.email = "Test@example.org"
pp = pprint.PrettyPrinter(indent = 5)

year = datetime.datetime.now().year 
manual_year = 2021 # the user can manually set the year and use this variable in order to get data from the server

number = 0
results = []

start = time.perf_counter()

while(True):
    handle = Entrez.esearch(db = "pmc", term = "", mindate = str(year), maxdate = str(year), retstart = number, retmax = 1)
    record = Entrez.read(handle)
    handle.close()
    pp.pprint(record)
    break
    if(len(record['IdList']) == 0):
        print("DONE")
        break
    number += 100000
    results.extend(record['IdList'])

finish = time.perf_counter()
#print(len(results))
#print(results)
print(f'Finished in {round(finish - start, 2)} seconds(s)')



h = Entrez.efetch(db = 'pmc', id = '7829055', retmode = 'xml')
r = Medline.parse(h)
for i in r:
    pp.pprint(i)
h.close()
