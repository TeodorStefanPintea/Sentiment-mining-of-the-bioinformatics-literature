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

def get_pmc_ids(number = 0, requested_year = datetime.datetime.now().year):
    '''
        Method that returns a list of ids from Pub Med Central
        Starts from number = 0 as the first item is at position 0
    '''
    results = []
    while(True):
        handle = Entrez.esearch(db = "pmc", term = "", mindate = str(requested_year), maxdate = str(requested_year), retstart = number, retmax = 100000)
        record = Entrez.read(handle)
        handle.close()

        if(len(record['IdList']) == 0):
            print("DONE")
            break

        number += 100000
        results.extend(record['IdList'])
    
    return results

# in development
def get_full_text(id):
    '''
        Method that returns the full text of singular bodies such as: Abstract, Introduction, Methods, Results, Discussion
        Any additional fields that might be relevant will be introduced.


    '''

    handle = Entrez.efetch(db = 'pmc', id = id, retmode = 'xml')
    record = Medline.parse(handle)
    for i in record:
        pp.pprint(i)
    handle.close()
    
# set the email for entrez and configure the pretty printer
Entrez.email = "Test@example.org"
pp = pprint.PrettyPrinter(indent = 5)

#set the counter to check how long does it take for the program to run, in order to gather statistics
start = time.perf_counter()

results = get_pmc_ids()

finish = time.perf_counter()

print(len(results))

print(f'Finished in {round(finish - start, 2)} seconds(s)')




