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
from bs4 import BeautifulSoup
import datetime
import time
import pprint

bad = 0
good = 0

def get_pmc_ids(number = 0, requested_year = datetime.datetime.now().year):
    '''
        Method that returns a list of ids from Pub Med Central
        Starts from number = 0 as the first item is at position 0
    '''
    results = []
    while(True):
        try:
            handle = Entrez.esearch(db = "pmc", term = "", mindate = str(requested_year), maxdate = str(requested_year), retstart = number, retmax = 100000)
            record = Entrez.read(handle)
            handle.close()

            if(len(record['IdList']) == 0):
                print("DONE")
                break

            number += 100000
            results.extend(record['IdList'])
    
        except Exception as e:
            print(str(e) + " when retrieving list of ids")

    return results

def hasFullText(article):
    '''
        Method that cheks if the publishers of the article allow a user to download the full text of an article for free.
    '''

    if "The publisher of this article does not allow downloading of the full text in XML form." in article:
        return True

    return False

def data_cleaner(doc):
    '''
        Method that cleans the data. It uses BeautifulSoup in order to get the text from the XML file.
        It receives the string that contains document's XML text.

        Each text has two main setions and multiple subsections. The two main sections are ABSTRACT and BODY. The method keeps them into two name specific variables,
        so they can be processed with the sentiment analysis methods.

        The extraction of data is done in order to avoid  the binary graphic data that is at the end of the sequence by using "[-1].extract()"

        Then, the method creates two lists, containing the pragraphs of each sections. These will be used to analyze and compute the individual sentiment of both the paper
        and the sections.

        Due to the inconsistencies of papers, some may lack either abstract or body, so the method handles these exceptions. If a paper does not have both of the main sections
        then it will be considered as faulty and will not be processed by the sentiment analysis tool.

    '''
    abstract_paragraphs = list()
    body_paragraphs = list()
    try:
        soup = BeautifulSoup(doc, "html.parser")
        try:
            paper_abstract = soup.find_all('abstract')[-1].extract()
            abstract_paragraphs = list(paper_abstract.find_all('p'))
        except Exception as e:
             print(str(e) + " Paper does not have an abstract.")
        try:
            paper_body = soup.find_all('body')[-1].extract()
            body_paragraphs = list(paper_body.find_all('p'))
        except Exception as e:
             print(str(e) + " Paper does not have a body.")
    except Exception as e:
        print(str(e) + " when parsing the data with BeautifulSoup.")
        
   

    return abstract_paragraphs, body_paragraphs
    

def get_full_text(id):
    '''
        Method that returns the full text of singular bodies such as: Abstract, Introduction, Methods, Results, Discussion
        Any additional fields that might be relevant will be introduced.
        Returns the text as string or empty string if the paper is not available.
    '''
    global bad

    try:
        handle = Entrez.efetch(db = 'pmc', id = id, retmode = 'xml')
        record = handle.read().decode('utf-8')
        handle.close()

        if hasFullText(record): 
            bad = bad + 1
            return ""

    except Exception as e:
        print(str(e) + " for ID:" + str(id))

    return record
    
#in development
    
# set the email for entrez and configure the pretty printer
Entrez.email = "stevenpintea@gmail.com"
pp = pprint.PrettyPrinter(indent = 5)

#set the counter to check how long does it take for the program to run, in order to gather statistics
start = time.perf_counter()

results = get_pmc_ids()

finish = time.perf_counter()

print(f'Finished in {round(finish - start, 2)} seconds(s)')

for index, docID in enumerate(results):
    if index % 25 == 0 and index > 0:
        time.sleep(3)
        print("Avoid an error created by the script requesting information to often.")
        print("To process: " + str(len(results)))
        print("Processed: " + str(index))
        print("Unable: " + str(bad))
        print("Good   " +   str(good))

    document_text = get_full_text(id = str(docID))
    if document_text != "":
        good += 1
        abstract_paragraphs, body_paragraphs = data_cleaner(document_text)
        #if len(abstract_paragraphs) == 0 and len(body_paragraphs) != 0:
            #print(body_paragraphs)
            #break

    #if good == 12: break
    print(docID)


print("Unable   " +   str(bad))
print("Good   " +   str(good))
