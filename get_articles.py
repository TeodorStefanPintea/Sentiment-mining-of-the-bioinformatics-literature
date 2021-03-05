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
from nltk.probability import FreqDist
from nltk.corpus import stopwords  
from nltk.tokenize import word_tokenize
from nltk.stem import PorterStemmer
from wordcloud import WordCloud, STOPWORDS, ImageColorGenerator
from AnaliseSentiment.sentiment_AnalyseSentiment import AnalyseSentiment
from NLTKSentimentAnalyser.sentiment_nltk import SWN_Analyser
from TextBlobSentimentAnalyser.sentiment_TextBlob import TextBlobAnalyser
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import time
import re
import nltk

all_words_used = []
sentimental_sentences = []
nr_positive_AnalyseSentiment = 0
nr_negative_AnalyseSentiment = 0
nr_neutral_AnalyseSentiment = 0
nr_positive_NLTK = 0
nr_negative_NLTK = 0
nr_neutral_NLTK = 0
nr_positive_TextBlob = 0
nr_negative_TextBlob = 0
nr_neutral_TextBlob = 0



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
            para = paper_abstract.find_all('p')
            for p in para:
                abstract_paragraphs.append(str(p))
        except Exception as e:
             print(str(e) + " Paper does not have an abstract.")
        try:
            paper_body = soup.find_all('body')[-1].extract()
            para = paper_body.find_all('p')
            for p in para:
                body_paragraphs.append(str(p))
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
    try:
        handle = Entrez.efetch(db = 'pmc', id = id, retmode = 'xml')
        record = handle.read().decode('utf-8')
        handle.close()

        if hasFullText(record): 
            return ""

    except Exception as e:
        print(str(e) + " for ID:" + str(id))

    return record
 
def most_used_words(sentence):
    '''
        Method used to gather data about the words used by scientists.
    '''
    global all_words_used

    stop_words = set(stopwords.words('english'))
    porter = PorterStemmer()
    wt_words = word_tokenize(sentence)
    all_words_used += [porter.stem(w.lower()) for w in wt_words if not w.lower() in stop_words and len(porter.stem(w.lower())) >= 4]

    return all_words_used

def check_extreme_or_limited_sentence(score1, score2, score3):
    """
        This method receives 3 scores and applies a threshhold in order to decide if the sentence is either at extremes or might neede mor attention.
        Returns True if it needs to be kept in memory.
        It uses a popular vote - if 2 or three scores are in the given threshhold, then the sentance needs to be kept.
    """
    vote = 0
    if abs(score1) >= 0.35: vote += 1
    if abs(score2) >= 0.35: vote += 1
    if abs(score3) >= 0.35: vote += 1

    if vote >= 2:
        return True

    return False
    
def process_paragraphs(year, id, abstract, body):
    '''
        This method takes in the id of the paper and the two arrays which contain the paragraphs of the text.
        For each paragraph, it splits it and applies different methods in order to determine the sentiment.
        The method returns an element which will be added to the specific year's data.
    '''
    pattern = "<[^>]*>"
    regex = re.compile('[^a-zA-Z]')
    paper_sections_data = [year, id]
    one_paper_sentimental_sentences = [year, id]
    global sentimental_sentences, nr_positive_AnalyseSentiment, nr_negative_AnalyseSentiment, nr_neutral_AnalyseSentiment, nr_positive_NLTK, \
        nr_negative_NLTK, nr_neutral_NLTK, nr_positive_TextBlob, nr_negative_TextBlob, nr_neutral_TextBlob

    abstract_sentiment_by_sum_AnalyseSentiment = 0
    abstract_sentiment_by_sum_NLTK = 0
    abstract_sentiment_by_sum_TextBlob = 0
    abstract_paragraphs_AnalyseSentiment = 0
    abstract_paragraphs_NLTK = 0
    abstract_sentiment_by_average_TextBlob = 0

    #clean the data
    cleaned_paragraph = re.sub(pattern, "", str(abstract))
    paper_sections_data.append(cleaned_paragraph)
    cleaned_paragraph = re.sub(pattern, "", str(body))
    paper_sections_data.append(cleaned_paragraph)

    for paragraph in abstract:
        #clean the data
        cleaned_paragraph = re.sub(pattern, "", str(paragraph))
        paragraph_sentences = list(cleaned_paragraph.split("."))

        for sentence in paragraph_sentences:
            cleaned_sentence = regex.sub(" ", sentence)
            #get sentences score
            sentence_score_AnalyseSentiment = AnalyseSentiment().Analyse(cleaned_sentence)
            sentence_score_TextBlob = TextBlobAnalyser().analyse_sentence(cleaned_sentence)
            sentence_score_NLTK = SWN_Analyser().swn_polarity(sentence)

            #check if the sentence is worthy to be memorised and append it if that is the case
            if check_extreme_or_limited_sentence(sentence_score_AnalyseSentiment, sentence_score_TextBlob, sentence_score_NLTK):
                one_paper_sentimental_sentences.append(sentence)
                one_paper_sentimental_sentences.append(sentence_score_AnalyseSentiment)
                one_paper_sentimental_sentences.append(sentence_score_NLTK)
                one_paper_sentimental_sentences.append(sentence_score_TextBlob)
                sentimental_sentences.append(one_paper_sentimental_sentences)
                one_paper_sentimental_sentences = [year, id]

            #add the scores to another variable that will be used to compute the average of the paragraph
            abstract_sentiment_by_sum_AnalyseSentiment += sentence_score_AnalyseSentiment
            abstract_sentiment_by_sum_NLTK += sentence_score_NLTK
            abstract_sentiment_by_sum_TextBlob += sentence_score_TextBlob
        
        #calculate the paragraphs average
        abstract_paragraphs_AnalyseSentiment += (abstract_sentiment_by_sum_AnalyseSentiment / len(paragraph_sentences))
        abstract_paragraphs_NLTK += (abstract_sentiment_by_sum_NLTK / len(paragraph_sentences))
        abstract_sentiment_by_average_TextBlob += (abstract_sentiment_by_sum_TextBlob / len(paragraph_sentences))

    #calculate the abstract average score
    abstract_sentiment_by_average_AnalyseSentiment = "%.5f" % (abstract_paragraphs_AnalyseSentiment / len(abstract))
    abstract_sentiment_by_average_NLTK = "%.5f" % (abstract_paragraphs_NLTK / len(abstract))
    abstract_sentiment_by_average_TextBlob = "%.5f" % (abstract_sentiment_by_average_TextBlob / len(abstract))


    #declare what is required for body
    body_sentiment_by_sum_AnalyseSentiment = 0
    body_sentiment_by_sum_NLTK = 0
    body_sentiment_by_sum_TextBlob = 0
    body_paragraphs_AnalyseSentiment = 0
    body_paragraphs_NLTK = 0
    body_sentiment_by_average_TextBlob = 0

    for paragraph in body:
        #clean the data
        cleaned_paragraph = re.sub(pattern, "", str(paragraph))
        paragraph_sentences = list(cleaned_paragraph.split("."))

        for sentence in paragraph_sentences:
            cleaned_sentence = regex.sub(" ", sentence)
            #get sentences score
            sentence_score_AnalyseSentiment = AnalyseSentiment().Analyse(cleaned_sentence)
            sentence_score_TextBlob = TextBlobAnalyser().analyse_sentence(cleaned_sentence)
            sentence_score_NLTK = SWN_Analyser().swn_polarity(sentence)

            #check if the sentence is worthy to be memorised and append it if that is the case
            if check_extreme_or_limited_sentence(sentence_score_AnalyseSentiment, sentence_score_TextBlob, sentence_score_NLTK):
                one_paper_sentimental_sentences.append(sentence)
                one_paper_sentimental_sentences.append(sentence_score_AnalyseSentiment)
                one_paper_sentimental_sentences.append(sentence_score_NLTK)
                one_paper_sentimental_sentences.append(sentence_score_TextBlob)
                sentimental_sentences.append(tuple(one_paper_sentimental_sentences))
                one_paper_sentimental_sentences = [year, id]
            
            #add the scores to another variable that will be used to compute the average of the paragraph
            body_sentiment_by_sum_AnalyseSentiment += sentence_score_AnalyseSentiment
            body_sentiment_by_sum_NLTK += sentence_score_NLTK
            body_sentiment_by_sum_TextBlob += sentence_score_TextBlob

        #calculate the paragraphs average
        body_paragraphs_AnalyseSentiment += (body_sentiment_by_sum_AnalyseSentiment / len(paragraph_sentences))
        body_paragraphs_NLTK += (body_sentiment_by_sum_NLTK / len(paragraph_sentences))
        body_sentiment_by_average_TextBlob += (body_sentiment_by_sum_TextBlob / len(paragraph_sentences))

    #calculate the body average score
    body_sentiment_by_average_AnalyseSentiment = "%.5f" % (body_paragraphs_AnalyseSentiment / len(body))
    body_sentiment_by_average_NLTK = "%.5f" % (body_paragraphs_NLTK / len(body))
    body_sentiment_by_average_TextBlob = "%.5f" % (body_sentiment_by_average_TextBlob / len(body))


    #append the scores
    paper_sections_data.append(abstract_sentiment_by_average_AnalyseSentiment)
    paper_sections_data.append(abstract_sentiment_by_average_NLTK)
    paper_sections_data.append(abstract_sentiment_by_average_TextBlob)
    paper_sections_data.append(body_sentiment_by_average_AnalyseSentiment)
    paper_sections_data.append(body_sentiment_by_average_NLTK)
    paper_sections_data.append(body_sentiment_by_average_TextBlob)
    #compute and apend the overall for every method used
    paper_sections_data.append("%.5f" % ((float(abstract_sentiment_by_average_AnalyseSentiment) + float(body_sentiment_by_average_AnalyseSentiment)) / 2))
    paper_sections_data.append("%.5f" % ((float(abstract_sentiment_by_average_NLTK) + float(body_sentiment_by_average_NLTK)) / 2))
    paper_sections_data.append("%.5f" % ((float(abstract_sentiment_by_average_TextBlob) + float(body_sentiment_by_average_TextBlob)) / 2))


    #classify the papers into either positive, negative or neutral
    if float(paper_sections_data[10]) >= 0.5: nr_positive_AnalyseSentiment += 1
    elif float(paper_sections_data[10]) <= -0.5: nr_negative_AnalyseSentiment += 1
    else: nr_neutral_AnalyseSentiment += 1

    if float(paper_sections_data[11]) >= 0.5: nr_positive_NLTK += 1
    elif float(paper_sections_data[11]) <= -0.5: nr_negative_NLTK += 1
    else: nr_neutral_NLTK += 1

    if float(paper_sections_data[12]) >= 0.5: nr_positive_TextBlob += 1
    elif float(paper_sections_data[12]) <= -0.5: nr_negative_TextBlob += 1
    else: nr_neutral_TextBlob += 1

    return tuple(paper_sections_data), tuple(one_paper_sentimental_sentences)

def run(wanted_year):
    '''
        This is the main method that will run the data gathering script and return the results.
        It requires the year that needs to be processed and will locally output save an archive with 2 yearexcel files which contain:
            year_sections - data about the papers that could be processed, with their individual scores per section
            sentence_overall - the sentences with the most bias and the sentences which might be tuned manually (score with the absolute 0.35 or higher)
        Each excel has columns for each method used to gather the sentiment. 
        Anothe excel containing 
            year_overall - data about the overall nr of papers out of which how many are positive, negative
        The output will be processed in another script that will create the UI and plot the data for better, visual representation.

    '''
    #prepare the data that will be introduced in the dataframes at the end
    processed_papers = 0
    sentence_overall = []
    columns_sentences = ['year', 'id', 'sentence', 'analyseSentiment', 'nltkSentiment', 'textBlobSentiment']
    columns_papers_sections = ['year', 'id', 'abstract', 'body', 'analyseSentiment_abstract', 'analyseSentiment_body', \
     'nltkSentiment_abstract', 'nltkSentiment_body',  'textBlobSentiment_abstract', 'textBlobSentiment_body', \
     'overallAnalsyeSentiment', 'overallNLTK', 'overallTextBlob']
    columns_year_overall = ['year', 'nr_papers', 'processed', 'nr_pos_analyseSentiment', 'nr_neg_analyseSentiment', 'nr_neutral_analyseSentiment',\
    'nr_pos_nltkSentiment', 'nr_neg_nltkSentiment', 'nr_neutral_nltkSentiment',\
         'nr_pos_textBlobSentiment', 'nr_neg_textBlobSentiment', 'nr_neutral_textBlobSentiment']

    df_sentences = pd.DataFrame(columns = columns_sentences)
    df_year = pd.DataFrame(columns = columns_year_overall)
    df_papers = pd.DataFrame(columns = columns_papers_sections)

    # set the email for entrez and configure the pretty printer
    Entrez.email = "stevenpintea@gmail.com"
    #get all records for the specific year
    results = get_pmc_ids(requested_year= wanted_year)

    all_available_papers = len(results)
    year_overall = [wanted_year, all_available_papers]


    for index, docID in enumerate(results):
        if index % 25 == 0 and index > 0:
            time.sleep(1)
            print("Avoid an error created by the script requesting information to often.")

        document_text = get_full_text(id = str(docID))
        if document_text != "":
            abstract_paragraphs, body_paragraphs = data_cleaner(document_text)
            if len(abstract_paragraphs) != 0 and len(body_paragraphs) != 0:
                #only do something if we have both sections of a text
                paper_data, paper_sentimental_sentences = process_paragraphs(wanted_year, docID, abstract_paragraphs, body_paragraphs)

                df_papers = df_papers.append(pd.Series(paper_data, index = df_papers.columns), ignore_index = True)

                if len(paper_sentimental_sentences) > 2:
                    df_sentences = df_sentences.append(pd.Series(paper_sentimental_sentences, index = df_sentences.columns), ignore_index = True)
                
                processed_papers += 1

        if processed_papers == 10: break

    year_overall.append(processed_papers)
    year_overall.append(nr_positive_AnalyseSentiment)
    year_overall.append(nr_negative_AnalyseSentiment)
    year_overall.append(nr_neutral_AnalyseSentiment)
    year_overall.append(nr_positive_NLTK)
    year_overall.append(nr_negative_NLTK)
    year_overall.append(nr_neutral_NLTK)
    year_overall.append(nr_positive_TextBlob)
    year_overall.append(nr_negative_TextBlob)
    year_overall.append(nr_neutral_TextBlob)

    df_year = df_year.append(pd.Series(tuple(year_overall), index = df_year.columns), ignore_index = True)

    df_year.to_csv(str(wanted_year), index = None, header=True)
    df_sentences.to_csv(str(wanted_year) + '_sentences', index = None, header=True)
    df_papers.to_csv(str(wanted_year) + '_papers', index = None, header=True)

import timeit
a = timeit.default_timer()
run(2021)
b = timeit.default_timer()
print(b-a)
