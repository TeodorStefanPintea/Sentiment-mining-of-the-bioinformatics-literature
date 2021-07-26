'''
    This script is used to assign extreme scores and create a small data set that is used in the testing phase.
'''

import pandas as pd
from AnaliseSentiment.sentiment_AnalyseSentiment import AnalyseSentiment
from nltk.tokenize import word_tokenize

#years = [ 2021, 2020, 2019, 2018, 2017, 2016, 2015, 2014, 2013, 2012, 2011, 2010, 2009, 2008, 2007]
year = 2021
df = pd.read_csv(str(year)+"_papers")
top_sentimental_sentences = []
sentences = set()

for row in df.itertuples():
    abstract = row[3]
    body = row[4]
    abs_sentences = list(abstract.split("."))
    body_sentences = list(body.split("."))

    for sentence in abs_sentences:
        if sentence in sentences:
            continue
        sentences.add(sentence)
        if AnalyseSentiment().Analyse(sentence) >= 0.90:
            continue
            top_sentimental_sentences.append((sentence, 1))
        elif AnalyseSentiment().Analyse(sentence) <= -0.90:
            top_sentimental_sentences.append((sentence, 0))

    for sentence in body_sentences:
        if sentence in sentences:
            continue
        sentences.add(sentence)
        if AnalyseSentiment().Analyse(sentence) >= 0.90:
            continue
            top_sentimental_sentences.append((sentence, 1))
        elif AnalyseSentiment().Analyse(sentence) <= -0.90:
            top_sentimental_sentences.append((sentence, 0))

    print(len(top_sentimental_sentences))

    if len(top_sentimental_sentences) >= 20:
         df_sents = pd.DataFrame.from_records(top_sentimental_sentences, columns = ["Sentence", "Score"])
         break


print(df_sents)


df_sents.to_excel("sents2.xlsx")

