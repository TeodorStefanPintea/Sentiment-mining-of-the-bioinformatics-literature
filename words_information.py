'''
    This script extracts different information about words.
'''

from nltk import word_tokenize, pos_tag
from nltk.corpus import wordnet, stopwords 
import pandas as pd
import re


target_ngrams = pd.read_csv("targeted_ngrams")


def prepare_for_split(text):
    '''
        Method that removes begining and end of string, but dot.
        Replaces all the ;!? with dot, to split the text into sentences using only one separator.
        Removes or numbers
    '''
    cleaned_text = re.sub(r"^\W+", "", text)
    cleaned_text = re.sub(r"\W+$", "", cleaned_text)
    cleaned_text += '.'  
    cleaned_text = cleaned_text.replace(";", ".")
    cleaned_text = cleaned_text.replace("!", ".")
    cleaned_text = cleaned_text.replace("?", ".")
    cleaned_text = re.sub(r"[^a-zA-Z .-]", "", cleaned_text)
    return cleaned_text

stopWords = set(stopwords.words('english'))
years = [2021, 2020, 2019, 2018, 2017, 2016, 2015, 2014, 2013, 2012, 2011, 2010, 2009, 2008, 2007]

year_dict = dict()
for year in years:
    year_info = pd.read_csv(str(year) + "_papers")
    word_count = dict()
    for row in year_info.itertuples():
        abstract = prepare_for_split(row[3])
        body = prepare_for_split(row[4])

        paper_sentences = [sentence.lower() for sentence in abstract.split(".")] + [sentence.lower() for sentence in body.split(".")]

        for sentence in paper_sentences:
            tokenized_sentence = word_tokenize(sentence)
            tagged_sentence = (pos_tag(tokenized_sentence))
            if tagged_sentence == []:
                continue
            for word in tagged_sentence:
                if word[0] in stopWords or len(word[0]) < 3:
                    continue
                if "JJ" in word[1] or "MD" in word[1] or "NN" in word[1] or "RB" in word[1] or "VB" in word[1]:
                    if (word[0], word[1][:2]) not in word_count:
                        word_count[(word[0], word[1][:2])] = 1
                    else:
                        word_count[(word[0], word[1][:2])] += 1
    
    year_dict[year] = word_count
    print("DONE " + str(year))

pd.DataFrame.from_dict(year_dict).to_csv("years_words")
