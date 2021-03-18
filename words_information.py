from nltk import word_tokenize, pos_tag
from nltk.stem import WordNetLemmatizer
from nltk.corpus import wordnet, stopwords 
from spacy import displacy
from scispacy.abbreviation import AbbreviationDetector
from scispacy.umls_linking import UmlsEntityLinker
import pandas as pd
import re
import scispacy
import spacy


target_ngrams = pd.read_csv("targeted_ngrams")
nlp = spacy.load("en_core_sci_sm")

pos_dict = {'JJ' : dict(), "MD" : dict(), "NN" : dict(), "RB" : dict(), "VB" : dict()}
word_count = dict()

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
    cleaned_text = re.sub(r"[^a-zA-Z.-]", " ", cleaned_text)
    return cleaned_text

stopWords = set(stopwords.words('english'))
years = [2021, 2020, 2019, 2018, 2017, 2016, 2015, 2014, 2013, 2012, 2011, 2010, 2009, 2008, 2007]
years_words = pd.DataFrame()
medical_terms = set()

for year in years:
    year_info = pd.read_csv(str(year) + "_papers")
    year_dict = dict()

    for row in year_info.itertuples():
        abstract = prepare_for_split(row[3])
        body = prepare_for_split(row[4])
        doc_abs = nlp(abstract)
        for med_term in doc_abs.ents:
            medical_terms.add(med_term)

        doc_body = nlp(body)
        for med_term in doc_body.ents:
            medical_terms.add(med_term)

        paper_sentences = [sentence.lower() for sentence in abstract.split(".")] + [sentence.lower() for sentence in body.split(".")]

        for sentence in paper_sentences:
            tokenized_sentence = word_tokenize(sentence)
            tagged_sentence = (pos_tag(tokenized_sentence))
            if tagged_sentence == []:
                continue
            for word in tagged_sentence:
                if word[0] in medical_terms or word[0] in stopWords or len(word[0]) < 3:
                    continue
                if "JJ" in word[1] or "MD" in word[1] or "NN" in word[1] or "RB" in word[1] or "VB" in word[1]:
                    lemmatizer = WordNetLemmatizer()
                    if "MD" not in word[1]:
                        if "JJ" in word[1]:
                            pos = wordnet.ADJ
                        elif "NN" in word[1]:
                            pos = wordnet.NOUN
                        elif "VB" in word[1]:
                            pos = wordnet.VERB
                        else:
                            pos = wordnet.ADV
                        word_lemma = lemmatizer.lemmatize(word[0], pos)
                    else:
                        word_lemma = lemmatizer.lemmatize(word[0])
                    if len(word_lemma) < 3 or word_lemma in stopWords or word_lemma in medical_terms:
                        continue
                    if word_lemma not in word_count:
                        word_count[word_lemma] = 1
                    else:
                        word_count[word_lemma] += 1
                    pos_dict[word[1][:2]] = word_count

    year_dict[year] = pos_dict
    print("DONE " + str(year))

pd.DataFrame(list(medical_terms)).to_csv("medical_terms")
pd.DataFrame.from_dict(year_dict).to_csv("years_words")
