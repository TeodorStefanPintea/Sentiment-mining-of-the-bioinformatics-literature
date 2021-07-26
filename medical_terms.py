import scispacy
import spacy
import pandas as pd
import re

nlp = spacy.load("en_core_sci_sm")
years = [2021, 2020, 2019, 2018, 2017, 2016, 2015, 2014, 2013, 2012, 2011, 2010, 2009, 2008, 2007]

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
    cleaned_text = re.sub(r"[^a-zA-Z .-]", " ", cleaned_text)
    return cleaned_text

def cont_medical_terms_per_year():
    for year in years:
        year_info = pd.read_csv(str(year) + "_papers")
        medical_terms = {"term" : list()}
        set_words = set()
        for row in year_info.itertuples():
            abstract = prepare_for_split(row[3])
            body = prepare_for_split(row[4])
            doc_abs = nlp(abstract)
            for med_term in doc_abs.ents:
                set_words.add(med_term)
            
            doc_body = nlp(body)
            for med_term in doc_body.ents:
                set_words.add(med_term)
            
            

        medical_terms["term"] = list(set_words)
        print(medical_terms)
        

def visualize():
    import pandas as pd

    data = pd.read_excel("sents4.xlsx")
    positive_words = open("positive-words.txt").read().split()[174:]
    negative_words = open("negative-words.txt").read().split()[175:]
    prediction = list()
    actual = list()
    medical_terms = {"term" : list()}
    set_words = set()

    for i in data.itertuples():
        actual.append(i[3])
        med_text = nlp(i[2])
        for med_term in med_text.ents:
            set_words.add(med_term)
            
    for i in data.itertuples():
            a = i[2][:-1]
            a = a.split()
            sc = 0
            for w in a:
                if w.lower() not in set_words:
                    if w.lower() in positive_words:
                        sc += 1
                    elif w.lower() in negative_words:
                        sc -= 1
            #if sc < 0:
                #prediction.append(-1)
            #elif sc > 0:
                #prediction.append(1)
            #else:
                #prediction.append(0) 
            if abs(sc) > 0:
                prediction.append(1)
            else:
                prediction.append(-1)      

    print(prediction)
    print(actual)

    from sklearn import metrics

    print(metrics.confusion_matrix(actual, prediction))
    print(metrics.classification_report(actual, prediction))

