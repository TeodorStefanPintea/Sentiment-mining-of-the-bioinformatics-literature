'''
    Script used for filtering the medical terms based on some criteria.
'''
import pandas as pd 
from nltk.corpus import wordnet, stopwords 

years = [2021, 2020, 2019, 2018, 2017, 2016, 2015, 2014, 2013, 2012, 2011, 2010, 2009, 2008, 2007]
stopWords = set(stopwords.words('english'))
all_medical_terms = list()
for year in years:
    med_terms_year = pd.read_csv("medical_terms_" + str(year))["term"].tolist()
    medical_terms = [str(term) for term in  med_terms_year]
    medical_terms = [term for term in  medical_terms if len(term) > 3]
    medical_terms = [term for term in  medical_terms if " " not in term]
    medical_terms = [term for term in  medical_terms if term not in stopWords]
    medical_terms = [term for term in  medical_terms if not any(c.isupper() for c in term)]
    medical_terms = [term for term in  medical_terms if "." not in term]
    medical_terms = [term for term in  medical_terms if "-" not in term]
    medical_terms.sort()
    all_medical_terms += list(set(medical_terms))
    all_medical_terms = list(set(all_medical_terms))
    print("DONE " + str(year))

print(len(all_medical_terms))
pd.DataFrame(list(all_medical_terms)).to_csv("filtered_med_terms")


