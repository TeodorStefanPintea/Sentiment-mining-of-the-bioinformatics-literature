'''
    This script is used to measure identify the subjective sentences in a paper.
'''

import pandas as pd 
import re

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

l = ['i prove', 'i discover', 'i proved', 'i discovered', 'i invented', 'i approach', 
     'we thus propose a novel', 'we thus propose', 'our paper proposes', 'our paper presents' ,'we prove', 'we discover', 'we proved',
     'we invented', 'we discovered', 'we used a novel', 
     'our study proves', 'our study dscovered', 'our study invented', 'our study proved', 
     'our study created', 'our study demolished',
     'our paper proves', 'our paper dscovered', 'our paper invented', 'our paper proved',
     'our paper created', 'our paper demolished',
     'our paper proves', 'our paper dscovered', 'our paper invented', 'our paper proved',
     'our paper created', 'our paper demolished'
     'our measurement proves',  'our measurement dscovered', 
     'our measurement invented', 'our measurement proved', 
     'our measurement created', 'our measurement demolished',
     'our results invented', 'our results proved',
     'our results created', 'our results demolished',
     'our finding proves',  'our finding dscovered', 
     'our finding invented', 'our finding proved',
     'our findings prove', 'our findings make', 'our findings dscovered', 
     'our findings invented', 'our findings proved',
     'the paper proves',  'the paper dscovered', 'the paper invented', 'the paper proved',
     'the paper created', 'the paper demolished',
     'our results prove', 'our results discovered', 
     "we did not", "our values are", "our results", "our data", "he has obtained", "the comparison of values shows", "the data presented", "this value", "this series"
     ]

def now():
    all_years = list()
    for year in years:
        count = 0
        sents = set()
        year_info = pd.read_csv(str(year) + "_papers")
        for row in year_info.itertuples():
            count2 = 0
            abstract = prepare_for_split(row[3])
            body = prepare_for_split(row[4])
            abs_sentences = list(abstract.lower().split("."))
            body_sentences = list(body.lower().split("."))
            #count += len(abs_sentences) + len(body_sentences)
            for i in abs_sentences:
                for word in l:
                    if word in i:
                        #sents.add(i)
                        count2 += 1
                        break
            for i in body_sentences:
                for word in l:
                    if word in i:
                        #sents.add(i)
                        count2 += 1
                        break
            if count2 > 0: count +=1

        all_years.append((year, count))
        #all_years.append((year, len(sents), count))
    all_years.reverse()
    df = pd.DataFrame(all_years, columns=['years', 'papers with subjective remarks'])
    df.set_index('years', inplace=True)

    return df


df2 = now()
import numpy as np
import matplotlib.pyplot as plt
c=["#5cb85c","#5bc0de","#d9534f"]
fig, ax = plt.subplots()
df2.plot(ax = plt.gca(), kind='bar',figsize=(12,8),color=c)     # plotting bar plot
plt.xlabel('Years',fontsize=14)
plt.ylabel('Number of subjective papers',fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()

for p in ax.patches:
    ax.annotate(str(p.get_height()), (p.get_x()+0.08, p.get_height()),
                va='bottom', ha='center', weight='bold')

plt.show()


