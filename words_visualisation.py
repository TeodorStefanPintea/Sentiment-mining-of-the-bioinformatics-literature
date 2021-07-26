'''
    This program is used to create visual representations of the collected data.
'''

import pandas as pd
from nltk.corpus import words

englsish_words = set(words.words()) #englishwords

sentimental_medical_terms = pd.read_csv("sentimental_medical_terms")
print(sentimental_medical_terms.columns)
sentimental_medical_terms.drop('Unnamed: 0', axis = 'columns', inplace=True)
#sentimentalwords = list(sentimental_medical_terms['word'])
#sentimentaltags = list(sentimental_medical_terms['POS'])
smt = list(sentimental_medical_terms['term'])
print(smt[0])
smt2 = [str(w).lower() for w in smt]


def pieplot_POS_distribution():
    import matplotlib.pyplot as plt

    labels = 'Nouns', 'Verbs', 'Modals', 'Adjectives', 'Adverbs'
    sizes = [sentimentaltags.count('NN')
    , sentimentaltags.count('VB')
    , sentimentaltags.count('MD')
    , sentimentaltags.count('JJ'),
    sentimentaltags.count('RB')
    ]
    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', 'red']

    # Plot
    plt.pie(sizes, labels=labels, colors=colors,
    autopct='%1.1f%%', shadow=False, startangle=140)

    plt.axis('equal')
    plt.show()

print(sentimental_medical_terms.columns)

def wordCloud_used_terms():
    sentence = ""
    for i in smt2:
        if i in a:
            sentence += i
            sentence += " "


    from wordcloud import WordCloud, STOPWORDS, ImageColorGenerator

    import matplotlib.pyplot as plt
    wordcloud = WordCloud( background_color="black", max_font_size=50, max_words=50).generate(s)

    # Display the generated image:
    plt.imshow(wordcloud, interpolation='bilinear')
    plt.axis("off")
    plt.show()


def POS_count_per_year():
    d = { 2007 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2008 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2009 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2010 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2011 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2012 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2013 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2014 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2015 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2016 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2017 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2018 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2019 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2020 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    2021 : {'VB' : 0, 'NN' : 0, 'MD': 0, 'WR' : 0, 'JJ' : 0, 'RB': 0},
    }
    #years_words = pd.read_csv("years_words")
    #years_words.drop('Unnamed: 0', axis = 'columns', inplace=True)
    sentimental_medical_terms = pd.read_csv("common_sentimental_words_filtred")
    sentimental_medical_terms.drop('Unnamed: 0', axis = 'columns', inplace=True)
    

    for i in sentimental_medical_terms.itertuples():
        #if i[1] not in a:
        #    continue
        d[2007][i[2]] += i[17]
        d[2008][i[2]] += i[16]
        d[2009][i[2]] += i[15]
        d[2010][i[2]] += i[14]
        d[2011][i[2]] += i[13]
        d[2012][i[2]] += i[12]
        d[2013][i[2]] += i[11]
        d[2014][i[2]] += i[10]
        d[2015][i[2]] += i[9]
        d[2016][i[2]] += i[8]
        d[2017][i[2]] += i[7]
        d[2018][i[2]] += i[6]
        d[2019][i[2]] += i[5]
        d[2020][i[2]] += i[4]
        d[2021][i[2]] += i[3]

    df = pd.DataFrame.from_dict(d, orient='index')
    df.drop('WR', axis = 'columns', inplace=True)
    df.drop('MD', axis = 'columns', inplace=True)
    df = df.rename(columns = {"VB": 'verb', "NN" : 'noun', 'MD' : 'modal', 'JJ': 'adjective', 'RB' : 'adverb'})



    import matplotlib.pyplot as plt


    fig, ax = plt.subplots()
    df.plot(ax = plt.gca(), kind='bar',figsize=(12,8))     # plotting bar plot
    plt.xlabel('Years',fontsize=14)
    plt.ylabel('Number of POS Tags',fontsize=14)
    plt.legend(fontsize=14)
    plt.tight_layout()


    plt.show()


def most_used_sentimental_terms_plot():
    results = dict()
    sent_words = pd.read_csv("common_sentimental_words_filtred")
    sent_words.drop('Unnamed: 0', axis = 'columns', inplace=True)
    sent_words.drop('POS', axis = 'columns', inplace=True)
    for row in sent_words.itertuples():
        if row[1] == 'lag' or row[1] == 'pep' or row[1] == 'pan' or row[1] == 'mar' or row[1] == 'led' or row[1] == 'pig' or row[1] == 'spite' or row[1] == 'fat': continue
        s = (row[2] + row[3] + row[4] + row[5] + row[6] + row[7] + row[8] +row[9] + row[10] + row[11] + row[12] + row[13] + row[14] + row[15] + row[16])
        if row[1] in dict():
            results[row[1]] += int(s)
        else:
            results[row[1]] = int(s)
    print(len(results))
    results = dict(sorted(results.items(), key = lambda item: item[1], reverse= True))

    import matplotlib.pyplot as plt

    plt.bar(results.keys(), results.values())
    plt.xticks(rotation = 270)
    plt.yticks(rotation = 270)
    plt.tick_params(axis='x', which = 'major', labelsize = 10)
    
    plt.show()


def count_words_per_domain():
    aw = pd.read_csv("years_words")
    mw = set(pd.read_csv("filtered_med_terms")['0'])
    print(len(mw))
    w = set(aw['Unnamed: 0.1'])
    count = 0
    count2 = 0
    a = []
    count3 = 0
    for i in w:
        if i in a:
            count+=1

    for i in mw:
        if i in a:
            count2+=1

    for i in mw:
        if i in aw:
            count3+=1

    print(count)
    print(count2)
    print(count3)
    print(count + len(mw))

    a = count - count3
    plot_word_count_distribution()

def plot_word_count_distribution():
    import matplotlib.pyplot as plt

    # Data to plot
    labels = 'Common English Terms', 'Domain-Specific Terms'
    count2 = 0
    sizes = [a, count2]
    colors = ['green', 'lightskyblue']
    explode = (0, 0.1)  # explode 1st slice

    # Plot
    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
    autopct='%1.1f%%', shadow=True, startangle=140)

    plt.axis('equal')
    plt.show()


def plot_number_sent_words_per_year():
    sentwordsyear = pd.read_csv("sentimental_words_per_year")
    print(sentwordsyear.columns)
    sentwordsyear.set_index('0', inplace = True)
    sentwordsyear.drop('Unnamed: 0', axis = 'columns', inplace=True)
    sentwordsyear.drop('2', axis = 'columns', inplace=True)
    sentwordsyear = sentwordsyear.iloc[::-1]
    sentwordsyear = sentwordsyear.rename(columns = {"1": 'sentimental words'})
    print(sentwordsyear)

    import matplotlib.pyplot as plt

    c=["#5cb85c","#5bc0de","#d9534f"]
    df = pd.DataFrame(c, columns=['sentimental words'])

    fig, ax = plt.subplots()
    sentwordsyear.plot(ax = plt.gca(), kind='bar',figsize=(12,8),color='b')     # plotting bar plot
    plt.xlabel('Years',fontsize=14)
    plt.ylabel('Number of sentimenal words',fontsize=14)
    plt.legend(fontsize=14)
    plt.tight_layout()

    for p in ax.patches:
        ax.annotate(str(p.get_height()), (p.get_x()+0.08, p.get_height()),
                    va='bottom', ha='center', weight='bold')

    plt.show()