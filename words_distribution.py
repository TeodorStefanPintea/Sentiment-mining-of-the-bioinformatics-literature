from nltk.corpus import words
import pandas as pd


def general_word_distribution():
        
    all_words = pd.read_csv("years_words")
    med_words = pd.read_csv("filtered_med_terms")
    word_list = set(set(all_words['Unnamed: 0.1']) & set(med_words['0']))
    english_words = words.words()
    english_words = set(english_words)
    common_words = set(set(all_words['Unnamed: 0.1']).difference(set(med_words['0'])))
    common_english_words = set(english_words & common_words)
    positive_words = open("positive-words.txt").read().split()[174:]
    negative_words = open("negative-words.txt").read().split()[175:]
    sentimental_words = set(positive_words + negative_words) #only 3 words are both positive and negative
    english_positives = (len(positive_words))
    english_negatives = (len(negative_words))
    #above contain medical words like cancer
    sentimental_words_no_med = set(sentimental_words.difference(set(med_words['0'])))
    papers_sentimental_words = set(common_english_words & sentimental_words_no_med)
    common_english_words_no_med = set(common_english_words.difference(set(med_words['0'])))

    nr_words_intersection = len(word_list)
    nr_all_words_used = len(set(all_words['Unnamed: 0.1']))
    nr_medical_words = len(set(med_words['0']))
    nr_common_words_collected = len(common_words)
    nr_english_words = len(common_english_words)
    nr_en_sentimental = len(sentimental_words)
    nr_en_sent_no_med = len(sentimental_words_no_med)
    nr_sentimental_in_paper = len(papers_sentimental_words)

    return papers_sentimental_words, common_english_words_no_med


def sentimental_words_dist_per_year(years):
    data = list()
    for year in years:
        all_words = pd.read_csv("years_words")
        nr_words = 0 
        nr_sentimental_words_year = 0
        for word, times in zip(all_words['Unnamed: 0.1'], all_words[str(year)]):
            nr_words += times
            if word in papers_sentimental_words:
                nr_sentimental_words_year += times
        data.append((year, int(nr_sentimental_words_year), int(nr_words)))
        print("DONE " + str(year))
    pd.DataFrame(data).to_csv("sentimental_words_per_year")


def most_common_words_and_sentimental_with_pos_differentiation():
    common_words_filtred = list()
    all_words = pd.read_csv("years_words")

    for row in all_words.itertuples():
        if 0 in row:
            continue
        elif row[2] in papers_sentimental_words:
            common_words_filtred.append(row[2:])

    pd.DataFrame(common_words_filtred, columns= cols).to_csv("common_sentimental_words_filtred")

papers_sentimental_words, com_words = general_word_distribution()
years = [2021, 2020, 2019, 2018, 2017, 2016, 2015, 2014, 2013, 2012, 2011, 2010, 2009, 2008, 2007]
cols = ["word", "POS", "2021", "2020", "2019", "2018", "2017", "2016", "2015", "2014", "2013", "2012", "2011", "2010", "2009", "2008", "2007"]


