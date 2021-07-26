'''
    This program is used to visualise year-to-year data change.
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


def get_info():
    all_years = list()
    for year in years:
        count = 0
        year_info = pd.read_csv(str(year) + "_papers")
        for row in year_info.itertuples():
            ok = True
            abstract = prepare_for_split(row[3])
            body = prepare_for_split(row[4])
            abs_sentences = list(abstract.lower().split("."))
            body_sentences = list(body.lower().split("."))

            for i in abs_sentences:
                if ok == False: break
                if "publication bias" in i:
                    count += 1
                    ok = False
                    break
            
            for i in body_sentences:
                if ok == False: break
                if "publication bias" in i:
                    count += 1
                    ok = False
                    break

        all_years.append((year, count))

    pd.DataFrame(all_years, columns=["Year", "Nr"]).to_csv("pub_bias2")

def growing():
    data = pd.read_csv("pub_bias").drop(columns = 'Unnamed: 0')
    data2 = pd.read_csv("pub_bias2").drop(columns = 'Unnamed: 0')

    avg0 = data['Nr'].iloc[0:5].mean(axis=0)
    avg1 = data['Nr'].iloc[5:10].mean(axis=0)
    avg2 = data['Nr'].iloc[10:15].mean(axis=0)


    avg3 = data2['Nr'].iloc[0:5].mean(axis=0)
    avg4 = data2['Nr'].iloc[5:10].mean(axis=0)
    avg5 = data2['Nr'].iloc[10:15].mean(axis=0)

    all_y = [str(year) for year in years]

    avg = [avg2, avg1, avg0]
    avg2 = [avg5, avg4, avg3]

    y2 = '2007-2011'
    y1 = '2012-2016'
    y0 = '2017-2021'

    y = [y2, y1, y0]

    import matplotlib.pyplot as plt

    f, ax = plt.subplots(nrows = 1, ncols=2, sharex=True, figsize = (15,5))
    ax[0].set_title("Average number of times 'publication bias' appears in text every 5 years")
    ax[1].set_title("Avarage number of papers that refer to 'publication bias' every 5 years")
    ax[0].bar(y, avg)
    ax[1].bar(y, avg2, color = 'r')
    #f.suptitle("The evolution of 'publication bias' in bioinformatics literature during the past 15 years")
    f.tight_layout(pad = 3.0)
    plt.show()

def get_info2():
    all_years = list()
    for year in years:
        count = 0
        year_info = pd.read_csv(str(year) + "_papers")
        for row in year_info.itertuples():
            ok = True
            abstract = prepare_for_split(row[3])
            body = prepare_for_split(row[4])
            abs_sentences = list(abstract.lower().split("."))
            body_sentences = list(body.lower().split("."))

            for i in abs_sentences:
                if ok == False: break
                if "genes" in i:
                    count += 1
                
            for i in body_sentences:
                if ok == False: break
                if "genes" in i:
                    count += 1

        all_years.append((year, count))

    return all_years


def get_word_count():
    covid_2021 = get_info2()

    covid_2021.reverse()
    print(covid_2021)
    nr_papers = [value[1] for value in covid_2021]
    years.reverse()

    import matplotlib.pyplot as plt
    plt.xticks(years)
    plt.tick_params(axis = 'x', which = 'major', labelsize = 8)
    plt.margins(x=0)
    plt.plot(years, nr_papers)

    plt.xlabel("Years")
    plt.ylabel("Number of times the term 'genes' has appeared")
    plt.show()
