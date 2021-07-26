'''
    This script uses frequency distribution to plot the top 50 used words in a text.
'''
import nltk
from nltk.probability import FreqDist
from nltk.corpus import stopwords  
from nltk.tokenize import word_tokenize

stop_words = set(stopwords.words('english'))
print("the" in stop_words)
 
all_words_used = []
def most_used_words(sentence):
    global all_words_used
    wt_words = word_tokenize(sentence)
    all_words_used += [w.lower() for w in wt_words if not w in stop_words]
    
    return all_words_used

#data_analysis = nltk.FreqDist(all_words_used)
#all_words_used = dict([(count, word) for count, word in data_analysis.items() if len(count) > 3])
#data_analysis = nltk.FreqDist(all_words_used)
#data_analysis.plot(50, cumulative=False)


 