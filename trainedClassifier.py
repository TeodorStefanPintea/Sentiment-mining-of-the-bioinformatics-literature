'''
    This is a classifier which was trained on a movie review data set and applied in the bioinformatics domain.
'''
import pandas as pd 
import random
from nltk import word_tokenize
from nltk.sentiment.util import mark_negation
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC
from sklearn.base import TransformerMixin

data = pd.read_csv("labeledTrainData.tsv", header=0, delimiter="\t", quoting=3) 
# 25000 movie reviews in data

sentiment_data = list(zip(data["review"], data["sentiment"]))
random.shuffle(sentiment_data)
 
# 80% for training
train_X, train_y = zip(*sentiment_data[:20000])
 
# Keep 20% for testing
test_X, test_y = zip(*sentiment_data[20000:])

unigram_bigram_clf = Pipeline([
    ('vectorizer', CountVectorizer(analyzer="word",
                                   ngram_range=(1, 2),
                                   tokenizer=word_tokenize,
                                   # tokenizer=lambda text: mark_negation(word_tokenize(text)),
                                   preprocessor=lambda text: text.replace("<br />", " "),)),
    ('classifier', LinearSVC())
])
 
#unigram_bigram_clf.fit(train_X, train_y)
#print(unigram_bigram_clf.score(test_X, test_y))



