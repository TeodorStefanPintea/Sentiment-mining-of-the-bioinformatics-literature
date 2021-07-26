'''
    This method is an implementation of the open source TextBlob library.
'''

from textblob import TextBlob
from collections import OrderedDict

class TextBlobAnalyser:
    def __init__(self):
        pass

    def analyse_sentence(self, sentence):
        # parameter can be a paragraph or a sentence
        blob = TextBlob(sentence)
        return float(blob.sentiment.polarity)

    def words_assessed(self, sentence):
        blob = TextBlob(sentence)
        return blob.sentiment_assessments
