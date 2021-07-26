'''
    SentiWordNet.txt contains a list with thw words from WordNet. Each word has a positive or negative score.
    This classifier uses SentiWordNet to get the most used sense of a word and calculates the overall sentiment of a sentence.
'''

from nltk.stem import WordNetLemmatizer
from nltk.corpus import wordnet as wn
from nltk.corpus import sentiwordnet as swn
from nltk import word_tokenize, pos_tag
 
class SWN_Analyser():
    def __init__(self):
        pass 

    def convert_tag(self, tag):
        """
            Convert PennTreebank tags to simple Wordnet tags
        """
        if tag.startswith('J'):
            return 'a'
        elif tag.startswith('N'):
            return 'n'
        elif tag.startswith('R'):
            return 'r'
        elif tag.startswith('V'):
            return 'v'
        return None

    def swn_polarity(self, sentence):
        """
            Method used to calculate the sentiment of a sentence.
            0 is neutral, -1 is negative and +1 is positive
            The threshold can be manually modified.
        """
        lemmatizer = WordNetLemmatizer()
        sentiment = 0.0
        tokens_count = 0
        tagged_sentence = pos_tag(word_tokenize(sentence))
        for word, tag in tagged_sentence:
            wn_tag = self.convert_tag(tag)
            if wn_tag not in ('n', 'a', 'r'):
                continue
            lemma = lemmatizer.lemmatize(word, pos=wn_tag)
            if not lemma:
                continue
            synsets = wn.synsets(lemma, pos=wn_tag)
            if not synsets:
                continue
            # Use the most common sense
            synset = synsets[0]
            swn_synset = swn.senti_synset(synset.name())
            sentiment += swn_synset.pos_score() - swn_synset.neg_score()
            tokens_count += 1

        if not tokens_count:
            #default is neutral
            return 0


        if sentiment > 0.1: return "POSITIVE"
        elif sentiment <= 0.1: return "NEGATIVE"

