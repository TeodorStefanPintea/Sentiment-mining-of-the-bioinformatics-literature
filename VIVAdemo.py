'''
    This is the demo for the VIVA presentation.
'''


print("Importing the required packages.")

from NLTKSentimentAnalyser import sentiment_nltk
import pandas as pd

print("Creating the analyser object.")

sa = sentiment_nltk.SWN_Analyser()

sent = ["He has a terminal disease." ,
        "The mortality rate has dropped in the past 5 years." ,
        "Our discovery can cure cancer."
        ]

print()
print("Loading the Subjective Lexicon")

df = list(pd.read_csv("SubjectiveLexicon")['0'])

print()



print("Processing the sentences:")
print()
for sentence in range(3):
    print()
    print("Sentiment identified: ")
    print("Sentence: '" + sent[sentence] + "' is: ", end=" ")
    print(sa.swn_polarity(sent[sentence]))

for sentence in range(3):
    ok = True
    print()
    print("Subjectivity level: ")
    print("Sentence: '" + sent[sentence] + "' is: ", end=" ")
    for term in df:
        if term in sent[sentence].lower():
            print("SUBJECTIVE")
            ok = False
            break
    if ok:
        print("NOT SUBJECTIVE")
    
print()
