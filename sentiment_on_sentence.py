from textblob import TextBlob
from collections import OrderedDict

text = """BACKGROUND: Twitter is a microblogging service where users can send and read 
short 140-character messages called "tweets." There are several unstructured, 
free-text tweets relating to health care being shared on Twitter, which is 
becoming a popular area for health care research. Sentiment is a metric commonly 
used to investigate the positive or negative opinion within these messages. 
Exploring the methods used for sentiment analysis in Twitter health care 
research may allow us to better understand the options available for future 
research in this growing field.
OBJECTIVE: The first objective of this study was to understand which tools would 
be available for sentiment analysis of Twitter health care research, by 
reviewing existing studies in this area and the methods they used. The second 
objective was to determine which method would work best in the health care 
settings, by analyzing how the methods were used to answer specific health care 
questions, their production, and how their accuracy was analyzed.
METHODS: A review of the literature was conducted pertaining to Twitter and 
health care research, which used a quantitative method of sentiment analysis for 
the free-text messages (tweets). The study compared the types of tools used in 
each case and examined methods for tool production, tool training, and analysis 
of accuracy.
RESULTS: A total of 12 papers studying the quantitative measurement of sentiment 
in the health care setting were found. More than half of these studies produced 
tools specifically for their research, 4 used open source tools available 
freely, and 2 used commercially available software. Moreover, 4 out of the 12 
tools were trained using a smaller sample of the study's final data. The 
sentiment method was trained against, on an average, 0.45% (2816/627,024) of the 
total sample data. One of the 12 papers commented on the analysis of accuracy of 
the tool used.
CONCLUSIONS: Multiple methods are used for sentiment analysis of tweets in the 
health care setting. These range from self-produced basic categorizations to 
more complex and expensive commercial software. The open source and commercial 
methods are developed on product reviews and generic social media messages. None 
of these methods have been extensively tested against a corpus of health care 
messages to check their accuracy. This study suggests that there is a need for 
an accurate and tested tool for sentiment analysis of tweets trained using a 
health care setting-specific corpus of manually annotated tweets first. """

blob = TextBlob(text)
textSentences = blob.sentences
sentanceSentiment = dict()

for index, sentence in enumerate(textSentences):
    sentanceSentiment[index] =  (sentence, sentence.sentiment.polarity, sentence.sentiment.subjectivity)

result = sorted(sentanceSentiment, key = lambda x : sentanceSentiment[x][1], reverse= True)

for index in result:
    print(sentanceSentiment[index], end = '\n\n')