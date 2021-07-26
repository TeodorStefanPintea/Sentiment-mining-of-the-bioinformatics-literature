'''
    This method uses VADER to measure the sentiment score.
    The thresholds can be changed.
'''

from vaderSentiment.vaderSentiment import SentimentIntensityAnalyzer

class AnalyseSentiment:
    def __init__(self):
        pass

    def Analyse(self, sentence):
        # parameter can be a paragraph or a sentence
        sentiment = {
            "sentence or paragraph":sentence,
            "overall_sentiment":"",
            "overall_sentiment_score":0.00,
            "scores":[]
        }
        sid_obj = SentimentIntensityAnalyzer()
        sentiment_dict = sid_obj.polarity_scores(sentence)
        if sentiment_dict['compound'] >= 0.50:
            sentiment["overall_sentiment"] = "Positive"
            sentiment["overall_sentiment_score"] = sentiment_dict['compound']
            sentiment["scores"].append({"positive":sentiment_dict['pos'], "negative":sentiment_dict['neg'],"neutral":sentiment_dict['neu']})
        elif sentiment_dict['compound'] <= - 0.50:
            sentiment["overall_sentiment"] = "Negative"
            sentiment["overall_sentiment_score"] = sentiment_dict['compound']
            sentiment["scores"].append({"positive":sentiment_dict['pos'], "negative":sentiment_dict['neg'],"neutral":sentiment_dict['neu']})
        else:
            sentiment["overall_sentiment"] = "Neutral"
            sentiment["overall_sentiment_score"] = sentiment_dict['compound']
            sentiment["scores"].append({"positive":sentiment_dict['pos'], "negative":sentiment_dict['neg'],"neutral":sentiment_dict['neu']})
        return sentiment['overall_sentiment_score']