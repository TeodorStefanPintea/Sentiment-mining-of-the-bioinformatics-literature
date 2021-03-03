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


#analyser = AnalyseSentiment()

#result = analyser.Analyse("Despite numerous experimental and analytical differences across studies, the literature does not clearly support a beneficial role in exposure to negative air ions and respiratory function or asthmatic symptom alleviation.")

#result2 = analyser.Analyse("We did not confirm that waiting time was associated with worse long-term psychosocial consequences but type II error (failure to detect a true difference) might be a plausible explanation for our results. ")
#result3 = analyser.Analyse("A novel discovery ")

#print(result)
#print("------------------------")
#print(result2)
#print("------------------------")
#print(result3)