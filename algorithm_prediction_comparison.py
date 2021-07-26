'''
    This script is used to visualise the difference between the number of papers that are seen by each algorithm as sentimental.
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

u = [2021, 2020, 2019, 2018, 2017, 2016, 2015, 2014, 2013, 2012, 2011, 2010, 2009, 2008, 2007].reverse()
years = pd.read_csv("years")
years.drop("processed", axis = 'columns', inplace=True)
years.drop("nr_papers", axis = 'columns', inplace=True)
years.drop("nr_pos_analyseSentiment", axis = 'columns', inplace=True)
years.drop("nr_neg_analyseSentiment", axis = 'columns', inplace=True)
years.drop("nr_neutral_analyseSentiment", axis = 'columns', inplace=True)
#years.drop("nr_pos_nltkSentiment", axis = 'columns', inplace=True)
#years.drop("nr_neg_nltkSentiment", axis = 'columns', inplace=True)
#years.drop("nr_neutral_nltkSentiment", axis = 'columns', inplace=True)
years.drop("nr_pos_textBlobSentiment", axis = 'columns', inplace=True)
years.drop("nr_neg_textBlobSentiment", axis = 'columns', inplace=True)
years.drop("nr_neutral_textBlobSentiment", axis = 'columns', inplace=True)

papers = 1000
x = 0
y = 0
z = 0
a = []
c = []
for i in years.itertuples():
    b = []
    d = []
    b.append(i[1])
    b.append(i[2] - x)
    b.append(i[3] - y) 
    b.append(i[4] - z)
    d.append(i[1])
    d.append(i[2] - x + i[4] - z - (i[3] - y))
    #d.append(i[3] - y) 
    if i[0] != 0:
        x = i[2]
        y = i[3]
        z = i[4]
    a.append((b))
    c.append((d))

a.reverse()
c.reverse()
df = pd.DataFrame(a, columns=['years', 'positive', 'neutral', 'negative'])
#df2 = pd.DataFrame(c, columns=['years', 'sentimental', 'neutral'])
df2 = pd.DataFrame(c, columns=['years', 'sentimental - neutral'])
df.set_index('years', inplace=True)
df2.set_index('years', inplace=True)
print(df)


# Generate sample data
c=["#5cb85c","#5bc0de","#d9534f"]
fig, ax = plt.subplots()
df2.plot(ax = plt.gca(), kind='bar',figsize=(12,8),color='y')     # plotting bar plot
plt.xlabel('Years',fontsize=14)
plt.ylabel('Difference between the papers with sentiment and the ones without',fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()

for p in ax.patches:
    ax.annotate(str(p.get_height()), (p.get_x()+0.08, p.get_height()),
                va='bottom', ha='center', weight='bold')

plt.show()
