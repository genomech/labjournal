from lib.blister import *

df = pd.DataFrame(index=range(10), columns=['spam', 'eggs'])

print(Blister.GitHubTable(df, True))
