import urllib.request,urllib.parse,json
for term in ['AFE_0631','WP_012536247.1','RU820_RS03045']:
    u='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':term,'format':'json','pageSize':3,'resultType':'core'})
    try:
        d=json.load(urllib.request.urlopen(u,timeout=30))
        print(term,d['hitCount'],[(x.get('title'),x.get('doi'),x.get('year')) for x in d['resultList']['result']])
    except Exception as e:print(term,repr(e))
