import urllib2


def getRunXML( runNum ) :
	
	user_agent = 'Mozilla/5.0 (Macintosh; U; Intel Mac OS X 10_6_4; en-US) AppleWebKit/534.3 (KHTML, like Gecko) Chrome/6.0.472.63 Safari/534.3'
	headers = { 'User-Agent' : user_agent }
	req = urllib2.Request('http://online.star.bnl.gov/RunLog/index.php?r=' + str(runNum) + '&method=xml', None, headers)
	response = urllib2.urlopen(req)
	page = response.read()

	return page