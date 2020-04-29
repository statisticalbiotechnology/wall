import bs4
with open("sunburst.html") as inf:
    txt = inf.read()
    soup = bs4.BeautifulSoup(txt)

json = soup.find("select", {"id": "json_sources"})

new_link = soup.new_tag("option", value='null.json')
new_link.string = "hey"
json.append(new_link)
json.append("\n")

with open("sunburst.html", 'w') as outf:
    outf.write(str(soup))
