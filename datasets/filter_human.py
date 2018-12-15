#!/usr/bin/env python2.7
from lxml import etree

#parser = etree.XMLParser(remove_blank_text=True)
nsmap = {}
for event,elem in etree.iterparse("uniprot_sprot.xml",events=['start-ns']):
    ns,url = elem
    print ns,url
    nsmap[ns] = url
print nsmap
#tree = etree.parse("uniprot_sprot.xml",parser)
#root = tree.get_root()
#removes = list()
#for entry in root.findall("entry"):
#    remove = True
#    for name in entry.find("name"):
#        if name.text=="Human" or name.text=="Homo sapiens":
#            remove = False
#    for ref in  entry.find("organism").findall("dbReference"):
#        if ref.get("id")=="9606":
#           remove = False
#    if remove:
#        removes.append(entry)
#for tag in removes:
#    root.remove(tag)

#tree.write("uniprot_human.xml")
            
