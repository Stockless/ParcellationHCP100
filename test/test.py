palabras = ['uno','dos','tres']
extension = '.txt'
with open(r'hola'+extension,'w') as doc:
	for palabra in palabras:
		doc.write("%s\n" % palabra)
	print('Listoko')