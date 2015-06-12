filename = 'RedGalaxies.csv' 
from csv import DictReader
catFilename = 'sdss_catalog'
id_set = set()

catalog = []

with open(filename) as f:
    reader = DictReader(f)
    for row in reader:
        run = row['run']
        camcol = row['camcol']
        field = row['field'] 

        run = '0'*(6-len(run))+run
        field = '0'*(4-len(field))+field

        ID = ''.join([run, '-', camcol, '-', field])
        if ID in id_set:
            continue
        id_set.add(ID)

        catalog.append(' '.join([ID, row['colc'], row['rowc']]))

fullText = '\n'.join(catalog)

with open(catFilename, 'w') as f:
    f.write(fullText)
