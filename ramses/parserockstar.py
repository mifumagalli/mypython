def parserockstar(halocat):
    """Read a ascii rockstar halo catalogue and return a table"""

    
    from astropy.table import Table
    import re

    #load halo catalogues
    hc = open(halocat, 'r')
    halo=[]

    #read header
    header=hc.readline()
    header = re.sub('#', '', header)
    header=header.split(" ")
    
    for line in hc:
        if "#" not in line:
            
            lines=(line.rstrip('\n')).split(" ")
            fllist = [float(x) for x in lines]
            halo.append(fllist)
            
            
    hc.close()
    haltab = Table(rows=halo,names=header)

    
    return haltab





