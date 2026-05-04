import os
import argparse
import requests
import xmltodict
import numpy as np
import pyvo
from pyvo.dal import tap
from astropy.table import Table, vstack
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from email.message import Message

def download_file(url, local_filename):
    print('Downloading file from {}'.format(url))
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print('File saved to {}'.format(local_filename))
    
    
def parse_xml_links(xmllinksfile):

    with open(xmllinksfile, 'r') as f:
        doc = (xmltodict.parse(f.read()))

    links_data = doc['VOTABLE']['RESOURCE'][0]['TABLE']
    #print(links_data)

    fields = links_data['FIELD']
    field_names = [f['@name'] for f in fields]

    data_rows = links_data['DATA']['TABLEDATA']['TR']

    table = Table(names=field_names, dtype=[object]*len(field_names)) # Initialize with object dtype for all columns to handle mixed types and None

    #Add each row in list_of_rows
    for row_data in data_rows:
        table.add_row(row_data['TD'])

    #Cast the columns that you can into floats
    for nn in (table.colnames):
      try:
         table[nn]=table[nn].astype(float)
      except:
        table[nn]=table[nn].astype(str)

    return table
 
def get_eso_night(mjd_obs):
    """
    Calculates the calendar date (Year-Month-Day) corresponding to the
    start of the observing night (Noon UT).

    Args:
        mjd_obs (float): The Modified Julian Date of the observation.

    Returns:
        str: The calendar date (YYYY-MM-DD) of the night's start.
        float: The MJD integer corresponding to the night's start (MJD_night).
    """
    # 1. Create an astropy Time object from the MJD observation value
    t_obs = Time(mjd_obs, format='mjd', scale='utc')

    # 2. Apply the half-day shift (0.5 days)
    # This shifts the 'midnight' boundary to 'noon'
    t_shifted = t_obs + 0.5 * u.day

    # 3. Truncate to the integer MJD (This finds the MJD of the night's start)
    mjd_night = np.floor(t_shifted.mjd)

    # 4. Convert the integer MJD back to a calendar date string
    # Create a new Time object from the MJD_night integer
    t_night = Time(mjd_night, format='mjd', scale='utc')

    # Use the 'iso' format and only keep the date part
    date_str = np.array([t.iso.split()[0] for t in t_night])

    return date_str


def download_archive_asset(url, datadir, filename=None, type=None, force_download=False):

        if datadir[-1] != '/':
            datadir += '/'

        # Set a User Agent (modify as you like, but please let intact the python version used for our usage statistics):
        headers = {}
        headers={'User-Agent': '%s (ESO script drc %s)'%(requests.utils.default_headers()['User-Agent'], "DOWNLOAD_MUSE")}

        response = requests.get(url, stream=True, headers=headers)

        if filename == None:
            contentdisposition = response.headers.get('Content-Disposition')
            if contentdisposition != None:
                m = Message()
                m['content-type'] = contentdisposition
                filename = m.get_param('filename')

            if filename == None:
                # last chance: get anything after the last '/'
                filename = url[url.rindex('/')+1:]

        if filename.startswith('association'):
            filename  = filename.split('=')[1][:-5]+'.xml'
            print(filename)

        download = True
        if os.path.exists(datadir+filename) and not force_download:
            download = False

        if response.status_code == 200 and download:
            with open(datadir+filename, 'wb') as f:
                for chunk in response.iter_content(chunk_size=50000):
                    f.write(chunk)

            if type is not None:
               print("File {} ({}) downloaded".format(filename,type))
            else:
               print("File {} downloaded".format(filename))

        return (response.status_code, filename)


#Use this if you want to work per night batches
def download_batch(ob_info, night, datadir):

        # Decide what you want to download:
        mode_requested = "raw2raw"  # other choice: raw2raw

        rawframes = ob_info[ob_info['night'] == night]

        # --- Downloading those science raw frames and the calibration files necessary to calibrate them

        nfiles = 0
        calib_association_trees = []

        print("Downloading the {} science files and runing cal selector in {} mode.".format(len(rawframes), mode_requested))

        #  ... 1.- RETRIEVE THE SCIENCE RAW FRAME, AND RUN CALSELECTOR FOR EACH OF THEM:
        # ---------------------------------------
        caltables = []

        for raw in rawframes:
            sci_url = "http://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?{}&eso_download=file".format(raw["dp_id"])
            status, sci_filename = download_archive_asset(sci_url, datadir, force_download=False)

            calselector_url = "http://archive.eso.org/calselector/v1/associations?dp_id={}&mode={}&responseformat=votable".format(raw["dp_id"], mode_requested)
            datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(calselector_url)

            caltables.append(Table(datalink.to_table()))

        #Stack all the caltables together
        full_cal_table = vstack(caltables)
        #Keep only unique entries of the access_url column
        full_cal_table = full_cal_table.group_by('access_url').groups.aggregate(lambda x: x[0])
        #Remove rows that have #this in semantics
        full_cal_table = full_cal_table[full_cal_table['semantics'] != '#this']
        #Write the full calibration table to file
        #full_cal_table.write('FULL_CAL_TABLE.txt', format='ascii.fixed_width', delimiter=' ', bookend=False, overwrite=True)

        #Count how many entries of the eso_category column are ASSOCIATION_TREE
        n_association_trees = len(full_cal_table[full_cal_table['eso_category'] == 'ASSOCIATION_TREE'])
        print("Found {} XML ASSOCIATION FILES".format(n_association_trees))

        for type in ['BIAS', 'DARK', 'FLAT', 'SKYFLAT', 'STD', 'ILLUM', 'ARC']:
            n_type = len(full_cal_table[full_cal_table['eso_category'] == type])
            print("Found {} {} calibration files".format(n_type, type))

        #  DOWNLOAD EACH CAL FILE FROM THE LIST
        # ---------------------------------------

        for url,category in zip(full_cal_table['access_url'],full_cal_table['eso_category']):
            status, filename = download_archive_asset(url, datadir, type=category, force_download=False)

#Use this if you want to download a single OB
def download_ob(ob_id, datadir, ESO_TAP_OBS = "http://archive.eso.org/tap_obs"):

        # Decide what you want to download:
        mode_requested = "raw2raw"  # other choice: raw2raw

           # --- instantiate the ESO TAP service for raw and processed data:


        tapobs = tap.TAPService(ESO_TAP_OBS)

        # --- Query TAP for science raw
        print()
        print("Looking for SCIENCE frames for a given OB ID." )
        print("Querying the ESO TAP service at %s" %(ESO_TAP_OBS))

        query="SELECT dp_id, mjd_obs, exp_start "
        query+="FROM dbo.raw "
        query+="WHERE ob_id={} ".format(ob_id)
        query+="AND dp_cat='SCIENCE'"

        #print("")
        #print(query)
        #print("")

        rawframes = tapobs.search(query=query)
        print(rawframes.to_table())
        print("")

        # --- Downloading those science raw frames and the calibration files necessary to calibrate them

        nfiles = 0
        calib_association_trees = []

        print("Downloading the {} science files and runinng cal selector in {} mode.".format(len(rawframes.to_table()), mode_requested))

        #  ... 1.- RETRIEVE THE SCIENCE RAW FRAME, AND RUN CALSELECTOR FOR EACH OF THEM:
        # ---------------------------------------
        caltables = []

        for raw in rawframes:
            sci_url = "http://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?{}&eso_download=file".format(raw["dp_id"])
            status, sci_filename = download_archive_asset(sci_url, datadir, force_download=False)

            calselector_url = "http://archive.eso.org/calselector/v1/associations?dp_id={}&mode={}&responseformat=votable".format(raw["dp_id"], mode_requested)
            datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(calselector_url)

            caltables.append(Table(datalink.to_table()))

        #Stack all the caltables together
        full_cal_table = vstack(caltables)
        #Keep only unique entries of the access_url column
        full_cal_table = full_cal_table.group_by('access_url').groups.aggregate(lambda x: x[0])
        #Remove rows that have #this in semantics
        full_cal_table = full_cal_table[full_cal_table['semantics'] != '#this']
        #Write the full calibration table to file
        #full_cal_table.write('FULL_CAL_TABLE.txt', format='ascii.fixed_width', delimiter=' ', bookend=False, overwrite=True)

        #Count how many entries of the eso_category column are ASSOCIATION_TREE
        n_association_trees = len(full_cal_table[full_cal_table['eso_category'] == 'ASSOCIATION_TREE'])
        print("Found {} XML ASSOCIATION FILES".format(n_association_trees))

        for type in ['BIAS', 'DARK', 'FLAT', 'SKYFLAT', 'STD', 'ILLUM', 'ARC']:
            n_type = len(full_cal_table[full_cal_table['eso_category'] == type])
            print("Found {} {} calibration files".format(n_type, type))

        #  DOWNLOAD EACH CAL FILE FROM THE LIST
        # ---------------------------------------

        for url,category in zip(full_cal_table['access_url'],full_cal_table['eso_category']):
            status, filename = download_archive_asset(url, datadir, type=category, force_download=False)


def query_ob_info(ob_ids, ESO_TAP_OBS = "http://archive.eso.org/tap_obs"):

        tapobs = tap.TAPService(ESO_TAP_OBS)

        obtables = []

        for ob_id in ob_ids:

            # --- Query TAP for science raw
            print()
            print("Looking for SCIENCE frames for a given OB ID." )
            print("Querying the ESO TAP service at %s" %(ESO_TAP_OBS))

            query="SELECT dp_id, ob_id, mjd_obs, exp_start "
            query+="FROM dbo.raw "
            query+="WHERE ob_id={} ".format(ob_id)
            query+="AND dp_cat='SCIENCE'"

            #print("")
            #print(query)
            #print("")

            rawframes = tapobs.search(query=query)
            obtables.append(Table(rawframes.to_table()))

        obtable = vstack(obtables)
        obtable['night'] = get_eso_night(obtable['mjd_obs'])
        #print(obtable)

        return obtable


#------------ MAIN PROGRAM ----------

def main(input_table, row_index, datadir="./"):
    #Now process each object independently
    dl_data = True

    #Focus on one system only
    onerow = [input_table[row_index]]
    print(onerow)

    for obj in onerow: #input_table:

        #print(obj)

        #Turn ra in sexagesimal
        coord = SkyCoord(ra=obj['ra']*u.degree, dec=obj['dec']*u.degree, frame='icrs')
        rastring = coord.ra.to_string(unit=u.hour, sep='', pad=True, precision=2)
        destring = coord.dec.to_string(unit=u.degree, sep='', pad=True, precision=1)

        #Define Objname
        targname = 'MUSE_J{}{}_{}'.format(rastring, destring, obj['obsmode'])
        print('Target name: {}'.format(targname))

        #open a directory in datadir with targname
        outdir = datadir + targname
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        #Download the file listed under access_url
        xmllinksfile = outdir+'/'+targname+'_links.xml'
        if not(os.path.exists(xmllinksfile)):
          download_file(obj['access_url'], xmllinksfile)

        #Parse the file and save it as a fits table
        linkstable = parse_xml_links(xmllinksfile)
        linkstable.write(xmllinksfile.replace('.xml', '.fits'), overwrite=True)

        #Find whitelight image and download it
        whitelight_row = linkstable[linkstable['eso_category']=='ANCILLARY.IMAGE.WHITELIGHT']
        if len(whitelight_row)==1:
            whitelight_url = whitelight_row['access_url'][0]
            if not(os.path.exists(outdir+'/'+targname+'_whitelight.fits')):
                download_file(whitelight_url, outdir+'/'+targname+'_whitelight.fits')
        else:
            print('No whitelight image found for target {}'.format(targname))

        #Find PDF preview file and download it
        pdf_row = linkstable[linkstable['description'] == 'PDF preview of the requested file']
        if len(pdf_row)==1:
            pdf_url = pdf_row['access_url'][0]
            if not(os.path.exists(outdir+'/'+targname+'_preview.pdf')):
                download_file(pdf_url, outdir+'/'+targname+'_preview.pdf')
        else:
            print('No PDF preview image found for target {}'.format(targname))

        #Find QC preview image and download  it
        qc_row = linkstable[(linkstable['eso_category'] == 'ANCILLARY.PREVIEW') & (np.char.count(linkstable['eso_origfile'], 'dpc')>0)]
        if len(qc_row)==1:
            qc_url = qc_row['access_url'][0]
            if not(os.path.exists(outdir+'/'+targname+'_QCpreview.png')):
                download_file(qc_url, outdir+'/'+targname+'_QCpreview.png')
        else:
            print('No QC preview image found for target {}'.format(targname))

        if dl_data:
            #Find how many OBs are in this target
            ob_ids = np.unique((obj['obs_id'][1:-1]).replace('"','').split(','))
            print('Number of OBs in this target: {}'.format(len(ob_ids)))
            #Now run a query to get basic info about each science frame in each OB
            ob_info = query_ob_info(ob_ids)
            #Find unique nights
            unique_nights = np.unique(ob_info['night'])
            print('Number of unique nights in this target: {}'.format(len(unique_nights)))
            #Create a dictionary with numbers as key and night as value
            night_dict = {}
            for i, night in enumerate(unique_nights):
                night_dict[i+1] = night

            #Download each OB
            for i in range(1,len(unique_nights)+1):
                print("Running Night {}: {}".format(i, night_dict[i]))
                batch_outdir = outdir + '/Batch{}/Raw'.format(i)
                if not os.path.exists(batch_outdir):
                    os.makedirs(batch_outdir)
                download_batch(ob_info, night_dict[i], batch_outdir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download MUSE raw data for a given target.")
    parser.add_argument("input_file", help="Path to the input FITS table containing targets.")
    parser.add_argument("row_index", type=int, help="Index of the row to process in the input table (0-indexed).")
    parser.add_argument("--datadir", default="./", help="Directory where data will be downloaded.")
    args = parser.parse_args()

    # Read the input table
    input_table = Table.read(args.input_file)
    
    # Run main routine
    main(input_table, args.row_index, datadir=args.datadir)
