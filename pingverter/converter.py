
import os, sys
from pingverter import hum, low
import time, datetime
import pandas as pd
from glob import glob

from joblib import Parallel, delayed, cpu_count

# =========================================================
# Humminbird to PINGMapper
# =========================================================

def hum2pingmapper(input: str, out_dir: str, nchunk: int=500, tempC: float=10, exportUnknown: bool=False):
    '''
    
    '''
    # Make sure input exists
    assert os.path.isfile(input), "{} does not exist.".format(input)

    # Create the class
    humminbird = hum(humFile=input, nchunk=nchunk, exportUnknown=exportUnknown)

    # Store temperatue
    humminbird.tempC = float(tempC)/10

    #################
    # Decode DAT File
    #################

    start_time = time.time()
    print("\nGetting DAT Metadata...")
    print(input)
    humminbird._getHumDatStruct()

    # Read in the humdat data
    if humminbird.isOnix == 0:
        humminbird._getHumdat()

        # Determine epsg code and transformation (if we can, ONIX doesn't have
        ## lat/lon in DAT, so will determine at a later processing step).
        humminbird._getEPSG()
    else:
        humminbird._decodeOnix()

    # Create 'meta' directory if it doesn't exist
    metaDir = os.path.join(out_dir, 'meta')
    try:
        os.mkdir(metaDir)
    except:
        pass
    humminbird.metaDir = metaDir #Store metadata directory in sonObj

    # Save DAT metadata to file (csv)
    outFile = os.path.join(metaDir, 'DAT_meta.csv') # Specify file directory & name
    pd.DataFrame.from_dict(humminbird.humDat, orient='index').T.to_csv(outFile, index=False) # Export DAT df to csv
    humminbird.datMetaFile = outFile # Store metadata file path in sonObj
    del outFile

    print("\nDone!")
    print("Time (s):", round(time.time() - start_time, ndigits=1))

    #####################################
    # Generate sonar file meta/attributes
    #####################################

    # Find available SON files
    sonFiles = glob(os.path.join(humminbird.sonFile, '*.SON'))
    
    # Determine which sonar beams are present (B000,B001,..)
    chanAvail = {}
    for s in sonFiles:
        beam = os.path.split(s)[-1].split('.')[0] #Get beam number (B000,B001,..)
        chanAvail[beam] = s
    del s, beam

    humminbird.chanAvail = chanAvail

    # Store sonar file meta in humminbird object
    ## Dictionary of dictionaries
    beamMeta = humminbird.beamMeta = {}
    for chan, file in chanAvail.items():
        
        beamMeta[chan] = {}

        # Get beam name
        beamMeta[chan]['beamName'] = humminbird._getBeamName(chan)

        beamMeta[chan]['beam'] = chan
        beamMeta[chan]['sonFile'] = file

        # Output csv name
        csv = '{}_{}_{}'.format(chan, beamMeta[chan]['beamName'], 'meta.csv')
        beamMeta[chan]['metaCSV'] = os.path.join(metaDir, csv)


    ##############################################################
    # Determine ping header structure (varies by Humminbird model)
    ##############################################################
    start_time = time.time()
    print("\nGetting Header Structure...")
    

    gotHeader = False # Flag indicating if length of header is found

    for beam, meta in beamMeta.items():
        
        # Get SON file
        son = meta['sonFile']

        # Count headbytes
        headbytes = humminbird._cntHead(son)

        if headbytes > 0: # Header length found
            print("Header Length: {}".format(headbytes))
            gotHeader = True

            # Add headbytes to humminbird object
            humminbird.frame_header_size = headbytes

            break

    
    # Consider adding decodeHeadStruct() function back into PINGVerter....


    if not gotHeader:
        sys.exit("\n#####\nERROR: Out of SON files... \n"+
                "Unable to determine header length.")
        
    
    #############################################
    # Get the SON header structure and attributes
    #############################################

    humminbird._getHeadStruct()

    ##################
    # Parse son header
    ##################

    _ = Parallel(n_jobs = len(beamMeta), verbose=10 )(delayed(humminbird._parsePingHeader)(meta['sonFile'], meta['metaCSV']) for beam, meta in beamMeta.items())

    print("\nDone!")
    print("Time (s):", round(time.time() - start_time, ndigits=1))
    
    return humminbird


# =========================================================
# Lowrance to PINGMapper
# =========================================================

def low2pingmapper(low_file: str, out_dir: str):
    
    return


# =========================================================
# Lowrance to Humminbird
# =========================================================

def low2hum(input: str, output: str):

    # Make sure input exists
    assert os.path.isfile(input), "{} does not exist.".format(input)

    # Create the classes
    lowrance = low(input)
    humminbird = hum(output)

    # Make output files
    humminbird._makeOutFiles()

    # Start the decode
    start_time = time.time()
    print('\n\nDecoding Lowrance File...')

    # Get Lowrance file length
    lowrance._getFileLen()

    # Parse file header
    lowrance._parseFileHeader()

    # Parse ping headers (attributes)
    lowrance._parsePingHeader()

    # Convert ping attributes to known units
    lowrance._convertPingAttributes()

    # lowrance.header_dat.to_csv('lowrance_test.csv')

    print("Time (s):", round(time.time() - start_time, ndigits=1))

    #######################
    # Convert to Humminbird
    #######################

    start_time = time.time()
    print('\n\nConverting to Humminbird...')

    # Convert to Humminbird attributes
    humminbird._convertLowHeader(lowrance)

    # Drop unknown beams
    humminbird._removeUnknownBeams()

    # Split sidescan, if necessary
    beams = humminbird.header_dat['beam'].unique()
    if 5 in beams:
        humminbird._splitLowSS()
        flip_port = True
    else:
        flip_port = False

    # Recalculate record number
    humminbird._recalcRecordNum()

    # humminbird.header_dat.to_csv('hum_converted.csv')

    # Get necessary data for DAT file
    humminbird._convertLowDAT(lowrance)

    # Write DAT to file
    humminbird._writeDAT()

    # Save b001 and b002
    beams = humminbird.header_dat['beam'].unique()

    if 0 in beams:
        humminbird._writeSonfromLow(0, lowrance.frame_header_size, lowrance.path)

    if 1 in beams:
        humminbird._writeSonfromLow(1, lowrance.frame_header_size, lowrance.path)

    if 2 in beams:
        humminbird._writeSonfromLow(2, lowrance.frame_header_size, lowrance.path, flip_port)

    if 3 in beams:
        humminbird._writeSonfromLow(3, lowrance.frame_header_size, lowrance.path)

    if 4 in beams:
        humminbird._writeSonfromLow(4, lowrance.frame_header_size, lowrance.path)

    # Split b005 (lowrance sidescan) into port (2) and star (3)


    humminbird.header_dat.to_csv('hum_converted.csv')

    print("Time (s):", round(time.time() - start_time, ndigits=1))

    return


