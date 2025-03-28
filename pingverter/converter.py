'''
Dependency of PINGMapper: https://github.com/CameronBodine/PINGMapper

Repository: https://github.com/CameronBodine/PINGVerter
PyPi: https://pypi.org/project/pingverter/ 

Developed by Cameron S. Bodine

###############
Acknowledgments
###############

None of this work would have been possible without the following repositories:

PyHum: https://github.com/dbuscombe-usgs/PyHum
SL3Reader: https://github.com/halmaia/SL3Reader
sonarlight: https://github.com/KennethTM/sonarlight


MIT License

Copyright (c) 2024 Cameron S. Bodine

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import os, sys
from pingverter import hum, low, cerul
import time
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
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

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

            auto_decode = False

            break

        else:
            auto_decode = True

    
    # Consider adding decodeHeadStruct() function back into PINGVerter....
    if not gotHeader:
        # sys.exit("\n#####\nERROR: Out of SON files... \n"+
        #         "Unable to determine header length.")
        print("\n#####\nERROR: Out of SON files... \n\n"+
                "Trying to automatically determine header length...")
        for beam, meta in beamMeta.items():
        
            # Get SON file
            son = meta['sonFile']

            # Autodecode
            headbytes = humminbird._decodeHeadStruct(son)

            if headbytes > 0:
                print("\n######################\nSLAMMA-JAMMA-DING-DONG \n"+
                      "Header Length Determined: {}\n\n".format(headbytes))
                print("As you were....\n\n")
                gotHeader = True

                # Add headbytes to humminbird object
                humminbird.frame_header_size = headbytes
                break
    if not gotHeader:
        sys.exit("\n#####\nERROR: Out of SON files... \n"+
                "Unable to automatically decode sonar header.")
        
    
    #############################################
    # Get the SON header structure and attributes
    #############################################
    if not auto_decode:
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

def low2pingmapper(input: str, out_dir: str, nchunk: int=500, tempC: float=10, exportUnknown: bool=False):

    # Make sure input exists
    assert os.path.isfile(input), "{} does not exist.".format(input)

    # Create the class
    lowrance = low(inFile=input, nchunk=nchunk, exportUnknown=exportUnknown)

    # Store temperature
    lowrance.tempC = float(tempC)/10

    ######################
    # Decode Lowrance File
    ######################

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Create 'meta' directory if it doesn't exist
    metaDir = os.path.join(out_dir, 'meta')
    try:
        os.mkdir(metaDir)
    except:
        pass
    lowrance.metaDir = metaDir # Store metadata directory

    # Get Lowrance file length
    lowrance._getFileLen()

    # Parse file header ***Probably not needed***
    lowrance._parseFileHeader()

    # Parse ping headers (attributes) and do conversions
    lowrance._parsePingHeader()

    # Remove unknown beams
    lowrance._removeUnknownBeams()

    # Drop Beam 0 (83kHz) or 1 (200kHz) if necessary
    lowrance._removeDownBeams()

    # Split sidescan, if necessary
    beams = lowrance.header_dat['beam'].unique()
    if 5 in beams:
        lowrance._splitLowSS()
        flip_port = True
    else:
        flip_port = False

    # Recalculate record number
    lowrance._recalcRecordNum()

    # Drop unknown
    if not exportUnknown:
        cols = lowrance.header_dat.columns
        cols = [c for c in cols if 'unknown' in c]
        
        lowrance.header_dat.drop(columns=cols, inplace=True)

    # Save ping metadata to csv based on beam
    lowrance._splitBeamsToCSV()

    # Store headBytes
    lowrance.headBytes = lowrance.frame_header_size

    # Not Humminbird Onix
    lowrance.isOnix = 0
    
    return lowrance


# =========================================================
# Lowrance to PINGMapper
# =========================================================

def cerul2pingmapper(input: str, out_dir: str, nchunk: int=500, tempC: float=10, exportUnknown: bool=False):
    '''
    '''
    # Make sure input exists
    assert os.path.isfile(input), "{} does not exist.".format(input)

    # Create the class
    cerulean = cerul(svlog = input, nchunk=nchunk, exportUnknown=exportUnknown)

    # Store Temperature
    cerulean.tempC = float(tempC)/10

    ######################
    # Decode Cerulean File
    ######################

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Create 'meta' directory if it doesn't exist
    metaDir = os.path.join(out_dir, 'meta')
    try:
        os.mkdir(metaDir)
    except:
        pass
    cerulean.metaDir = metaDir # Store metadata directory

    # Get Cerulean file length
    cerulean._getFileLen()

    # Parse the file header
    cerulean._parseFileHeader()

    # Parse all packet headers
    if exportUnknown:
        cerulean._locatePacketsRaw()

    # Locate Packet Headers
    cerulean._locatePackets()

    # Set beam
    cerulean._convertBeam()

    # Set frequency
    cerulean._convertFrequency()
    
    # Recalculate record num
    cerulean._recalcRecordNum()

    # Save to file
    cerulean._splitBeamsToCSV()

    # print(cerulean)

    return cerulean


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



