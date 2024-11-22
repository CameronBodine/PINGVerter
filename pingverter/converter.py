
import os, sys
from pingverter import hum, low
import time

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



