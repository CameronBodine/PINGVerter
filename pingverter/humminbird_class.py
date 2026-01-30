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


import os, sys, struct
import numpy as np
import pandas as pd
from array import array as arr
import pyproj
import datetime

from .lowrance_class import low

class hum(object):

    #===========================================================================
    def __init__(self, humFile: str, nchunk: int=0, exportUnknown: bool=False):
        
        self.humFile = humFile
        self.sonFile = humFile.split('.DAT')[0]
        self.nchunk = nchunk
        self.exportUnknown = exportUnknown

        self.head_start_val = 3235818273
        self.head_end_val = 33

        self.son8bit = True

        return
    #===========================================================================
    # Humminbird to PINGMapper
    #===========================================================================

    def _fread_dat(self,
                   infile: str,
                   num: int,
                   typ: str):
        '''
        Helper function that reads binary data in a file.

        ----------------------------
        Required Pre-processing step
        ----------------------------
        Called from self._getHumDat(), self._cntHead(), self._decodeHeadStruct(),
        self._getSonMeta(), self._loadSonChunk()

        ----------
        Parameters
        ----------
        infile : file
            DESCRIPTION - A binary file opened in read mode at a pre-specified
                          location.
        num : int
            DESCRIPTION - Number of bytes to read.
        typ : type
            DESCRIPTION - Byte type

        -------
        Returns
        -------
        List of decoded binary data

        --------------------
        Next Processing Step
        --------------------
        Returns list to function it was called from.
        '''
        dat = arr(typ)
        dat.fromfile(infile, num)
        return(list(dat))

    def _getHumDatStruct(self, ):
        '''
        Determines .DAT file structure for a sonObj instance.

        ----------------------------
        Required Pre-processing step
        ----------------------------
        self.__init__()

        -------
        Returns
        -------
        A dictionary with .DAT file structure stored in self.humDatStruct with
        the following format:

        self.humDatStruct = {name : [byteIndex, offset, dataLen, data],
                             .... : ....} where:
            name == Name of attribute;
            byteIndex == Index indicating position of name;
            offset == Byte offset for the actual data;
            dataLen == number of bytes for data (i.e. utm_e is 4 bytes long);
            data = actual value of the attribute.

        --------------------
        Next Processing Step
        --------------------
        self._getHumDat()
        '''

        # Get class attributes as variables
        humFile = self.humFile
        nchunk = self.nchunk

        # DAT structure dependent on file size
        self.datLen = datLen = os.path.getsize(humFile)
        
        ############################
        # Set DAT struct from datLen
        ############################

        # The length (in bytes, `datLen`) of the .DAT file can indicate which
        ## Humminbird model a .DAT file is from.  That means we know where certain
        ## attributes are stored in the .DAT file.
        #1199, Helix
        if datLen == 64:
            self.isOnix = 0

            humDic = {
            'endianness':'>i', #>=big endian; I=unsigned Int
            'SP1':[0, 0, 1, -1], #Unknown (spacer)
            'water_code':[1, 0, 1, -1], #Water code: 0=fresh,1=deep salt, 2=shallow salt
            'SP2':[2, 0, 1, -1], #Unknown (spacer)
            'unknown_1':[3, 0, 1, -1], #Unknown (gps flag?)
            'sonar_name':[4, 0, 4, -1], #Sonar name
            'unknown_2':[8, 0, 4, -1], #Unknown
            'unknown_3':[12, 0, 4, -1], #Unknown
            'unknown_4':[16, 0, 4, -1], #Unknown
            'unix_time':[20, 0, 4, -1], #Unix Time
            'utm_e':[24, 0, 4, -1], #UTM X
            'utm_n':[28, 0, 4, -1], #UTM Y
            'filename':[32, 0, 10, -1], #Recording name
            'unknown_5':[42, 0, 2, -1], #Unknown
            'numrecords':[44, 0, 4, -1], #Number of records
            'recordlens_ms':[48, 0, 4, -1], #Recording length milliseconds
            'linesize':[52, 0, 4, -1], #Line Size (?)
            'unknown_6':[56, 0, 4, -1], #Unknown
            'unknown_7':[60, 0, 4, -1], #Unknown
            }

        #Solix (Little Endian)
        elif datLen == 96:
            self.isOnix = 0
            humDic = {
            'endianness':'<i', #<=little endian; I=unsigned Int
            'SP1':[0, 0, 1, -1], #Unknown (spacer)
            'water_code':[1, 0, 1, -1], #Need to check if consistent with other models (1=fresh?)
            'SP2':[2, 0, 1, -1], #Unknown (spacer)
            'unknown_1':[3, 0, 1, -1], #Unknown (gps flag?)
            'sonar_name':[4, 0, 4, -1], #Sonar name
            'unknown_2':[8, 0, 4, -1], #Unknown
            'unknown_3':[12, 0, 4, -1], #Unknown
            'unknown_4':[16, 0, 4, -1], #Unknown
            'unix_time':[20, 0, 4, -1], #Unix Time
            'utm_e':[24, 0, 4, -1], #UTM X
            'utm_n':[28, 0, 4, -1], #UTM Y
            'filename':[32, 0, 12, -1], #Recording name
            'numrecords':[44, 0, 4, -1], #Number of records
            'recordlens_ms':[48, 0, 4, -1], #Recording length milliseconds
            'linesize':[52, 0, 4, -1], #Line Size (?)
            'unknown_5':[56, 0, 4, -1], #Unknown
            'unknown_6':[60, 0, 4, -1], #Unknown
            'unknown_7':[64, 0, 4, -1], #Unknown
            'unknown_8':[68, 0, 4, -1], #Unknown
            'unknown_9':[72, 0, 4, -1], #Unknown
            'unknown_10':[76, 0, 4, -1], #Unknown
            'unknown_11':[80, 0, 4, -1], #Unknown
            'unknown_12':[84, 0, 4, -1], #Unknown
            'unknown_13':[88, 0, 4, -1], #Unknown
            'unknown_14':[92, 0, 4, -1]
            }

        #### TESTING ######
        elif datLen == 100:
            self.isOnix = 0
            humDic = {
            'endianness':'<i', #<=little endian; I=unsigned Int
            'SP1':[0, 0, 1, -1], #Unknown (spacer)
            'water_code':[1, 0, 1, -1], #Need to check if consistent with other models (1=fresh?)
            'SP2':[2, 0, 1, -1], #Unknown (spacer)
            'unknown_1':[3, 0, 1, -1], #Unknown (gps flag?)
            'sonar_name':[4, 0, 4, -1], #Sonar name
            'unknown_2':[8, 0, 4, -1], #Unknown
            'unknown_3':[12, 0, 4, -1], #Unknown
            'unknown_4':[16, 0, 4, -1], #Unknown
            'unix_time':[20, 0, 4, -1], #Unix Time
            'utm_e':[24, 0, 4, -1], #UTM X
            'utm_n':[28, 0, 4, -1], #UTM Y
            'filename':[32, 0, 12, -1], #Recording name
            'numrecords':[44, 0, 4, -1], #Number of records
            'recordlens_ms':[48, 0, 4, -1], #Recording length milliseconds
            'linesize':[52, 0, 4, -1], #Line Size (?)
            'unknown_5':[56, 0, 4, -1], #Unknown
            'unknown_6':[60, 0, 4, -1], #Unknown
            'unknown_7':[64, 0, 4, -1], #Unknown
            'unknown_8':[68, 0, 4, -1], #Unknown
            'unknown_9':[72, 0, 4, -1], #Unknown
            'unknown_10':[76, 0, 4, -1], #Unknown
            'unknown_11':[80, 0, 4, -1], #Unknown
            'unknown_12':[84, 0, 4, -1], #Unknown
            'unknown_13':[88, 0, 4, -1], #Unknown
            'unknown_14':[92, 0, 4, -1], #Unknown
            'unknown_15':[96, 0, 4, -1]
            }

        #Onix
        else:
            humDic = {}
            self.isOnix = 1

        self.humDatStruct = humDic
        return
    
    def _getHumdat(self):

        '''
        Decode .DAT file using known DAT file structure.

        ----------------------------
        Required Pre-processing step
        ----------------------------
        self._getHumDatStruct()

        -------
        Returns
        -------
        A dictionary stored in self.humDat containing data from .DAT file.

        --------------------
        Next Processing Step
        --------------------
        self._getEPSG()
        '''

        # Get necessary class attributes
        humDic = self.humDatStruct # Dictionary to store .DAT file structure.
        humFile = self.humFile # Path to .DAT file.
        datLen = self.datLen # Number of bytes in DAT file.
        t = self.tempC # Water temperature (Celcius) during survey divided by 10.

        humDat = {} #defaultdict(dict) # Empty dict to store DAT contents
        endian = humDic['endianness'] # Is data big or little endian?
        file = open(humFile, 'rb') # Open the file

        # Search for humDic items in DAT file
        for key, val in humDic.items(): # Iterate each item in humDic
            if key == 'endianness':
                pass
            else:
                file.seek(val[0]) # Move to correct byte offset
                # If the expected data is 4 bytes long
                if val[2] == 4:
                    byte = struct.unpack(endian, arr('B', self._fread_dat(file, val[2], 'B')).tobytes())[0] # Decode byte
                # If the expected data is less than 4 bytes long
                elif val[2] < 4:
                    byte = self._fread_dat(file, val[2], 'B')[0] # Decode byte
                # If the expected data is greater than 4 bytes long
                elif val[2] > 4:
                    byte = arr('B', self._fread_dat(file, val[2], 'B')).tobytes().decode() # Decode byte
                # Something went wrong...
                else:
                    byte = -9999
                humDat[key] = byte # Store the data

        file.close() # Close the file

        # Determine Humminbird water type setting and update (S)alinity appropriately
        waterCode = humDat['water_code']
        if datLen == 64:
            if waterCode == 0:
                humDat['water_type'] = 'fresh'
                S = 1
            elif waterCode == 1:
                humDat['water_type'] = 'deep salt'
                S = 35
            elif waterCode == 2:
                humDat['water_type'] = 'shallow salt'
                S = 30
            else:
                humDat['water_type'] = 'unknown'
        # #Need to figure out water code for solix
        # elif datLen == 96:
        #     if waterCode == 1:
        #         humDat['water_type'] = 'fresh'
        #         S = 1
        #     else:
        #         humDat['water_type'] = 'unknown'
        #         c = 1475

        # ###### TESTING ######
        # elif datLen == 100:
        #     if waterCode == 1:
        #         humDat['water_type'] = 'fresh'
        #         S = 1
        #     else:
        #         humDat['water_type'] = 'unknown'
        #         c = 1475

        elif datLen >= 96:
            if waterCode == 1:
                humDat['water_type'] = 'fresh'
                S = 1
            elif waterCode == 2:
                humDat['water_type'] = 'shallow salt'
                S = 30
            elif waterCode == 3:
                humDat['water_type'] = 'deep salt'
                S = 35
            else:
                humDat['water_type'] = 'unknown'

        # Calculate speed of sound based on temp & salinity
        c = 1449.05 + 45.7*t - 5.21*t**2 + 0.23*t**3 + (1.333 - 0.126*t + 0.009*t**2)*(S - 35)

        # Calculate time varying gain
        self.tvg = ((8.5*10**-5)+(3/76923)+((8.5*10**-5)/4))*c
        self.c = c
        self.S = S

        humDat['nchunk'] = self.nchunk
        self.humDat = humDat # Store data in class attribute for later use

        # Calculate pixel size in [m]
        t = 0.108 # Transducer length
        f = 455
        # theta at 3dB in the horizontal
        theta3dB = np.arcsin(c/(t*(f*1000)))
        #resolution of 1 sidescan pixel to nadir
        ft = (np.pi/2)*(1/theta3dB)
        # size of pixel in meters
        pix_m = (1/ft)

        humDat['pixM'] = pix_m
        self.pixM = pix_m

        return
    
    def _getEPSG(self, utm_e: int=None, utm_n: int=None):
        '''
        Determines appropriate UTM zone based on location (EPSG 3395 Easting/Northing)
        provided in .DAT file.  This is used to project coordinates from
        EPSG 3395 to local UTM zone.

        ----------------------------
        Required Pre-processing step
        ----------------------------
        self._getHumdat()

        -------
        Returns
        -------
        self.trans which will re-poroject Humminbird easting/northings to local UTM zone.

        --------------------
        Next Processing Step
        --------------------
        self._cntHead()
        '''

        # Convert easting/northing to latitude/longitude
        lat = np.arctan(np.tan(np.arctan(np.exp(utm_n/ 6378388.0)) * 2.0 - 1.570796326794897) * 1.0067642927) * 57.295779513082302
        lon = (utm_e * 57.295779513082302) / 6378388.0

        # Determine epsg code
        self.humDat['epsg'] = "epsg:"+str(int(float(self._convert_wgs_to_utm(lon, lat))))
        self.humDat['wgs'] = "epsg:4326"

        # Configure re-projection function
        self.trans = pyproj.Proj(self.humDat['epsg'])

        return
    
    def _convert_wgs_to_utm(self, lon: float, lat: float):
        """
        This function estimates UTM zone from geographic coordinates
        see https://stackoverflow.com/questions/40132542/get-a-cartesian-projection-accurate-around-a-lat-lng-pair
        """
        utm_band = str((np.floor((lon + 180) / 6 ) % 60) + 1)
        if len(utm_band) == 1:
            utm_band = '0'+utm_band
        if lat >= 0:
            epsg_code = '326' + utm_band
        else:
            epsg_code = '327' + utm_band
        return epsg_code

    def _decodeOnix(self):
        '''
        Decodes .DAT file from Onix Humminbird models.  Onix has a significantly
        different .DAT file structure compared to other Humminbird models,
        requiring a specific function to decode the file.

        -------
        Returns
        -------
        A dictionary stored in self.humDat containing data from .DAT file.
        '''
        t = self.tempC # Water temperature (Celcius) during survey divided by 10.

        fid2 = open(self.humFile, 'rb') # Open file

        dumpstr = fid2.read() # Store file contents
        fid2.close() # Close the file

        if sys.version.startswith('3'):
          dumpstr = ''.join(map(chr, dumpstr))

        humdat = {}
        hd = dumpstr.split('<')[0]
        tmp = ''.join(dumpstr.split('<')[1:])
        humdat['NumberOfPings'] = int(tmp.split('NumberOfPings=')[1].split('>')[0])
        humdat['TotalTimeMs'] = int(tmp.split('TotalTimeMs=')[1].split('>')[0])
        humdat['linesize'] = int(tmp.split('PingSizeBytes=')[1].split('>')[0])
        humdat['FirstPingPeriodMs'] = int(tmp.split('FirstPingPeriodMs=')[1].split('>')[0])
        humdat['BeamMask'] = int(tmp.split('BeamMask=')[1].split('>')[0])
        humdat['Chirp1StartFrequency'] = int(tmp.split('Chirp1StartFrequency=')[1].split('>')[0])
        humdat['Chirp1EndFrequency'] = int(tmp.split('Chirp1EndFrequency=')[1].split('>')[0])
        humdat['Chirp2StartFrequency'] = int(tmp.split('Chirp2StartFrequency=')[1].split('>')[0])
        humdat['Chirp2EndFrequency'] = int(tmp.split('Chirp2EndFrequency=')[1].split('>')[0])
        humdat['Chirp3StartFrequency'] = int(tmp.split('Chirp3StartFrequency=')[1].split('>')[0])
        humdat['Chirp3EndFrequency'] = int(tmp.split('Chirp3EndFrequency=')[1].split('>')[0])
        humdat['SourceDeviceModelId2D'] = int(tmp.split('SourceDeviceModelId2D=')[1].split('>')[0])
        humdat['SourceDeviceModelIdSI'] = int(tmp.split('SourceDeviceModelIdSI=')[1].split('>')[0])
        humdat['SourceDeviceModelIdDI'] = int(tmp.split('SourceDeviceModelIdDI=')[1].split('>')[0])
        humdat['water_type'] = 'fresh' #'shallow salt' #'deep salt'
        self.humDat = humdat # Store data in class attribute for later use

        # Unsure about salinity for Onix
        S = 1

        # Calculate speed of sound based on temp & salinity
        c = 1449.05 + 45.7*t - 5.21*t**2 + 0.23*t**3 + (1.333 - 0.126*t + 0.009*t**2)*(S - 35)

        # Calculate time varying gain
        self.tvg = ((8.5*10**-5)+(3/76923)+((8.5*10**-5)/4))*c
        self.c = c
        self.S = S

        self.humDat['nchunk'] = self.nchunk
        self.humDat = self.humDat # Store data in class attribute for later use

        # Calculate pixel size in [m]
        t = 0.108 # Transducer length
        f = 455
        # theta at 3dB in the horizontal
        theta3dB = np.arcsin(c/(t*(f*1000)))
        #resolution of 1 sidescan pixel to nadir
        ft = (np.pi/2)*(1/theta3dB)
        # size of pixel in meters
        pix_m = (1/ft)

        self.humDat['pixM'] = pix_m
        self.pixM = pix_m

        return

    def _getBeamName(self, beam: str):

        '''
        '''

        if beam == 'B000':
            beamName = 'ds_lowfreq'
        elif beam == 'B001':
            beamName = 'ds_highfreq'
        elif beam == 'B002':
            beamName = 'ss_port'
        elif beam == 'B003':
            beamName = 'ss_star'
        elif beam == 'B004':
            beamName = 'ds_vhighfreq'
        else:
            beamName = 'unknown'
        return beamName
    
    def _cntHead(self, sonFile: str):

        '''
        Determine .SON ping header length based on known Humminbird
        .SON file structure.  Humminbird stores sonar records in packets, where
        the first x bytes of the packet contain metadata (record number, northing,
        easting, time elapsed, depth, etc.), proceeded by the sonar/ping returns
        associated with that ping.  This function will search the first
        ping to determine the length of ping header.

        ----------------------------
        Required Pre-processing step
        ----------------------------
        self.__init__()

        -------
        Returns
        -------
        self.headBytes, indicating length, in bytes, of ping header.

        --------------------
        Next Processing Step
        --------------------
        self._getHeadStruct()
        '''
        
        # file = open(sonFile, 'rb') # Open sonar file
        with open(sonFile, 'rb') as file:
            i = 0 # Counter to track sonar header length
            foundEnd = False # Flag to track if end of sonar header found
            while foundEnd is False and i < 200:
                lastPos = file.tell() # Get current position in file (byte offset from beginning)
                byte = self._fread_dat(file, 1, 'B') # Decode byte value

                # Check if we found the end of the ping.
                ## A value of 33 **may** indicate end of ping.
                ## 67 is the shortest known Humminbird ping header length.
                if byte[0] == self.head_end_val and lastPos > 65:
                    # Double check we found the actual end by moving backward -6 bytes
                    ## to see if value is 160 (spacer preceding number of ping records)
                    file.seek(-6, 1)
                    byte = self._fread_dat(file, 1, 'B')
                    if byte[0] == 160:
                        foundEnd = True
                    else:
                        # Didn't find the end of header
                        # Move cursor back to lastPos+1
                        file.seek(lastPos+1)
                else:
                    # Haven't found the end
                    pass
                i+=1

        # i reaches 200, then we have exceeded known Humminbird header length.
        ## Set i to 0, then the next sonFile will be checked.
        if i == 200:
            i = 0

        file.close()
        self.headBytes = i # Store data in class attribute for later use
        return i
    
    def _getHeadStruct(self):
        '''
        '''

        # Get frame header size
        header_len = self.frame_header_size

        if header_len == 67:
            headStruct = np.dtype([
                ('head_start', '>u4'),
                ('SP128', '>u1'),
                ('record_num', '>u4'),
                ('SP129', '>u1'),
                ('time_s', '>u4'),
                ('SP130', '>u1'),
                ('utm_e', '>i4'),
                ('SP131', '>u1'),
                ('utm_n', '>i4'),
                ('SP132', '>u1'),
                ('gps1', '>u2'),
                ('instr_heading', '>u2'),
                ('SP133', '>u1'),
                ('gps2', '>u2'),
                ('speed_ms', '>u2'),
                ('SP135', '>u1'),
                ('inst_dep_m', '>u4'),
                ('SP80', '>u1'),
                ('beam', '>u1'),
                ('SP81', '>u1'),
                ('volt_scale', '>u1'),
                ('SP146', '>u1'),
                ('f', '>u4'),
                ('SP83', '>u1'),
                ('unknown_83', '>u1'),
                ('SP84', '>u1'),
                ('unknown_84', '>u1'),
                ('SP149', '>u1'),
                ('unknown_149', '>u4'),
                ('SP86', '>u1'),
                ('e_err_m', '>u1'),
                ('SP87', '>u1'),
                ('n_err_m', '>u1'),
                ('SP160', '>u1'),
                ('ping_cnt', '>u4'),
                ('head_end', '>u1')
            ])

        elif header_len == 72:
            headStruct = np.dtype([
                ('head_start', '>u4'),
                ('SP128', '>u1'),
                ('record_num', '>u4'),
                ('SP129', '>u1'),
                ('time_s', '>u4'),
                ('SP130', '>u1'),
                ('utm_e', '>i4'),
                ('SP131', '>u1'),
                ('utm_n', '>i4'),
                ('SP132', '>u1'),
                ('gps1', '>u2'),
                ('instr_heading', '>u2'),
                ('SP133', '>u1'),
                ('gps2', '>u2'),
                ('speed_ms', '>u2'),
                ('SP134', '>u1'),
                ('unknown_134', '>u4'),
                ('SP135', '>u1'),
                ('inst_dep_m', '>u4'),
                ('SP80', '>u1'),
                ('beam', '>u1'),
                ('SP81', '>u1'),
                ('volt_scale', '>u1'),
                ('SP146', '>u1'),
                ('f', '>u4'),
                ('SP83', '>u1'),
                ('unknown_83', '>u1'),
                ('SP84', '>u1'),
                ('unknown_84', '>u1'),
                ('SP149', '>u1'),
                ('unknown_149', '>u4'),
                ('SP86', '>u1'),
                ('e_err_m', '>u1'),
                ('SP87', '>u1'),
                ('n_err_m', '>u1'),
                ('SP160', '>u1'),
                ('ping_cnt', '>u4'),
                ('head_end', '>u1')
            ])

        elif header_len == 152:
            headStruct = np.dtype([
                ('head_start', '>u4'),
                ('SP128', '>u1'),
                ('record_num', '>u4'),
                ('SP129', '>u1'),
                ('time_s', '>u4'),
                ('SP130', '>u1'),
                ('utm_e', '>i4'),
                ('SP131', '>u1'),
                ('utm_n', '>i4'),
                ('SP132', '>u1'),
                ('gps1', '>u2'),
                ('instr_heading', '>u2'),
                ('SP133', '>u1'),
                ('gps2', '>u2'),
                ('speed_ms', '>u2'),
                ('SP134', '>u1'),
                ('unknown_134', '>u4'),
                ('SP135', '>u1'),
                ('inst_dep_m', '>u4'),
                ('SP136', '>u1'),
                ('unknown_136', '>u4'),
                ('SP137', '>u1'),
                ('unknown_137', '>u4'),
                ('SP138', '>u1'),
                ('unknown_138', '>u4'),
                ('SP139', '>u1'),
                ('unknown_139', '>u4'),
                ('SP140', '>u1'),
                ('unknown_140', '>u4'),
                ('SP141', '>u1'),
                ('unknown_141', '>u4'),
                ('SP142', '>u1'),
                ('unknown_142', '>u4'),
                ('SP143', '>u1'),
                ('unknown_143', '>u4'),
                ('SP80', '>u1'),
                ('beam', '>u1'),
                ('SP81', '>u1'),
                ('volt_scale', '>u1'),
                ('SP146', '>u1'),
                ('f', '>u4'),
                ('SP83', '>u1'),
                ('unknown_83', '>u1'),
                ('SP84', '>u1'),
                ('unknown_84', '>u1'),
                ('SP149', '>u1'),
                ('unknown_149', '>u4'),
                ('SP86', '>u1'),
                ('e_err_m', '>u1'),
                ('SP87', '>u1'),
                ('n_err_m', '>u1'),
                ('SP152', '>u1'),
                ('unknown_152', '>u4'),
                ('SP153', '>u1'),
                ('unknown_153', '>u4'),
                ('SP154', '>u1'),
                ('unknown_154', '>u4'),
                ('SP155', '>u1'),
                ('unknown_155', '>u4'),
                ('SP156', '>u1'),
                ('unknown_156', '>u4'),
                ('SP157', '>u1'),
                ('unknown_157', '>u4'),
                ('SP158', '>u1'),
                ('unknown_158', '>u4'),
                ('SP159', '>u1'),
                ('unknown_159', '>u4'),
                ('SP160', '>u1'),
                ('ping_cnt', '>u4'),
                ('head_end', '>u1')
            ])

        self.son_struct = headStruct

        return

    def _decodeHeadStruct(self, sonFile: str):
        '''
        This function attempts to automatically decode the sonar return header
        structure if self.headValid == FALSE as determined by self._checkHeadStruct().
        This function will iterate through each byte at the beginning of the
        sonar file, decode the byte, determine if it matches any known or unknown
        spacer value (ping attribute 'name') and if it does, store the
        byte offset.

        ----------
        Parameters
        ----------
        exportUnknown : bool
            DESCRIPTION - Flag indicating if unknown attributes in ping
                          should be exported or not.  If a user of PING Mapper
                          determines what an unkown attribute actually is, please
                          report using a github issue.

        ----------------------------
        Required Pre-processing step
        ----------------------------
        self._cntHead()

        -------
        Returns
        -------
        A dictionary with .SON file structure stored in self.headStruct with
        the following format:

        self.headStruct = {byteVal : [byteIndex, offset, dataLen, name],
                           ...} where:
            byteVal == Spacer value (integer) preceding attribute values (i.e. depth);
            byteIndex == Index indicating position of byteVal;
            offset == Byte offset for the actual data;
            dataLen == number of bytes for data (i.e. utm_e is 4 bytes long);
            name = name of attribute.

        --------------------
        Next Processing Step
        --------------------
        self._checkHeadStruct()
        '''

        header_len = 0
        headBytes = 200
        headStruct = []
        toCheck = {
            128:[('SP128', '>u1'), ('record_num', '>u4')], #Record Number (Unique for each ping)
            129:[('SP129', '>u1'), ('time_s', '>u4')], #Time Elapsed milliseconds
            130:[('SP130', '>u1'), ('utm_e', '>i4')], #UTM X
            131:[('SP131', '>u1'), ('utm_n', '>i4'),], #UTM Y
            132:[('SP132', '>u1'), ('gps1', '>u2'), ('instr_heading', '>u2')], #GPS quality flag (?) and heading
            133:[('SP133', '>u1'), ('gps2', '>u2'), ('speed_ms', '>u2'),], #GPS quality flag (?) & speed in meters/second
            134:[('SP134', '>u1'), ('unknown_134', '>u4'),], #Unknown
            135:[('SP135', '>u1'), ('inst_dep_m', '>u4'),], #Depth in centimeters, then converted to meters
            136:[('SP136', '>u1'), ('unknown_136', '>u4'),], #Unknown
            137:[('SP137', '>u1'), ('unknown_137', '>u4')], #Unknown
            138:[('SP138', '>u1'), ('unknown_138', '>u4'),], #Unknown
            139:[('SP139', '>u1'), ('unknown_139', '>u4'),], #Unkown
            140:[('SP140', '>u1'), ('unknown_140', '>u4'),], #Unknown
            141:[('SP141', '>u1'), ('unknown_141', '>u4'),], #Unknown
            142:[('SP142', '>u1'), ('unknown_142', '>u4'),], #Unknown
            143:[('SP143', '>u1'), ('unknown_143', '>u4'),], #Unknown
            80:[('SP80', '>u1'), ('beam', '>u1'),], #Beam number: 0 (50 or 83 kHz), 1 (200 kHz), 2 (SI Poort), 3 (SI Starboard)
            81:[('SP81', '>u1'), ('volt_scale', '>u1'),], #Volt Scale (?)
            146:[('SP146', '>u1'), ('f', '>u4'),], #Frequency of beam in hertz
            83:[('SP83', '>u1'), ('unknown_83', '>u1'),], #Unknown (number of satellites???)
            84:[('SP84', '>u1'), ('unknown_84', '>u1'),], #Unknown
            149:[('SP149', '>u1'), ('unknown_149', '>u4'),], #Unknown (magnetic deviation???)
            86:[('SP86', '>u1'), ('e_err_m', '>u1'),], #Easting variance (+-X error)
            87:[('SP87', '>u1'), ('n_err_m', '>u1'),], #Northing variance
            152:[('SP152', '>u1'), ('unknown_152', '>u4'),], #Unknown
            153:[('SP153', '>u1'), ('unknown_153', '>u4'),], #Unknown
            154:[('SP154', '>u1'), ('unknown_154', '>u4'),], #Unknown
            155:[('SP155', '>u1'), ('unknown_155', '>u4'),], #Unknown
            156:[('SP156', '>u1'), ('unknown_156', '>u4'),], #Unknown
            157:[('SP157', '>u1'), ('unknown_157', '>u4'),], #Unknown
            158:[('SP158', '>u1'), ('unknown_158', '>u4'),], #Unknown
            159:[('SP159', '>u1'), ('unknown_159', '>u4'),], #Unknown
            160:[('SP160', '>u1'), ('ping_cnt', '>u4'),] #Number of ping values (in bytes)
            }

        file = open(sonFile, 'rb') # Open the file
        lastPos = 0 # Track last position in file
        head = self._fread_dat(file, 4,'B') # Get first 4 bytes of file

        # If first four bytes match known Humminbird ping header
        gotHeader = False
        if head[0] == 192 and head[1] == 222 and head[2] == 171 and head[3] == 33:
            headStruct.append(('head_start', '>u4'))
            while lastPos < headBytes - 1:
                lastPos = file.tell() # Get current position in file
                byte = self._fread_dat(file, 1, 'B')[0] # Decode the spacer byte
               
                if byte == 33:
                    file.seek(-6, 1)
                    nextByte = self._fread_dat(file, 1, 'B')[0]
                    if nextByte == 160:
                        header_len = lastPos + 1
                        lastPos = headBytes
                        headStruct.append(('head_end', '>u1'))
                
                else:
                    if byte in toCheck:
                        for v in toCheck[byte]:
                            headStruct.append(v)
                    else:
                        print('{} not in sonar header. Terminating.'.format(byte))
                        print('Offset: {}'.format(file.tell()))
                        sys.exit()
                    if byte < 100:
                        nextByte = 1
                    else:
                        nextByte = 4
                    lastPos = file.tell() + nextByte # Update with current position

                
                file.seek(lastPos)

        file.close() # Close the file

        self.son_struct = np.dtype(headStruct) # Store data in class attribute for later use
        return header_len

    def _parsePingHeader(self, in_file: str, out_file: str):
        '''
        '''

        # Get file length
        file_len = os.path.getsize(in_file)

        # Initialize counter
        i = 0
        chunk_i = 0
        chunk = 0


        header_dat_all = []

        frame_offset = []

        chunk_id = []

        file = open(in_file, 'rb')

        # Decode ping header
        while i < file_len:

            # Get header data at offset i
            header_dat, cpos = self._getPingHeader(file, i)

            # Add frame offset
            frame_offset.append(i)

            header_dat_all.append(header_dat)

            chunk_id.append(chunk)

            # update counter with current position
            i = cpos

            if chunk_i == self.nchunk:
                chunk_i = 0
                chunk += 1
            else:
                chunk_i += 1

        header_dat_all = pd.DataFrame.from_dict(header_dat_all)

        # Add in the frame offset
        header_dat_all['index'] = frame_offset

        # Add in the son_offset (headBytes for Humminbird)
        header_dat_all['son_offset'] = self.headBytes

        # Add chunk id
        header_dat_all['chunk_id'] = chunk_id

        # Do unit conversions
        header_dat_all = self._doUnitConversion(header_dat_all)
        

        # Drop spacer and unknown columns
        for col in header_dat_all.columns:
            if 'SP' in col:
                header_dat_all.drop(col, axis=1, inplace=True)

            if not self.exportUnknown and 'unknown' in col:
                header_dat_all.drop(col, axis=1, inplace=True)

        # Drop head_start
        header_dat_all.drop('head_start', axis=1, inplace=True)
        header_dat_all.drop('head_end', axis=1, inplace=True)

        # Update last chunk if too small (for rectification)
        lastChunk = header_dat_all[header_dat_all['chunk_id'] == chunk]
        if len(lastChunk) <= self.nchunk/2:
            header_dat_all.loc[header_dat_all['chunk_id'] == chunk, 'chunk_id'] = chunk-1


        # Save to csv
        header_dat_all.to_csv(out_file, index=False)

        return self.trans, self.humDat
    
    def _getPingHeader(self, file, i: int):

        # Get necessary attributes
        head_struct = self.son_struct
        length = self.frame_header_size # Account for start and end header

        # Move to offset
        file.seek(i)

        # Get the data
        buffer = file.read(length)

        # Read the data
        header = np.frombuffer(buffer, dtype=head_struct)

        out_dict = {}
        for name, typ in header.dtype.fields.items():
            out_dict[name] = header[name][0].item()

        # Next ping header is from current position + ping_cnt
        next_ping = file.tell() + header[0][-2]

        return out_dict, next_ping

    def _doUnitConversion(self, df: pd.DataFrame):

        '''
        '''

        # Calculate range
        df['max_range'] = df['ping_cnt'] * self.pixM

        # Easting and northing variances appear to be stored in the file
        ## They are reported in cm's so need to convert
        df['e_err_m'] = np.abs(df['e_err_m'])/100
        df['n_err_m'] = np.abs(df['n_err_m'])/100

        # Now calculate hdop from n/e variances
        df['hdop'] = np.round(np.sqrt(df['e_err_m']+df['n_err_m']), 2)

        # Get epsg code
        self._getEPSG(df['utm_e'].iloc[0], df['utm_n'].iloc[0])
        print('\n\n', df['utm_e'].iloc[0], df['utm_n'].iloc[0])

        # Convert eastings/northings to latitude/longitude (from Py3Hum - convert using International 1924 spheroid)
        lat = np.arctan(np.tan(np.arctan(np.exp(df['utm_n']/ 6378388.0)) * 2.0 - 1.570796326794897) * 1.0067642927) * 57.295779513082302
        lon = (df['utm_e'] * 57.295779513082302) / 6378388.0

        df['lon'] = lon
        df['lat'] = lat

        # Reproject latitude/longitude to UTM zone
        e, n = self.trans(lon, lat)
        df['e'] = e
        df['n'] = n

        # Instrument heading, speed, and depth need to be divided by 10 to be
        ## in appropriate units.
        df['instr_heading'] = df['instr_heading']/10
        df['speed_ms'] = df['speed_ms']/10
        df['inst_dep_m'] = df['inst_dep_m']/10

        # Get units into appropriate format
        df['f'] = df['f']/1000 # Hertz to Kilohertz
        df['time_s'] = df['time_s']/1000 #milliseconds to seconds
        df['tempC'] = self.tempC*10
        # # Can we figure out a way to base transducer length on where we think the recording came from?
        # df['t'] = 0.108
        # Use recording unix time to calculate each sonar records unix time
        try:
            # starttime = float(df['unix_time'])
            starttime = float(self.humDat['unix_time'])
            df['caltime'] = starttime + df['time_s']

        except :
            df['caltime'] = 0

        # Update caltime to timestamp
        sonTime = []
        sonDate = []
        needToFilt = False
        for t in df['caltime'].to_numpy():
            try:
                t = datetime.datetime.fromtimestamp(t)
                sonDate.append(datetime.datetime.date(t))
                sonTime.append(datetime.datetime.time(t))
            except:
                sonDate.append(-1)
                sonTime.append(-1)
                needToFilt = True
            
        df = df.drop('caltime', axis=1)
        df['date'] = sonDate
        df['time'] = sonTime

        if needToFilt:
            df = df[df['date'] != -1]
            df = df[df['time'] != -1]

            df = df.dropna()

        df = df[df['e'] != np.inf]
        df = df[df['record_num'] >= 0]

        lastIdx = df['index'].iloc[-1]
        df = df[df['index'] <= lastIdx]

        # Calculate along-track distance from 'time's and 'speed_ms'. Approximate distance estimate
        df = self._calcTrkDistTS(df)

        # Add transect number (for aoi processing)
        df['transect'] = 0

        # Store pixM
        df['pixM'] = self.pixM

        # Other corrections Dan did, not implemented yet...
        # if sonHead['beam']==3 or sonHead['beam']==2:
        #     dist = ((np.tan(25*0.0174532925))*sonHead['inst_dep_m']) +(tvg)
        #     bearing = 0.0174532925*sonHead['instr_heading'] - (pi/2)
        #     bearing = bearing % 360
        #     sonHead['heading'] = bearing
        # print("\n\n", sonHead, "\n\n")

        return df
    
    def _calcTrkDistTS(self,
                       df: pd.DataFrame):
        '''
        Calculate along track distance based on time ellapsed and gps speed.
        '''

        ts = df['time_s'].to_numpy()
        ss = df['speed_ms'].to_numpy()
        ds = np.zeros((len(ts)))

        # Offset arrays for faster calculation
        ts1 = ts[1:]
        ss1 = ss[1:]
        ts = ts[:-1]

        # Calculate instantaneous distance
        d = (ts1-ts)*ss1
        ds[1:] = d

        # Accumulate distance
        ds = np.cumsum(ds)

        df['trk_dist'] = ds
        return df


    #===========================================================================
    # END Humminbird to PINGMapper
    #===========================================================================

    #===========================================================================
    # Lowrance file to Humminbird (BETA)
    #===========================================================================
    def _makeOutFiles(self):

        # Make DAT file
        f = open(self.humFile, 'w')
        f.close()

        # Make son directory
        try:
            os.mkdir(self.sonFile)
        except:
            pass

        # Make son and idx files
        beams = np.arange(0, 5)
        for b in beams:
            son = os.path.join(self.sonFile, 'B00{}.SON'.format(b))
            f = open(son, 'w')
            f.close

            idx = son.replace('SON', 'IDX')
            f = open(idx, 'w')
            f.close()

        self.b000 = os.path.join(self.sonFile, 'B000.SON') #83kHz sonar
        self.b001 = os.path.join(self.sonFile, 'B001.SON') #200kHz sonar
        self.b002 = os.path.join(self.sonFile, 'B002.SON') #Port sonar
        self.b003 = os.path.join(self.sonFile, 'B003.SON') #Star sonar
        self.b004 = os.path.join(self.sonFile, 'B004.SON') #DownImage
    
    def _convertLowHeader(self, lowrance: low):

        '''
        Convert lowrance ping header (attributes)
        to humminbird. Keeping PING-Mapper naming
        conventions for humminbird attributes even
        though some fields indicate units other than
        what they are:

        ex: time_s indicates time is is seconds but will
        be stored in milliseconds, etc.

        For unknown attributes, simply adding default value
        from a sample sonar recording...

        Using the latest (2024) file format.
        frame_header_size == headBytes == 152
        '''

        # Set headBytes
        self.frame_header_size = 152

        # Create empty df
        df = pd.DataFrame()

        # Get lowrance ping attributes
        dfLow = lowrance.header_dat

        # Get record_num from index
        df['record_num'] = dfLow.index

        # Get time as ms
        ## Lowrance time in seconds
        df['time_s'] = ( dfLow['time'] * 1000 ).astype(int)

        # # UTM Easting
        # df['utm_e'] = dfLow['utm_e']

        # # UTM Northing
        # df['utm_n'] = dfLow['utm_n']

        # Humminbird uses a strange projection based on International 1924 ellipsoid
        ## In order to convert Lowrance to Humminbird coords, first convert to
        ## lat / lon then to Humminbird coords.
        df = self._convertLowCoordinates(df, dfLow)

        # Add gps1 (flag of some sort, unknown.)
        df['gps1'] = 1

        # Add instrument heading [radians to degrees]
        ## And multiply by 10
        df['instr_heading'] = ( np.rad2deg(dfLow['track_cog']) * 10 ).astype(int)

        # Add gps2 (flag of some sort, unknown.)
        df['gps2'] = 1

        # Speed [m/s to decimeters/second]
        df['speed_ms'] = ( (dfLow['gps_speed']) * 10 ).astype(int)

        # Unknown 134
        df['unknown_134'] = 0

        # Instrument depth
        df['inst_dep_m'] = ( ( dfLow['depth_ft'] ) * 10 ).astype(int)

        # unknown_136
        df['unknown_136'] = 1814532

        # unknown_137
        df['unknown_137'] = -1582119980

        # unknown_138
        df['unknown_138'] = -1582119980

        # unknown_139
        df['unknown_139'] = -1582119980

        # unknown_140
        df['unknown_140'] = -1582119980

        # unknown_141
        df['unknown_141'] = -1582119980

        # unknown_142
        df['unknown_142'] = -1582119980

        # unknown_143
        df['unknown_143'] = -1582119980

        # Beam
        df = self._convertLowBeam(df, dfLow)

        # Volt scale (?)
        df['volt_scale'] = 0#36

        # Frequency
        df = self._convertLowFrequency(df, dfLow)

        # unknown_83
        df['unknown_83'] = 18

        # unknown_84
        df['unknown_84'] = 1

        # unknown_149
        df['unknown_149'] = 26

        # Easting variance (+-X error);
        ## Unknown if this is actual value or if present in lowrance...setting to 0
        df['e_err_m'] = 0

        # Northing variance (+-X error);
        ## Unknown if this is actual value or if present in lowrance...setting to 0
        df['n_err_m'] = 0

        # unknown_152
        df['unknown_152'] = 4

        # unknown_155
        df['unknown_155'] = 3

        # unknown_156
        df['unknown_156'] = -1582119980

        # unknown_157
        df['unknown_157'] = -1582119980

        # unknown_158
        df['unknown_158'] = -1582119980

        # unknown_159
        df['unknown_159'] = -1582119980

        # ping_cnt
        df['ping_cnt'] = dfLow['packet_size']

        # Store frame offset
        df['frame_offset'] = dfLow['frame_offset']

        # Store son offset (from frame_offset)
        df['son_offset'] = lowrance.frame_header_size

        self.header_dat = df

        return

    def _convertLowCoordinates(self, df: pd.DataFrame, dfLow: pd.DataFrame):

        '''
        Humminbird uses International 1924 ellipsoid (epsg:4022????)
        Lowrance uses WGS 1984 ellipsoid (epsg:4326)
        '''

        ellip_1924 = 6378388.0

        # Convert eastings and northings into latitude and longitude based on wgs84 spheroid
        df['lat'] = ((2*np.arctan(np.exp(dfLow['utm_n']/6356752.3142)))-(np.pi/2))*(180/np.pi)
        df['lon'] = dfLow['utm_e']/6356752.3142*(180/np.pi)

        # # Get transformation epsg:7022
        # trans = pyproj.Proj('epsg:4022')

        # # Do transformation
        # df['utm_e'], df['utm_n'] = trans(df['lon'], df['lat'])

        # Conversion available in PING-Mapper and PyHum, but solved for northing / easting (sloppy...)
        df['utm_n'] = ellip_1924 * np.log( np.tan( ( np.arctan( np.tan( df['lat']/57.295779513082302 ) / 1.0067642927 ) + 1.570796326794897 ) / 2.0 ) )
        df['utm_e'] = ellip_1924 * (np.pi/180) * df['lon']


        return df
    
    def _convertLowBeam(self, dfHum: pd.DataFrame, dfLow: pd.DataFrame):

        '''
        Lowrance                Humminbird
        0 primary sounder       0 should be low frequency 83kHz
        1 secondary sounder     1 should be high frequency 200kHz
        2 downscan              4 downscan
        3 port ss               2 port ss
        4 star ss               3 star ss
        5 sidescan              NA Store as 5, convert in port star later
        '''

        # Store lowrance sidescan (5) as 5 and parse into port (2)
        ## and star (3) later..
        beam_xwalk = {0: 0, 1: 1, 2:4, 3:2, 4:3, 5:5}

        dfHum['beam'] = [beam_xwalk.get(i, "unknown") for i in dfLow['channel_type']]

        return dfHum
    
    def _convertLowFrequency(self, dfHum: pd.DataFrame, dfLow: pd.DataFrame):

        '''
        Crosswalk Lowrance frequency to Humminbird.
        Humminbird has slots for frequency, min-frequency, max-frequency

        {lowrance-frequency: [Humminbird Frequecy, min, max]}
        '''
        
        frequency_xwalk = {'200kHz': [200, 200, 200], '50kHz': [50, 50, 50],
                           '83kHz': [83, 83, 83], '455kHz': [455, 455, 455],
                           '800kHz': [800, 800, 800], '38kHz': [38, 38, 38],
                           '28kHz': [28, 28, 28], '130kHz_210kHz': [170, 130, 210],
                           '90kHz_150kHz': [120, 90, 150], '40kHz_60kHz': [50, 40, 60],
                           '25kHz_45kHz': [35, 25, 45]}
        
        frequency_min = {200: 200, 50: 50, 83: 83, 455: 455, 800: 800, 38: 38,
                         28: 28, 170: 130, 120:90, 50: 40, 35: 25}
        
        dfHum['f'] = [frequency_xwalk[i][0] for i in dfLow['frequency']]
        dfHum['f_min'] = [frequency_xwalk[i][1] for i in dfLow['frequency']]
        dfHum['f_max'] = [frequency_xwalk[i][2] for i in dfLow['frequency']]

        return dfHum

    def _removeUnknownBeams(self):

        df = self.header_dat

        # Drop unknown
        df = df[df['beam'] != 'unknown']

        self.header_dat = df
        return

    def _splitLowSS(self):
        '''
        If beam 5 present in lowrance, then port and starboard ss are merged.
        Must be split to export into their own files.
        '''

        # Get dataframe
        dfAll = self.header_dat

        # Get beam 5
        df = dfAll[dfAll['beam'] == 5]

        # Make copies, one for port, other for star
        port = df.copy()
        star = df.copy()

        # Re-label beam numbers
        port['beam'] = 2
        star['beam'] = 3

        # Divide ping_cnt in half
        port['ping_cnt'] = (port['ping_cnt'] / 2).astype(int)
        star['ping_cnt'] = (star['ping_cnt'] / 2).astype(int)

        # Assume left half are port returns and right are starboard
        # Add additional offset to star the account for this
        star['son_offset'] += star['ping_cnt']

        # Remove beam 5 from dfAll
        dfAll = dfAll[dfAll['beam'] != 5]

        # Concatenate df's
        dfAll = pd.concat([dfAll, port, star], ignore_index=True)

        dfAll.sort_values(by=['time_s', 'beam'], inplace=True)

        self.header_dat = dfAll

        return

    def _recalcRecordNum(self):

        df = self.header_dat

        # Reset index and recalculate record num
        ## Record num is unique for each ping across all sonar beams
        df = df.reset_index(drop=True)
        df['record_num'] = df.index

        self.header_dat = df
        return

    def _convertLowDAT(self, lowrance: low):

        '''
        Humminbird recordings need a DAT pointer file
        '''

        # Dictionary to store data
        dat = dict()

        # Get ping attributes
        dfHum = self.header_dat
        dfLow = lowrance.header_dat

        # Unknown spacer
        dat['SP1'] = 195

        # Water code; unsure if present in Lowrance, setting to 1 (freshwater) for now
        dat['water_code'] = 1

        # Unknown spacer
        dat['SP2'] = 125

        # unknown_1
        dat['unknown_1'] = 1

        # Sonar name (??)
        dat['sonar_name'] = 1029

        # unknown_2
        dat['unknown_2'] = 11

        # unknown_3
        dat['unknown_3'] = 0

        # unknown_4
        dat['unknown_4'] = 0

        # unix_time
        dat['unix_time'] = dfLow['creation_date_time'][0].item()

        # utm_e
        dat['utm_e'] = dfHum['utm_e'][0].item()

        # utm_n
        dat['utm_n'] = dfHum['utm_n'][0].item()

        # Filename
        dat['filename'] = os.path.basename(self.b002)

        # Number of records
        dat['numrecords'] = len(dfHum)

        # Recording length in milliseconds
        dat['recordlens_ms'] = dfHum.iloc[-1]['time_s'].item()

        # linesize: size of ping frame = frame_head_size + ping_cnt
        dat['linesize'] = self.frame_header_size + dfHum['ping_cnt'][0].item()

        # unknown_5
        dat['unknown_5'] = 5

        # unknown_6
        dat['unknown_6'] = 30

        # unknown_7
        dat['unknown_7'] = dat['sonar_name']

        # unknown_8
        dat['unknown_8'] = dat['sonar_name']

        # unknown_9
        dat['unknown_9'] = 0

        # unknown_10
        dat['unknown_10'] = -1582119980

        # unknown_11
        dat['unknown_11'] = -1582119980

        # unknown_12
        dat['unknown_12'] = -1582119980

        # unknown_13
        dat['unknown_13'] = -1582119980

        # unknown_14
        dat['unknown_14'] = -1582119980

        self.dat = dat

        return
    
    def _writeDAT(self):

        '''
        Write dat contents to DAT file
        '''

        # Get DAT struct
        if self.frame_header_size == 152:
            # humDic = {
            #             'endianness':'<i', #<=little endian; I=unsigned Int
            #             'SP1':[0, 0, 1, -1], #Unknown (spacer)
            #             'water_code':[1, 0, 1, -1], #Need to check if consistent with other models (1=fresh?)
            #             'SP2':[2, 0, 1, -1], #Unknown (spacer)
            #             'unknown_1':[3, 0, 1, -1], #Unknown (gps flag?)
            #             'sonar_name':[4, 0, 4, -1], #Sonar name
            #             'unknown_2':[8, 0, 4, -1], #Unknown
            #             'unknown_3':[12, 0, 4, -1], #Unknown
            #             'unknown_4':[16, 0, 4, -1], #Unknown
            #             'unix_time':[20, 0, 4, -1], #Unix Time
            #             'utm_e':[24, 0, 4, -1], #UTM X
            #             'utm_n':[28, 0, 4, -1], #UTM Y
            #             'filename':[32, 0, 12, -1], #Recording name
            #             'numrecords':[44, 0, 4, -1], #Number of records
            #             'recordlens_ms':[48, 0, 4, -1], #Recording length milliseconds
            #             'linesize':[52, 0, 4, -1], #Line Size (?)
            #             'unknown_5':[56, 0, 4, -1], #Unknown
            #             'unknown_6':[60, 0, 4, -1], #Unknown
            #             'unknown_7':[64, 0, 4, -1], #Unknown
            #             'unknown_8':[68, 0, 4, -1], #Unknown
            #             'unknown_9':[72, 0, 4, -1], #Unknown
            #             'unknown_10':[76, 0, 4, -1], #Unknown
            #             'unknown_11':[80, 0, 4, -1], #Unknown
            #             'unknown_12':[84, 0, 4, -1], #Unknown
            #             'unknown_13':[88, 0, 4, -1], #Unknown
            #             'unknown_14':[92, 0, 4, -1]
            #             }
            
            dat_dtype = ([
                ('SP1', '<u1'),
                ('water_code', '<u1'),
                ('SP2', '<u1'),
                ('unknown_1', '<u1'),
                ('sonar_name', '<u4'),
                ('unknown_2', '<u4'),
                ('unknown_3', '<u4'),
                ('unknown_4', '<u4'), 
                ('unix_time', '<u4'),
                ('utm_e', '<i4'),
                ('utm_n', '<i4'),
                ('filename', 12),
                ('numrecords', '<u4'),
                ('recordlens_ms', '<u4'), 
                ('linesize', '<u4'),
                ('unknown_5', '<u4'),
                ('unknown_6', '<u4'),
                ('unknown_7', '<u4'),
                ('unknown_8', '<u4'),
                ('unknown_9', '<u4'),
                ('unknown_10', '<i4'),
                ('unknown_11', '<i4'),
                ('unknown_12', '<i4'),
                ('unknown_13', '<i4'),
                ('unknown_14', '<i4'),
            ])
            
        for i in dat_dtype:
            name = i[0]
            dtype = i[1]

            if name != 'filename':
                val = np.array(self.dat[name], dtype=dtype)

                with open(self.humFile, 'ab') as file:
                    file.write(val)
            else:
                val = self.dat[name]
                topad = dtype - len(val)
                s = 0
                while s < topad:
                    val += ' '
                    s += 1
                
                with open(self.humFile, 'a') as f:
                    f.write(val)

        return
    
    def _writeSonfromLow(self, beam: int, header_size: int, lowrance_path: str, flip_port: bool = False):

        '''
        Each ping attribute in the header of a Humminbird SON file
        has a tag preceding the attribute. Conversely, Lowrance has
        one attribute followed by another. Therefore, the tag must 
        be inserted while writing the data from Lowrance to Humminbird.

        son_dtype: ([('attribute_name', tag value, attribute_dtype)])

        *** Big Endian > ***
        '''

        son_dtype = ([
            ('head_start', 3235818273, '>u4'),
            ('record_num', 128, '>u4'),
            ('time_s', 129, '>u4'),
            ('utm_e', 130, '>i4'),
            ('utm_n', 131, '>i4'),
            ('gps1', 132, '>u2'),
            ('instr_heading', 132.2, '>u2'),
            ('gps2', 133, '>u2'),
            ('speed_ms', 133.2, '>u2'),
            ('unknown_134', 134, '>u4'),
            ('inst_dep_m', 135, '>u4'),
            ('unknown_136', 136, '>i4'),
            ('unknown_137', 137, '>i4'),
            ('unknown_138', 138, '>i4'),
            ('unknown_139', 139, '>i4'),
            ('unknown_140', 140, '>i4'),
            ('unknown_141', 141, '>i4'),
            ('unknown_142', 142, '>i4'),
            ('unknown_143', 143, '>i4'),
            ('beam', 80, '>u1'),
            ('volt_scale', 81, '>u1'),
            ('f', 146, '>u4'),
            ('unknown_83', 83, '>u1'),
            ('unknown_84', 84, '>u1'),
            ('unknown_149', 149, '>u4'),
            ('e_err_m', 86, '>u1'),
            ('n_err_m', 87, '>u1'),
            ('unknown_152', 152, '>u4'),
            ('f_min', 153, '>u4'),
            ('f_max', 154, '>u4'),
            ('unknown_155', 155,'>u4'),
            ('unknown_156', 156,'>i4'),
            ('unknown_157', 157,'>i4'),
            ('unknown_158', 158,'>i4'),
            ('unknown_159', 159,'>i4'),
            ('ping_cnt', 160, '>u4'),
            ('head_end', 33, '>u1')
            ])
        
        if beam == 0:
            file_name = self.b000
        elif beam == 1:
            file_name = self.b001
        elif beam == 2:
            file_name = self.b002
        elif beam == 3:
            file_name = self.b003
        elif beam == 4:
            file_name = self.b004
        else:
            sys.exit('{} not a valid beam.')

        # Get the header_dat
        df = self.header_dat

        # Filter df based off beam
        df = df[df['beam'] == beam]

        # Track ping offset
        offset = 0

        # Get IDX file path
        idx_file = file_name.replace('SON', 'IDX')

        # Iterate df rows
        for i, row in df.iterrows():

            # # For IDX
            # idx = []

            # Convert row to a dictionary
            row = row.to_dict()

            with open(file_name, 'ab') as file:
                # Iterate son_dtype
                for i in son_dtype:

                    buffer = []

                    name = i[0]
                    tag_val = i[1]
                    dtype = i[2]

                    if name == 'head_start' or name == 'head_end':
                        spacer = np.array(tag_val, dtype=dtype)
                        val = -9999
                        a = 0
                    elif isinstance(tag_val, float):
                        spacer = -9999
                        val = np.array(row[name], dtype=dtype)
                        a = 1
                    else:
                        spacer = np.array(tag_val, '>u1')
                        val = np.array(row[name], dtype=dtype)
                        a = 2

                    if spacer != -9999:
                        file.write(spacer)
                        buffer.append(spacer)

                    if val != -9999:
                        file.write(val)
                        buffer.append(val)

                    del spacer, val

                # Get the ping returns
                ping_returns = self._getLowPingReturns(lowrance_path, row['frame_offset'], row['son_offset'], row['ping_cnt'])

                if flip_port:
                    ping_returns = ping_returns[::-1]

                # Write returns to file
                file.write(ping_returns)

            # Write time and offset to IDX
            with open(idx_file, 'ab') as file:

                time = np.array(row['time_s'], '>u4')
                offset = np.array(offset, '>u4')

                file.write(time)
                file.write(offset)

            # Update offset
            # Offset just size of IDX?????
            offset = os.path.getsize(file_name)

    def _getLowPingReturns(self, file: str, offset: int, son_offset: int, length: int):

        # Open file
        f = open(file, 'rb')

        # Move to position
        f.seek(offset + son_offset)

        # Get the data
        buffer = f.read(length)

        f.close()

        return buffer

    #===========================================================================
    # END Lowrance file to Humminbird
    #===========================================================================

    # ======================================================================
    def __str__(self):
        '''
        Generic print function to print contents of sonObj.
        '''
        output = "Humminbird Class Contents"
        output += '\n\t'
        output += self.__repr__()
        temp = vars(self)
        for item in temp:
            output += '\n\t'
            output += "{} : {}".format(item, temp[item])
        return output




        




        