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
import numpy as np
import pandas as pd
from datetime import datetime, timezone, timedelta
import pyproj

# Add 'pingmapper' to the path, may not need after pypi package...
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(PACKAGE_DIR)

from pingverter.verter_utils import filterGPS

# # RSD structur
# rsdStruct = np.dtype([
#     ("test", "<u4"),
# ])

garCols2PM = {
    'bottom_depth': 'inst_dep_m',
    'drawn_bottom_depth': 'keel_depth_m',
    'sample_cnt': 'ping_cnt', 
    'first_sample_depth': 'min_range',
    'last_sample_depth': 'max_range',
    'water_temp': 'tempC',
    'recording_time_ms': 'time_s',
}

class gar(object):

    #===========================================================================
    def __init__(self, inFile: str, nchunk: int=0, exportUnknown: bool=False):
        
        '''
        '''

        self.humFile = None
        self.sonFile = inFile
        self.nchunk = nchunk
        self.exportUnknown = exportUnknown

        self.magicNum = 3085556358



        self.extension = os.path.basename(inFile).split('.')[-1]

        # self.son_struct = rsdStruct

        self.garCols2PM = garCols2PM

        self.humDat = {} # Store general sonar recording metadata

        self.son8bit = False

        # Set Sonar beams
        # Neils beams: 1,3,4,5
        # UD beams: 0,2,4,6
        # Unsure and need to test
        # self.beam_set = {
        #     2: ['ds_lowfreq', 0],
        #     1: ['ds_hifreq', 1],
        #     5: ['ss_port', 2],
        #     6: ['ss_port', 2],
        #     4: ['ss_star', 3],
        #     0: ['ds_vhighfreq', 4],
        #     3: ['ds_vhighfreq', 5],
        # }

        return
    
    # ======================================================================
    def _getFileLen(self):
        self.file_len = os.path.getsize(self.sonFile)

        return
    
    # ======================================================================
    def _fread_dat(self,
            infile,
            num,
            typ):
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

        buffer = infile.read(num)
        data = np.frombuffer(buffer, dtype=typ)

        return data
    
    
    ### File Header ###
    # ======================================================================
    def _parseFileHeader(self):
        '''
        '''

        self.headBytes = 20480 # Hopefully a fixed value for all RSD files
        chanInfoLen = 1069 # It is not clear if there is helpful info in channel information...

        # Get the file header structure
        headStruct, firstHeadBytes = self._getFileHeaderStruct()

        # print('\n\n\nheadStruct:')
        # for v in headStruct:
        #     print(v)
        
        # Read the header
        file = open(self.sonFile, 'rb') # Open the file
        file.seek(0)

        # Get the data
        buffer = file.read(firstHeadBytes)

        # Read the data
        header = np.frombuffer(buffer=buffer, dtype=headStruct)

        out_dict = {}
        for name, typ in header.dtype.fields.items():
            out_dict[name] = header[name][0].item()

        # for k,v in out_dict.items():
        #     print(k, v)

        self.file_header = out_dict


        return
    
    # ======================================================================
    def _getFileHeaderStruct(self):
        '''

        ffh: field - file header
        ffi: field - file information
        fci: field - channel information
        fcnt: field count
        '''


        # headBytes = 20480
        # firstHeadBytes = 35
        headStruct = [] 

        toCheck = {
            6:[('header_fcnt', '<u1')], #06: number of fields in header structure
            4:[('ffh_0', '<u1'), ('magic_number', '<u4')], #04: field 0 "magic_number", length 4
            10:[('ffh_1', '<u1'), ('format_version', '<u2')], #0a: field 1 "format_version", length 2
            20:[('ffh_2', '<u1'), ('channel_count', '<u4')], #14: field 2 "channel_count", length 4
            25:[('ffh_3', '<u1'), ('max_channel_count', '<u1')], #19: field 3 "max_channel_count", length 1
            47:[('ffh_4', '<u1'), ('ffh_4_actlen', '<u1'), ('ffi_fcnt', '<u1')], #2f: field 4 "file information", length 7; #actual length; #number of "file information" field 
            55: [('ffh_5', '<u1'), ('ffh_5_actlen', '<u2'), ('chan_cnt', '<u1'), ('fci_acnt', '<u1')], #37: "channel information", length 7; 
        }

        fileInfoToCheck = {
            2: [('ffi_0', '<u1'), ('unit_software_version', '<u2')], #02: "file information" field 0
            12: [('ffi_1', '<u1'), ('unit_id_type', '<u4')], #0c: "file information" field 1
            18: [('ffi_2', '<u1'), ('unit_product_number', '<u2')], #12: "file information" field 2    
            28: [('ffi_3', '<u1'), ('date_time_of_recording', '<u4')], #1c: "file information" field 3
        }

        # chanInfoToCheck = {
        #     3: [('fci_0_data_info', '<u1'), ('fci_0_acnt', '<u1'), ('fci_0_actlen', '<u1'), ('fci_channel_id', '<u1')], #03: Channel data info
        #     15: [('fci_1', '<u1'), ('fci_1_actlen', '<u1'), ('fci_first_chunk_offset', '<u8')], #0f: first chunk offset
        #     23: [('fci_2', '<u1'), ('fci_2_actlen', '<u2'), ('fci_2_acnt', '<u1'), ('fci_2_actlen2', '<u2')], #17 prop_chan_info

        
        # }

        # chanDataInfoToCheck = {

        # }

        file = open(self.sonFile, 'rb') # Open the file
        lastPos = 0 # Track last position in file

        foundChanInfo = False

        # while lastPos < firstHeadBytes - 1:
        while not foundChanInfo:
            # lastPos = file.tell()
            # print('lastPos:', lastPos)
            byte = self._fread_dat(file, 1, 'B')[0] # Decode the spacer byte

            if byte == 47:
                # File Information
                structDict = fileInfoToCheck

                for v in toCheck[byte]:
                    headStruct.append(v)

                length = self._fread_dat(file, 1, 'B')[0] # Decode the spacer byte
                field_cnt = self._fread_dat(file, 1, 'B')[0] # Decode the spacer byte

                fidx = 0
                while fidx < field_cnt:
                    byte = self._fread_dat(file, 1, 'B')[0] # Decode the spacer byte
                    if byte in structDict:
                        elen = 0
                        for v in structDict[byte]:

                            headStruct.append(v)
                            
                            # Get length of element
                            elen += (np.dtype(v[-1]).itemsize)

                        # Move forward elen amount
                        cpos = file.tell()
                        npos = cpos + elen - 1

                        file.seek(npos)
                        fidx += 1
                    else:
                        print('{} not in sonar header. Terminating.'.format(byte))
                        print('Offset: {}'.format(file.tell()))
                        sys.exit()

                # lastPos = headBytes
            
            # elif byte == 55:
            #     # File Information

            #     for v in toCheck[byte]:
            #         headStruct.append(v)

            #     length = self._fread_dat(file, 2, '<u2')[0] # Decode the spacer byte
            #     chan_cnt = self._fread_dat(file, 1, 'B')[0] # Decode the spacer byte

            #     chanidx = 0

            #     # Iterate each channel
            #     while chanidx < chan_cnt:

            #         field_cnt_0 = self._fread_dat(file, 1, 'B')[0] # Decode the spacer byte
            #         fidx_0 = 0

            #         # Iterate each field
            #         while fidx_0 < field_cnt_0:
            #             byte = self._fread_dat(file, 1, 'B')[0] # Decode the spacer byte

            #             if byte in chanInfoToCheck:
            #                 for v in chanInfoToCheck[byte]:
            #                     # Add chanidx to field name
            #                     field_name = '{}_{}'.format(v[0], chanidx)
            #                     v = (field_name, v[1])

            #                     headStruct.append(v)

            #             fidx_0 += 1

            #         chanidx += 1
                
            #     lastPos = headBytes

            elif byte == 55:
                foundChanInfo = True                            

            else:
                if byte in toCheck:
                    elen = 0
                    for v in toCheck[byte]:
                        headStruct.append(v)

                        # Get length of element
                        elen += (np.dtype(v[-1]).itemsize)
                    
                    # Move forward elen amount
                    cpos = file.tell()
                    npos = cpos + elen - 1

                    file.seek(npos)

            lastPos = file.tell()
            # print('lastPos:', lastPos)        

        return headStruct, lastPos-1
    
    
    ### Ping Header ###
    # ======================================================================
    def _parsePingHeader(self,):
        '''
        '''

        # Get the header struct
        self.son_struct, self.son_header_struct, self.record_body_header_len = self._getPingHeaderStruct()

        # Get the file length
        file_len = self.file_len

        # Initialize offset after file header
        i = self.headBytes

        # Open the file
        file = open(self.sonFile, 'rb')

        # Store contents in list
        header_dat_all = []

        # Decode ping header
        while i < file_len:

            # Get header data at offset i
            header_dat, cpos = self._getPingHeader(file, i)

            if header_dat:
                header_dat_all.append(header_dat)

            i = cpos

        # Convert to dataframe
        df = pd.DataFrame.from_dict(header_dat_all)

        # Convert fields
        df = self._doUnitConversion(df)

        # Do column name conversions to PINGMapper units
        df.rename(columns=self.garCols2PM, inplace=True)

        # Calculate speed & track distance (based on coords and time)
        df = self._calcSpeedTrkDist(df)

        # Drop negative son_offset
        df = df[df['son_offset'] > 0]


        # Test file to see outputs
        out_test = os.path.join(self.metaDir, 'All-Garmin-Sonar-MetaData.csv')
        df.to_csv(out_test, index=False)

        # Store in class
        self.header_dat = df

        return
    
    # ======================================================================
    def _getPingHeaderStruct(self, ):
        '''
        fpf: field - ping field
        fps: field - ping state
        '''

        headBytes = self.headBytes # Header length

        headStruct = [] 

        # pingHeaderToCheck = {
        #     6:[('header_fcnt', '<u1')], #06: number of fields in header structure
        #     4:[('fpf_0', '<u1'), ('magic_number', '<u4')], #04: field 0 "magic_number", length 4
        #     15:[('fpf_1', '<u1'), ('fpf_1_len', '<u1'), ('fpf_1_fcnt', '<u1'),
        #         ('fps_0', '<u1'), ('state', '<u1'),
        #         ('fps_1', '<u1'), ('data_info', '<u1'), ('data_info_cnt', '<u1'), ('data_info_len', '<u1'), ('channel_id', '<u1')], #0f: state data structure
        #     20:[('SP14', '<u1'), ('sequence_cnt', '<u4')], #14: sequence_count
        #     28:[('SP1c', '<u1'), ('data_crc', '<u4')], #1c: data crc
        #     34:[('SP22', '<u1'), ('data_size', '<u4')], #22: data size
        #     44:[('SP2c', '<u1'), ('recording_time_ms', '<u4')], #2c recording time offset
        # }

        # Record header (len==37)
        self.pingHeaderLen = pingHeaderLen = 37
        self.pingHeaderLenFirst = pingHeaderLenFirst = 49
        pingHeader = [
            ('header_fcnt', '<u1'),
            ('fpf_0', '<u1'), 
            ('magic_number', '<u4'),
            ('fpf_1', '<u1'), 
            ('fpf_1_len', '<u1'), 
            ('fpf_1_fcnt', '<u1'),
            ('fps_0', '<u1'), 
            ('state', '<u1'),
            ('fps_1', '<u1'), 
            # ('data_info', '<u1'), 
            ('data_info_cnt', '<u1'), 
            ('data_info_len', '<u1'), 
            ('channel_id', '<u1'),
            ('SP14', '<u1'), 
            ('sequence_cnt', '<u4'),
            ('SP1c', '<u1'), 
            ('data_crc', '<u4'),
            ('SP22', '<u1'),
            ('data_size', '<u2'),
            ('SP2c', '<u1'), 
            ('recording_time_ms', '<u4'),
            ('record_crc', '<u4')
        ]

        # pingBodyHeaderToCheck = {
        #     -1:('record_body_fcnt', '<u1'),
        #     1:[('SP1', '<u1'), ('channel_id_1', '<u1')], #01 channel_id
        #     11:[('SP0b', '<u1'), ('bottom_depth_unknown', '<u1'), ('bottom_depth', '<u2')], #0b bottom depth
        #     13:[('SP0d', '<u1'), ('unknown_sp0d', '<u4'), ('unknown_sp0d_1', '<u1')],
        #     18:[('SP12', '<u1'), ('unknown_sp12', '<u2')],
        #     19:[('SP13', '<u1'), ('drawn_bottom_depth_unknown', '<u1'), ('drawn_bottom_depth', '<u2')], #13 drawn bottom depth
        #     21:[('SP15', '<u1'), ('unknown_sp15', '<u4'), ('unknown_sp15_1', '<u1')],
        #     25:[('SP19', '<u1'), ('first_sample_depth', '<u1')], #19 first sample depth
        #     35:[('SP23', '<u1'), ('last_sample_depth_unknown', '<u1'), ('last_sample_depth', '<u2')], #23 last sample depth
        #     41:[('SP29', '<u1'), ('gain', '<u1')], #29 gain
        #     49:[('SP31', '<u1'), ('sample_status', '<u1')], #31 sample status
        #     60:[('SP3c', '<u1'), ('sample_cnt', '<u4')], #3c sample count
        #     65:[('SP41', '<u1'), ('shade_avail', '<u1')], #41 shade available
        #     76:[('SP4c', '<u1'), ('scposn_lat', '<u4')], #4c latitude
        #     84:[('SP54', '<u1'), ('scposn_lon', '<u4')], #54 longitude
        #     92:[('SP5c', '<u1'), ('water_temp', '<f4')], #5c temperature
        #     97:[('SP61', '<u1'), ('beam', '<u1')], #61 beam
        # }

        # magic number 86 DA E9 B7 == 3085556358
        # Ping headers always start at 20480

        # First and last state is 1. First time is the pingHeaderToCheck, forllowed
        ## by 16 CRC values, followed by first ping.

        # Beam 2 & 3 have extra unknown beam info (length 63)
        # Beam 1 & 4 sonar data starts immediately after 'beam'

        start_pos = lastPos = headBytes

        # Open file and move to offset
        file = open(self.sonFile, 'rb') # Open the file
        file.seek(start_pos)

        foundChanInfo = False

        headStruct = []

        # Start reading
        while not foundChanInfo:

            # Get the ping header (should be fixed structure)
            buffer = file.read(pingHeaderLen)

            # Read the data
            header = np.frombuffer(buffer, dtype=np.dtype(pingHeader))

            # Parse the data
            out_dict = {}
            for name, typ in header.dtype.fields.items():
                out_dict[name] = header[name][0].item()

            # Check if there is a record body
            if out_dict['state'] == 1: # no record body
                lastPos += pingHeaderLenFirst
                file.seek(lastPos)

            else:

                # # # Add pingheader
                # # for i in pingHeader:
                # #     headStruct.append(i)

                # # Get field count
                # field_cnt = self._fread_dat(file, 1, 'B')[0]
                # headStruct.append(pingBodyHeaderToCheck[-1])

                # if field_cnt > 13: # Only 13 known fields. Some beams have up to 15
                #     field_cnt = 13

                # fidx = 0
                # record_body_header_len = 1

                # while fidx < field_cnt:

                #     byte = self._fread_dat(file, 1, 'B')[0]

                #     if byte in pingBodyHeaderToCheck:
                #         elen = 0
                #         for v in pingBodyHeaderToCheck[byte]:
                #             headStruct.append(v)

                #             # Get length of element
                #             elen += (np.dtype(v[-1]).itemsize)

                #         # Move forward elen amount
                #         cpos = file.tell()
                #         npos = cpos + elen - 1

                #         record_body_header_len += elen

                #         file.seek(npos)

                #         fidx += 1

                #     else:
                #         print('{} not in sonar body. Terminating.'.format(byte))
                #         print('Offset: {}'.format(file.tell()))
                #         sys.exit()

                    foundChanInfo = True

        # self.son_header_struct = pingHeader
        # self.son_struct = headStruct

        # return headStruct, pingHeader, record_body_header_len
        return headStruct, pingHeader, 0

    # ======================================================================
    def _getPingHeader(self, file, i: int):

        # print('\n\n\n', i)

        # Get necessary attributes
        son_header_struct = self.son_header_struct
        pingHeaderLen = self.pingHeaderLen

        # head_struct = self.son_struct
        # record_body_header_len = self.record_body_header_len

        # Move to offset
        file.seek(i)

        # Get the ping header
        buffer = file.read(pingHeaderLen)

        # Read the data
        header = np.frombuffer(buffer, dtype=np.dtype(son_header_struct))

        out_dict = {}
        for name, typ in header.dtype.fields.items():
            out_dict[name] = header[name][0].item()

        # Check if there is a record body
        if out_dict['state'] != 2: # no record bod
            return False, i + self.pingHeaderLenFirst
        
        # # Get record body
        # # Get the ping header
        # buffer = file.read(record_body_header_len)

        # # Read the data
        # header = np.frombuffer(buffer, dtype=np.dtype(head_struct))

        # for name, typ in header.dtype.fields.items():
        #     out_dict[name] = header[name][0].item()

        # Variable structure so above doesn't work
        # Must determine structure ping by ping

        pingBodyHeaderToCheck = {
            -1:('record_body_fcnt', '<u1'),
            1:[('SP1_bh', '<u1'), ('channel_id_1', '<u1')], #01 channel_id
            10:[('SP0a', '<u1'), ('bottom_depth', 'V2')],
            # 11:[('SP0b', '<u1'), ('bottom_depth_unknown', '<u1'), ('bottom_depth', '<u2')], #0b bottom depth
            11:[('SP0b', '<u1'), ('bottom_depth', 'V3')], #0b bottom depth
            13:[('SP0d', '<u1'), ('unknown_sp0d', 'V5')],
            18:[('SP12', '<u1'), ('drawn_bottom_depth', 'V2')],
            19:[('SP13', '<u1'), ('drawn_bottom_depth', 'V3')], #13 drawn bottom depth
            21:[('SP15', '<u1'), ('unknown_sp15', 'V5')],
            25:[('SP19', '<u1'), ('first_sample_depth', '<u1')], #19 first sample depth
            35:[('SP23', '<u1'), ('last_sample_depth', 'V3')], #23 last sample depth
            41:[('SP29', '<u1'), ('gain', '<u1')], #29 gain
            49:[('SP31', '<u1'), ('sample_status', '<u1')], #31 sample status
            60:[('SP3c', '<u1'), ('sample_cnt', '<u4')], #3c sample count
            65:[('SP41', '<u1'), ('shade_avail', '<u1')], #41 shade available
            76:[('SP4c', '<u1'), ('scposn_lat', '<u4')], #4c latitude
            84:[('SP54', '<u1'), ('scposn_lon', '<u4')], #54 longitude
            92:[('SP5c', '<u1'), ('water_temp', '<f4')], #5c temperature
            97:[('SP61', '<u1'), ('beam', '<u1')], #61 beam
        }

        beamInfoToCheck = {
            # 111:[('SP6f', '<u1'), ('bi_len', '<u1')],
            1:[('SP1_bi', '<u1'), ('port_star_beam_angle', '<u1')],
            9:[('SP9', '<u1'), ('fore_aft_beam_angle', '<u1')],
            17:[('SP11', '<u1'), ('port_star_elem_angle', '<u1')],
            25:[('SP19_bi', '<u1'), ('fore_aft_elem_angle', '<u1')],
            47:[('SP2f', '<u1'), ('su2_len', '<u1'), ('su2_fcnt', '<u1'),
                ('su2_f0', '<u1'), ('port_star_id', '<f4'),
                ('su2_f1', '<u1'), ('su2_f1_unkown', '<f4'),
                ],
            55:[('SP37', '<u1'), ('su3_len', '<u1'), ('su3_fcnt', '<u1'),
                ('su3_f0', '<u1'), ('su3_f0_unknown', '<u1'),
                ('su3_f1', '<u1'), ('su3_f1_unkown', '<f4'),
                ('su3_f2', '<u1'), ('su3_f2_unkown', '<f4'),
                ('su3_f3', '<u1'), ('su3_f3_unkown', '<f4'),
                ('su3_f4', '<u1'), ('su3_f4_unkown', '<f4'),
                ('su3_f5', '<u1'), ('su3_f5_unkown', '<f4'),
                ('su3_f6', '<u1'), ('su3_f6_unkown', '<f4'),
                ],
            115:[('SP73', '<u1'), ('interrogation_id', '<u2'), ('son_byte_len', '<u1')]
        
        }

        beam_info = False

        # Get field count
        rb_field_cnt = field_cnt = self._fread_dat(file, 1, 'B')[0]
        out_dict['record_body_fcnt'] = field_cnt

        if rb_field_cnt > 13: # Only 13 known fields. Some beams have up to 15
            field_cnt = 13
            beam_info = True

        fidx = 0
        record_body_header_len = 0

        while fidx < field_cnt:

            byte = self._fread_dat(file, 1, 'B')[0]

            if byte in pingBodyHeaderToCheck:
                # Add byte
                out_dict[pingBodyHeaderToCheck[byte][0][0]] = byte

                son_struct = pingBodyHeaderToCheck[byte][1:]

                elen = 0
                for v in son_struct:
                    elen += np.dtype(v[-1]).itemsize

                buffer = file.read(elen)

                # Read the data
                header = np.frombuffer(buffer, dtype=np.dtype(son_struct))

                for name, typ in header.dtype.fields.items(): # type: ignore
                    out_dict[name] = header[name][0].item()
                
                fidx += 1
                record_body_header_len += elen


        

        if beam_info:

            fid_beam_info = self._fread_dat(file, 1, 'B')[0]
            bi_len = self._fread_dat(file, 1, 'B')[0]
            bi_fcnt = self._fread_dat(file, 1, 'B')[0]

            fidx = 0

            while fidx < bi_fcnt:

                byte = self._fread_dat(file, 1, 'B')[0]

                if byte in beamInfoToCheck:
                    # Add byte
                    out_dict[beamInfoToCheck[byte][0][0]] = byte

                    son_struct = beamInfoToCheck[byte][1:]

                    elen = 0
                    for v in son_struct:
                        elen += np.dtype(v[-1]).itemsize

                    buffer = file.read(elen)

                    # Read the data
                    header = np.frombuffer(buffer, dtype=np.dtype(son_struct))

                    for name, typ in header.dtype.fields.items(): # type: ignore
                        out_dict[name] = header[name][0].item()
                    
                    fidx += 1
                    record_body_header_len += elen

            


        # Next ping header is from current position + ping_cnt
        # next_ping = file.tell() + out_dict['packet_size']
        next_ping = i + pingHeaderLen + out_dict['data_size'] + 12 #12 for magic number & crc

        out_dict['index'] = i

        out_dict['son_offset'] = (out_dict['data_size']) - (out_dict['sample_cnt']*2) + self.pingHeaderLen

        # out_dict['son_offset'] = record_body_header_len+1
 
        return out_dict, next_ping
    
    
    ### Ping Header Conversions ###
    # ======================================================================
    def _calcSpeedTrkDist(self, df: pd.DataFrame, jump_thresh: float=1.0):

        x = df['e'].to_numpy()
        y = df['n'].to_numpy()
        t = df['time_s'].to_numpy()# / 1000
        t = np.diff(t)
        ds = np.zeros((len(x)))

        x1 = x[:-1]
        y1 = y[:-1]
        x2 = x[1:]
        y2 = y[1:]

        d = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        ds[1:] = d
        df['dist'] = ds

        # Calculate speed
        s = np.where(t != 0, ds[1:] / t, np.nan)
        s = np.append(s[0], s)
        s = pd.Series(s).fillna(method='ffill').fillna(method='bfill').to_numpy()


        # Assume constant speed for nan's. Need to interpolate.
        s = pd.Series(s, index=df.index)
        s.replace(0, np.nan, inplace=True)  # Replace 0 with NaN for interpolation
        s.interpolate(method='linear', inplace=True)

        # Accumulate distance
        ds = np.cumsum(ds)

        # Store
        df['speed_ms'] = np.around(s, 1)
        df['trk_dist'] = ds
        # df['speed_ms'] = df['speed_ms'].fillna(0)



        return df
    
    # ======================================================================
    def _doUnitConversion(self, df: pd.DataFrame):

        #####################
        # Convert depth units
        def decode_varint(data):
            """Decodes a varint from a byte string.

            Args:
                data: A byte string containing the varint.

            Returns:
                The decoded integer value.
            """
            result = 0
            shift = 0
            for byte in data:
                result |= (byte & 0x7f) << shift
                shift += 7
                if not (byte & 0x80):  # If MSB is 0, end of varint
                    break
            return result

        # Convert varint values (Assume they are in thousandths of a unit)
        cols_to_convert = ['bottom_depth', 'drawn_bottom_depth', 'last_sample_depth']
        for col in cols_to_convert:

            if col in df.columns:
                df[col] = df[col].apply(lambda x: decode_varint(x) if isinstance(x, bytes) else x)
                df[col] = df[col].astype(float) / 1000.0
                df[col] /= 3.2808399

        df['first_sample_depth'] /= 3.2808399



        ##############
        # Convert time
        # Garmin uses 0000 December 31, 1989 as start time

        start_date = datetime(1989, 12, 31, 0, 0, 0, tzinfo=timezone.utc)

        custom_unix_time = start_date + timedelta(seconds=self.file_header['date_time_of_recording'])

        custom_unix_time = custom_unix_time.timestamp()

        df['recording_time_ms'] /= 1000  # Convert ms to seconds

        custom_unix_time = custom_unix_time + (df['recording_time_ms'])

        df['caltime'] = pd.to_datetime(custom_unix_time, unit='s', utc=True)

        df['date'] = df['caltime'].dt.date
        df['time'] = df['caltime'].dt.time # Time in utc

        df.drop(columns=['caltime'], inplace=True)


        ##################################
        # Calculate latitude and longitude
        # df['lat'] = df['scposn_lat'] * 360 / (1<<32)
        # df['lon'] = df['scposn_lon'] * 360 / (1<<32)

        # df['lat'] = df['scposn_lat'].astype('float64') * (180.0 / (2**31))
        # df['lon'] = df['scposn_lon'].astype('float64') * (180.0 / (2**31))

        df['lat'] = df['scposn_lat'] * 360 / (1 << 32)
        df['lon'] = df['scposn_lon'] * 360 / (1 << 32)

        df['lat'] = df['lat'].apply(lambda x: x - 360 if x > 180 else x)
        df['lon'] = df['lon'].apply(lambda x: x - 360 if x > 180 else x)

        # Do filtering
        df = filterGPS(df)




        # print('\n\nConverted lat lon:')
        # print(df[['lat', 'lon']].describe())

        # import matplotlib.pyplot as plt

        # plt.scatter(df['scposn_lon'], df['scposn_lat'], s=1, alpha=0.5)
        # plt.title('Vessel Track')
        # plt.xlabel('Longitude')
        # plt.ylabel('Latitude')
        # plt.grid(True)
        # plt.show()

        # plt.figure(figsize=(10, 6))
        # plt.scatter(df['lon'], df['lat'], s=1, alpha=0.5, label='Cleaned')
        # plt.title('GPS Track After Percentile Filtering')
        # plt.xlabel('Longitude')
        # plt.ylabel('Latitude')
        # plt.grid(True)
        # plt.legend()
        # plt.show()



        # Determine epsg code
        self.humDat['epsg'] = "EPSG:"+str(int(float(self._convert_wgs_to_utm(df['lon'][0], df['lat'][0]))))
        self.humDat['wgs'] = "EPSG:4326"

        # Configure re-projection function
        self.trans = pyproj.Proj(self.humDat['epsg'])

        # Reproject lat/lon to UTM zone
        e, n = self.trans(df['lon'], df['lat'])
        df['e'] = e
        df['n'] = n


        #########################
        # Calculate COG (heading)
        ## Garmin does not appear to store heading....
        heading = self._getCOG(df)
        # self._getBearing() returns n-1 values because last ping can't
        ## have a COG value.  We will duplicate the last COG value and use it for
        ## the last ping.
        last = heading[-1]
        heading = np.append(heading, last)
        df['instr_heading'] = heading # Store COG in sDF

        # Replace 0 with NaN, interpolate, then fill any remaining NaN (e.g., with nearest valid value)
        df['instr_heading'] = np.around(df['instr_heading'].replace(0, np.nan).interpolate().bfill().ffill(), 1)

        # Add transect number (for aoi processing)
        df['transect'] = 0

        # Calculate pixel size [m]  *** ....MAYBE.... ***
        df['pixM'] = (df['last_sample_depth'] - df['first_sample_depth']) / df['sample_cnt']

        return df
    
    #===========================================
    def _getCOG(self,
                df,
                lon = 'lon',
                lat = 'lat'):
        '''
        Calculates course over ground (COG) from a set of coordinates.  Since the
        last coordinate pair cannot have a COG value, the length of the returned
        array is len(n-1) where n == len(df).

        ----------
        Parameters
        ----------
        df : DataFrame
            DESCRIPTION - Pandas dataframe with geographic coordinates of sonar
                          records.
        lon : str : [Default='lons']
            DESCRIPTION - DataFrame column name for longitude coordinates.
        lat : str : [Default='lats']
            DESCRIPTION - DataFrame column name for latitude coordinates.

        ----------------------------
        Required Pre-processing step
        ----------------------------
        Called from self._interpTrack()

        -------
        Returns
        -------
        Numpy array of COG values.

        --------------------
        Next Processing Step
        --------------------
        Return to self._interpTrack()
        '''
        # COG calculation will be calculated on numpy arrays for speed.  Since
        ## COG is calculated from one point to another (pntA -> pntB), we need
        ## to store pntA values, beginning with the first value and ending at
        ## second to last value, in one array and pntB values, beginning at second
        ## value and ending at last value, in another array.  We can then use
        ## vector algebra to efficiently calculate COG.

        # Prepare pntA values [0:n-1]
        lonA = df[lon].to_numpy() # Store longitude coordinates in numpy array
        latA = df[lat].to_numpy() # Store longitude coordinates in numpy array
        lonA = lonA[:-1] # Omit last coordinate
        latA = latA[:-1] # Omit last coordinate
        pntA = [lonA,latA] # Store in array of arrays

        # Prepare pntB values [0+1:n]
        lonB = df[lon].to_numpy() # Store longitude coordinates in numpy array
        latB = df[lat].to_numpy() # Store longitude coordinates in numpy array
        lonB = lonB[1:] # Omit first coordinate
        latB = latB[1:] # Omit first coordinate
        pntB = [lonB,latB] # Store in array of arrays

        # Convert latitude values into radians
        lat1 = np.deg2rad(pntA[1])
        lat2 = np.deg2rad(pntB[1])

        diffLong = np.deg2rad(pntB[0] - pntA[0]) # Calculate difference in longitude then convert to degrees
        bearing = np.arctan2(np.sin(diffLong) * np.cos(lat2), np.cos(lat1) * np.sin(lat2) - (np.sin(lat1) * np.cos(lat2) * np.cos(diffLong))) # Calculate bearing in radians

        db = np.degrees(bearing) # Convert radians to degrees
        db = (db + 360) % 360 # Ensure degrees in range 0-360

        return db
    
    # ======================================================================
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
    
    
    ### Format to PINGMapper ###
    # ======================================================================
    def _recalcRecordNum(self):

        df = self.header_dat

        # Reset index and recalculate record num
        ## Record num is unique for each ping across all sonar beams
        df = df.reset_index(drop=True)
        df['record_num'] = df.index

        self.header_dat = df
        return
    

    # ======================================================================
    def _splitBeamsToCSV(self):

        beam_set = {}

        # Dictionary to store necessary attributes for PING-Mapper
        self.beamMeta = beamMeta = {}

        # Get df
        df = self.header_dat

        # Get all unique beam and port_star_id combinations
        try:
            dfBeams = df.drop_duplicates(subset=['channel_id', 'F'])
        except:
            dfBeams = df.drop_duplicates(subset=['channel_id'])
            dfBeams['port_star_id'] = np.nan

        # # Find channel_id for nan port_star_id
        # nanPortStar = dfBeams[dfBeams['port_star_id'].isna()]['channel_id'].unique()

        # if len(nanPortStar) > 1:
        #     beam_set[nanPortStar.min()] = ('ds_hifreq', 1) # Default beam
        #     beam_set[nanPortStar.max()] = ('ds_vhifreq', 4) # Default beam
        # else:
        #     beam_set[nanPortStar[0]] = ('ds_hifreq', 1)

        # try:
        #     # Get channel_id from dfBeams for port is port_star_id==60
        #     port_chan_id = dfBeams[dfBeams['port_star_id'] == 60]['channel_id'].iloc[0]
        #     star_chan_id = dfBeams[dfBeams['port_star_id'] == -60]['channel_id'].iloc[0]
        #     beam_set[port_chan_id] = ('ss_port', 2)
        #     beam_set[star_chan_id] = ('ss_star', 3)
        # except:
        #     pass

        # Nievely assign beams
        # if there are four beams, assume:
        ## min value is 2d
        ## second is down image
        ## third is port
        ## max is star
        print('\n\nBeams Available:')
        print(dfBeams['channel_id'])
        if len(dfBeams) == 4:
            # min, port, star, di = sorted(dfBeams['channel_id'].unique())
            min, di, port, star = sorted(dfBeams['channel_id'].unique())
            beam_set[min] = ('ds_hifreq', 1)
            beam_set[port] = ('ss_port', 2)
            beam_set[star] = ('ss_star', 3)
            beam_set[di] = ('ds_vhifreq', 4)
        
        # If only 1 beam, assign to high freq
        elif len(dfBeams) == 1:
            beam_set[dfBeams['channel_id'][0]] = ('ds_hifreq', 1)

        # Unknown return error
        else:
            print('\n\nERROR!')
            print('Unknown beam ids:')
            print(dfBeams['channel_id'].unique())
            sys.exit()


        # Iterate each beam
        for beam, group in df.groupby('channel_id'):
            meta = {}
            
            # Get Garmin beam to Humminbird beam
            humBeamName, humBeamint = beam_set[beam]
            humBeam = 'B00'+str(humBeamint)
            meta['beamName'] = humBeamName
            meta['beam'] = humBeam
            group['beam'] = humBeamint

            # # Set pixM based on side scan
            # if humBeamint == 2 or humBeamint == 3:
            #     self.pixM = group['pixM'].iloc[0]

            # Store sonFile
            meta['sonFile'] = self.sonFile

            # Drop columns
            cols2Drop = ['magic_number']
            cols = group.columns
            cols2Drop += [c for c in cols if 'fp' in c]
            cols2Drop += [c for c in cols if 'SP' in c]
            cols2Drop += [c for c in cols if 'su' in c]
            group.drop(columns=cols2Drop, inplace=True)

            # Add chunk_id
            group = self._getChunkID(group)

            # Save csv
            outCSV = '{}_{}_meta.csv'.format(humBeam, meta['beamName'])
            outCSV = os.path.join(self.metaDir, outCSV)
            group.to_csv(outCSV, index=False)

            meta['metaCSV'] = outCSV

            # Store the beams metadata
            beamMeta[humBeam] = meta

        return
    

    # ======================================================================
    def _getChunkID(self, df: pd.DataFrame):

        df.reset_index(drop=True, inplace=True)

        df['chunk_id'] = int(-1)

        chunk = 0
        start_idx = chunk
        end_idx = self.nchunk

        while start_idx < len(df):

            df.iloc[start_idx:end_idx, df.columns.get_loc('chunk_id')] = int(chunk)

            chunk += 1
            start_idx = end_idx
            end_idx += self.nchunk

        # Update last chunk if too small (for rectification)
        lastChunk = df[df['chunk_id'] == chunk]
        if len(lastChunk) <= self.nchunk/2:
            df.loc[df['chunk_id'] == chunk, 'chunk_id'] = chunk-1


        return df





    # ======================================================================
    def __str__(self):
        '''
        Generic print function to print contents of sonObj.
        '''
        output = "Lowrance Class Contents"
        output += '\n\t'
        output += self.__repr__()
        temp = vars(self)
        for item in temp:
            output += '\n\t'
            output += "{} : {}".format(item, temp[item])
        return output