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
from io import BytesIO

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

        # Parse channel_information_array (field 6) for per-channel metadata.
        self.channel_info = self._parseChannelInformation(firstHeadBytes)


        return

    # ======================================================================
    def _read_varuint(self, file_obj):
        """Read protobuf-style VarUInt32/VarUInt64 from a file-like object."""
        result = 0
        shift = 0

        while True:
            b = file_obj.read(1)
            if not b:
                raise EOFError('Unexpected EOF while reading varuint.')

            byte = b[0]
            result |= (byte & 0x7F) << shift

            if not (byte & 0x80):
                return result

            shift += 7
            if shift > 63:
                raise ValueError('Invalid varuint: too many continuation bytes.')

    # ======================================================================
    def _decode_varuint_bytes(self, data: bytes):
        stream = BytesIO(data)
        return self._read_varuint(stream)

    # ======================================================================
    def _decode_varint32_bytes(self, data: bytes):
        """Decode zig-zag encoded VarInt32 from bytes."""
        u = self._decode_varuint_bytes(data)
        return (u >> 1) ^ -(u & 1)

    # ======================================================================
    def _decode_depth_varint(self, data: bytes):
        """Decode Garmin depth-like varints.

        Some RSD files appear to store these values as unsigned varints while
        others align with zig-zag VarInt32 semantics from the reverse-engineered
        spec. If zig-zag yields a negative depth, fall back to unsigned.
        """
        u = self._decode_varuint_bytes(data)
        z = self._decode_varint32_bytes(data)

        candidates = [v for v in (u, z) if v >= 0]
        if not candidates:
            return u

        # Keep depth-like values within a broad physical bound (<= 5 km in mm)
        # when possible, otherwise pick the smaller non-negative decoding.
        plausible = [v for v in candidates if v <= 5_000_000]
        if plausible:
            return min(plausible)

        # Some files appear to serialize these fields as fixed-width integers.
        if len(data) in (2, 4):
            le = int.from_bytes(data, 'little', signed=False)
            if le <= 5_000_000:
                return le

        return min(candidates)

    # ======================================================================
    def _parse_var_struct_payload(self, payload: bytes):
        """Parse a variable-structure payload into (field_number, raw_value_bytes)."""
        stream = BytesIO(payload)
        fields = []

        try:
            field_cnt = self._read_varuint(stream)
        except EOFError:
            return fields

        for _ in range(field_cnt):
            key = self._read_varuint(stream)
            field_num = key >> 3
            val_len = key & 0x07

            if val_len == 7:
                val_len = self._read_varuint(stream)

            val = stream.read(val_len)
            fields.append((field_num, val))

        return fields

    # ======================================================================
    def _parseChannelInformation(self, chan_info_offset: int):
        """Parse header field 6 channel_information_array for each channel."""
        out = []

        with open(self.sonFile, 'rb') as file:
            file.seek(chan_info_offset)

            key = self._read_varuint(file)
            if key != 55:  # field 6 with length marker 7 -> 6<<3 + 7
                return out

            payload_len = self._read_varuint(file)
            payload = file.read(payload_len)

        stream = BytesIO(payload)

        try:
            chan_cnt = self._read_varuint(stream)
        except EOFError:
            return out

        for _ in range(chan_cnt):
            ch_payload = b''
            ch = {
                'channel_id': np.nan,
                'first_chunk_offset': np.nan,
                'transducer_port': np.nan,
                'start_freq_hz': np.nan,
                'end_freq_hz': np.nan,
                'channel_capabilities': np.nan,
            }

            # Each channel-info element is a variable structure.
            # Parse directly from stream by reading element field count and fields.
            try:
                elem_field_cnt = self._read_varuint(stream)
            except EOFError:
                break

            for _ in range(elem_field_cnt):
                fkey = self._read_varuint(stream)
                fnum = fkey >> 3
                flen = fkey & 0x07
                if flen == 7:
                    flen = self._read_varuint(stream)
                fval = stream.read(flen)

                if fnum == 0:
                    # data_info: varray of DataInfo structure(s)
                    info_stream = BytesIO(fval)
                    n = self._read_varuint(info_stream)
                    if n > 0:
                        for _ in range(n):
                            item_len = self._read_varuint(info_stream)
                            item = info_stream.read(item_len)
                            if item_len > 0:
                                ch['channel_id'] = self._decode_varuint_bytes(item)
                elif fnum == 1:
                    # first_chunk_offset: 8-byte little-endian ulong
                    if len(fval) == 8:
                        ch['first_chunk_offset'] = int.from_bytes(fval, 'little', signed=False)
                elif fnum == 2:
                    # prop_chan_info: varray of DpsChannelInformation struct(s)
                    dps_stream = BytesIO(fval)
                    dps_n = self._read_varuint(dps_stream)

                    for _ in range(dps_n):
                        dps_len = self._read_varuint(dps_stream)
                        dps_payload = dps_stream.read(dps_len)
                        dps_fields = self._parse_var_struct_payload(dps_payload)

                        for dnum, dval in dps_fields:
                            if dnum == 0:
                                ch['transducer_port'] = self._decode_varuint_bytes(dval)
                            elif dnum == 1:
                                tf_fields = self._parse_var_struct_payload(dval)
                                for tnum, tval in tf_fields:
                                    if tnum == 1:
                                        ch['start_freq_hz'] = self._decode_varuint_bytes(tval)
                                    elif tnum == 2:
                                        ch['end_freq_hz'] = self._decode_varuint_bytes(tval)
                            elif dnum == 2:
                                ch['channel_capabilities'] = self._decode_varuint_bytes(dval)

            out.append(ch)

        return out

    # ======================================================================
    def _parseBeamInfoPayload(self, payload: bytes):
        """Parse optional beam_info structure payload and return decoded values."""
        out = {}
        fields = self._parse_var_struct_payload(payload)

        for field_num, value in fields:
            if field_num == 0 and len(value) > 0:
                out['port_star_beam_angle'] = int(value[0])
            elif field_num == 1 and len(value) > 0:
                out['fore_aft_beam_angle'] = int(value[0])
            elif field_num == 2 and len(value) > 0:
                out['port_star_elem_angle'] = int(value[0])
            elif field_num == 3 and len(value) > 0:
                out['fore_aft_elem_angle'] = int(value[0])
            elif field_num == 5:
                # StructUnknown2 field 0 is strongly associated with port/star sign.
                su2_fields = self._parse_var_struct_payload(value)
                for su2_num, su2_val in su2_fields:
                    if su2_num == 0 and len(su2_val) == 4:
                        out['port_star_id'] = float(np.frombuffer(su2_val, dtype='<f4')[0])

        return out
    
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

        # Decode variable record body by field key/length instead of fixed count assumptions.
        rb_field_cnt = self._read_varuint(file)
        out_dict['record_body_fcnt'] = rb_field_cnt

        for _ in range(rb_field_cnt):
            key = self._read_varuint(file)
            field_num = key >> 3
            value_len = key & 0x07
            if value_len == 7:
                value_len = self._read_varuint(file)

            raw = file.read(value_len)

            if field_num == 0:
                out_dict['channel_id_1'] = self._decode_varuint_bytes(raw)
            elif field_num == 1:
                out_dict['bottom_depth'] = self._decode_depth_varint(raw)
            elif field_num == 2:
                out_dict['drawn_bottom_depth'] = self._decode_depth_varint(raw)
            elif field_num == 3:
                out_dict['first_sample_depth'] = self._decode_depth_varint(raw)
            elif field_num == 4:
                out_dict['last_sample_depth'] = self._decode_depth_varint(raw)
            elif field_num == 5 and len(raw) > 0:
                out_dict['gain'] = int(raw[0])
            elif field_num == 6:
                out_dict['sample_status'] = self._decode_varuint_bytes(raw)
            elif field_num == 7:
                out_dict['sample_cnt'] = int.from_bytes(raw, 'little', signed=False)
            elif field_num == 8 and len(raw) > 0:
                out_dict['shade_avail'] = int(raw[0])
            elif field_num == 9:
                out_dict['scposn_lat'] = int.from_bytes(raw, 'little', signed=True)
            elif field_num == 10:
                out_dict['scposn_lon'] = int.from_bytes(raw, 'little', signed=True)
            elif field_num == 11 and len(raw) == 4:
                out_dict['water_temp'] = float(np.frombuffer(raw, dtype='<f4')[0])
            elif field_num == 12:
                out_dict['beam'] = self._decode_varuint_bytes(raw)
            elif field_num == 13:
                out_dict.update(self._parseBeamInfoPayload(raw))
            elif field_num == 14:
                out_dict['interrogation_id'] = self._decode_varuint_bytes(raw)

            


        # Next ping header is from current position + ping_cnt
        # next_ping = file.tell() + out_dict['packet_size']
        next_ping = i + pingHeaderLen + out_dict['data_size'] + 12 #12 for magic number & crc

        out_dict['index'] = i

        sample_cnt = out_dict.get('sample_cnt', 0)
        out_dict['son_offset'] = (out_dict['data_size']) - (sample_cnt*2) + self.pingHeaderLen

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
        s = np.full_like(t, np.nan, dtype='float64')
        np.divide(ds[1:], t, out=s, where=t != 0)
        s = np.append(s[0], s)
        s = pd.Series(s).ffill().bfill().to_numpy()


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

        # Garmin depth-like fields are stored in millimeters. Convert to meters.
        cols_to_convert = ['bottom_depth', 'drawn_bottom_depth', 'last_sample_depth']
        for col in cols_to_convert:

            if col in df.columns:
                df[col] = df[col].apply(lambda x: decode_varint(x) if isinstance(x, bytes) else x)
                # Garmin uses sentinel-style out-of-range values in some files.
                df[col] = df[col].where((df[col] >= 0) & (df[col] < 0xFFFFFF00), np.nan)
                df[col] = df[col].astype(float) / 1000.0

        if 'first_sample_depth' in df.columns:
            df['first_sample_depth'] = df['first_sample_depth'].where(
                (df['first_sample_depth'] >= 0) & (df['first_sample_depth'] < 0xFFFFFF00),
                np.nan,
            )
            df['first_sample_depth'] = df['first_sample_depth'].astype(float) / 1000.0



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

        # Prefer explicit in-record beam values and use port_star_id sign to disambiguate sidescan.
        beam_name_by_id = {
            1: 'ds_hifreq',
            2: 'ss_port',
            3: 'ss_star',
            4: 'ds_vhifreq',
        }

        chan_ids = sorted(df['channel_id'].dropna().unique())

        for chan_id in chan_ids:
            g = df[df['channel_id'] == chan_id]

            beam_id = np.nan
            if 'beam' in g.columns:
                b = g['beam'].dropna()
                b = b[b.isin([1, 2, 3, 4])]
                if len(b) > 0:
                    beam_id = int(b.mode().iloc[0])

            if np.isnan(beam_id) and isinstance(getattr(self, 'channel_info', None), list):
                # Conservative fallback from channel header metadata when beam is not present.
                match = [c for c in self.channel_info if c.get('channel_id', np.nan) == chan_id]
                if len(match) == 1:
                    cinfo = match[0]
                    sf = cinfo.get('start_freq_hz', np.nan)
                    ef = cinfo.get('end_freq_hz', np.nan)
                    if np.isfinite(sf) and np.isfinite(ef):
                        fkhz = (sf + ef) / 2000.0
                        if fkhz >= 700:
                            beam_id = 4
                        elif fkhz >= 300:
                            beam_id = 2
                        else:
                            beam_id = 1

            if np.isnan(beam_id):
                beam_set[chan_id] = ('unknown', -1)
            else:
                beam_set[chan_id] = (beam_name_by_id.get(int(beam_id), 'unknown'), int(beam_id))

        # Use sign of port_star_id to force sidescan orientation when available.
        if 'port_star_id' in df.columns:
            for chan_id in chan_ids:
                g = df[df['channel_id'] == chan_id]
                ps = g['port_star_id'].dropna()
                if len(ps) == 0:
                    continue

                ps_med = float(np.median(ps))
                if ps_med < 0:
                    beam_set[chan_id] = ('ss_port', 2)
                elif ps_med > 0:
                    beam_set[chan_id] = ('ss_star', 3)

        # Final fallback if still unresolved.
        unresolved = [k for k, v in beam_set.items() if v[1] == -1]
        if len(unresolved) > 0:
            ordered = sorted(chan_ids)
            for idx, chan_id in enumerate(ordered):
                if chan_id not in unresolved:
                    continue
                if idx == 0:
                    beam_set[chan_id] = ('ds_hifreq', 1)
                elif idx == 1:
                    beam_set[chan_id] = ('ss_port', 2)
                elif idx == 2:
                    beam_set[chan_id] = ('ss_star', 3)
                else:
                    beam_set[chan_id] = ('ds_vhifreq', 4)


        # Iterate each beam
        for beam, group in df.groupby('channel_id'):
            meta = {}
            
            # Get Garmin beam to Humminbird beam
            humBeamName, humBeamint = beam_set[beam]
            humBeam = 'B00'+str(humBeamint)
            meta['beamName'] = humBeamName
            meta['beam'] = humBeam
            group = group.copy()
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