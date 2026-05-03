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

import json
import math
import os, sys
import numpy as np
import pandas as pd
from datetime import datetime, timezone, timedelta
from io import BytesIO

try:
    import pyproj
except ImportError:  # Garmin RSD parsing and waterfall validation can run without projection support.
    pyproj = None

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
        if not data:
            return np.nan
        stream = BytesIO(data)
        try:
            return self._read_varuint(stream)
        except (EOFError, ValueError):
            return np.nan

    # ======================================================================
    def _decode_varint32_bytes(self, data: bytes):
        """Decode zig-zag encoded VarInt32 from bytes."""
        if not data:
            return np.nan
        u = self._decode_varuint_bytes(data)
        if pd.isna(u):
            return np.nan
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

        if pd.isna(u):
            if len(data) in (1, 2, 4):
                le = int.from_bytes(data, 'little', signed=False)
                if le < 0xFFFFFF00:
                    return le
            return np.nan

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
        except (EOFError, ValueError):
            return fields

        for _ in range(field_cnt):
            try:
                key = self._read_varuint(stream)
            except (EOFError, ValueError):
                break

            field_num = key >> 3
            val_len = key & 0x07

            if val_len == 7:
                try:
                    val_len = self._read_varuint(stream)
                except (EOFError, ValueError):
                    break

            val = stream.read(val_len)
            if len(val) < val_len:
                break

            fields.append((field_num, val))

        return fields

    # ======================================================================
    def _read_var_struct_fields(self, file_obj, max_end: int=None):
        """Read a variable-structure from a file and return raw field values."""
        fields = []
        field_cnt = self._read_varuint(file_obj)

        for _ in range(field_cnt):
            if max_end is not None and file_obj.tell() >= max_end:
                break

            key = self._read_varuint(file_obj)
            field_num = key >> 3
            value_len = key & 0x07

            if value_len == 7:
                value_len = self._read_varuint(file_obj)

            if max_end is not None and file_obj.tell() + value_len > max_end:
                break

            raw = file_obj.read(value_len)
            if len(raw) < value_len:
                raise EOFError('Unexpected EOF while reading variable structure field.')

            fields.append((field_num, raw))

        return fields

    # ======================================================================
    def _parseChannelInformation(self, chan_info_offset: int):
        """Parse header field 6 channel_information_array for each channel."""
        out = []

        with open(self.sonFile, 'rb') as file:
            file.seek(chan_info_offset)

            try:
                key = self._read_varuint(file)
            except (EOFError, ValueError):
                return out

            if key != 55:  # field 6 with length marker 7 -> 6<<3 + 7
                return out

            try:
                payload_len = self._read_varuint(file)
            except (EOFError, ValueError):
                return out

            payload = file.read(payload_len)
            if len(payload) < payload_len:
                return out

        stream = BytesIO(payload)

        try:
            chan_cnt = self._read_varuint(stream)
        except (EOFError, ValueError):
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
            except (EOFError, ValueError):
                break

            for _ in range(elem_field_cnt):
                try:
                    fkey = self._read_varuint(stream)
                except (EOFError, ValueError):
                    break

                fnum = fkey >> 3
                flen = fkey & 0x07
                if flen == 7:
                    try:
                        flen = self._read_varuint(stream)
                    except (EOFError, ValueError):
                        break
                fval = stream.read(flen)
                if len(fval) < flen:
                    break

                if fnum == 0:
                    # data_info: varray of DataInfo structure(s)
                    info_stream = BytesIO(fval)
                    try:
                        n = self._read_varuint(info_stream)
                    except (EOFError, ValueError):
                        n = 0
                    if n > 0:
                        for _ in range(n):
                            try:
                                item_len = self._read_varuint(info_stream)
                            except (EOFError, ValueError):
                                break
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
                    try:
                        dps_n = self._read_varuint(dps_stream)
                    except (EOFError, ValueError):
                        dps_n = 0

                    for _ in range(dps_n):
                        try:
                            dps_len = self._read_varuint(dps_stream)
                        except (EOFError, ValueError):
                            break
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

        # Store contents in list
        header_dat_all = []

        # Decode ping header
        with open(self.sonFile, 'rb') as file:
            while i < file_len:

                # Get header data at offset i
                header_dat, cpos = self._getPingHeader(file, i)

                if header_dat:
                    header_dat_all.append(header_dat)

                i = cpos

        # Convert to dataframe
        df = pd.DataFrame.from_dict(header_dat_all)
        if len(df) == 0:
            self.header_dat = df
            return

        # Convert fields
        df = self._doUnitConversion(df)

        # Do column name conversions to PINGMapper units
        df.rename(columns=self.garCols2PM, inplace=True)

        # Calculate speed & track distance (based on coords and time)
        df = self._calcSpeedTrkDist(df)
        self._recomputeTrackSpeedFromWgs(df)

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

        # Move to offset
        file.seek(i)

        try:
            header_fields = self._read_var_struct_fields(file)
        except (EOFError, ValueError):
            return False, self.file_len

        out_dict = {}
        for field_num, raw in header_fields:
            if field_num == 0 and len(raw) == 4:
                out_dict['magic_number'] = int.from_bytes(raw, 'little', signed=False)
            elif field_num == 1:
                out_dict.update(self._parseStateDataPayload(raw))
            elif field_num == 2 and len(raw) == 4:
                out_dict['sequence_cnt'] = int.from_bytes(raw, 'little', signed=False)
            elif field_num == 3 and len(raw) == 4:
                out_dict['data_crc'] = int.from_bytes(raw, 'little', signed=False)
            elif field_num == 4 and len(raw) == 2:
                out_dict['data_size'] = int.from_bytes(raw, 'little', signed=False)
            elif field_num == 5 and len(raw) == 4:
                out_dict['recording_time_ms'] = int.from_bytes(raw, 'little', signed=False)

        if out_dict.get('magic_number') != self.magicNum or 'data_size' not in out_dict:
            return False, self._find_next_record(file, i + 1)

        # Skip record header CRC.
        header_crc = file.read(4)
        if len(header_crc) < 4:
            return False, self.file_len
        out_dict['record_crc'] = int.from_bytes(header_crc, 'little', signed=False)
        pingHeaderLen = file.tell() - i
        out_dict['ping_header_len'] = pingHeaderLen

        # Check if there is a record body
        if out_dict.get('state') != 2 or out_dict['data_size'] == 0: # no record body
            next_ping = i + pingHeaderLen + 12
            return False, self._align_next_record(file, next_ping)
        
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
        try:
            rb_field_cnt = self._read_varuint(file)
        except (EOFError, ValueError):
            return False, self.file_len

        out_dict['record_body_fcnt'] = rb_field_cnt
        record_body_start = i + pingHeaderLen
        record_body_end = record_body_start + out_dict['data_size']

        for _ in range(rb_field_cnt):
            if file.tell() >= record_body_end:
                break

            try:
                key = self._read_varuint(file)
            except (EOFError, ValueError):
                break

            field_num = key >> 3
            value_len = key & 0x07
            if value_len == 7:
                try:
                    value_len = self._read_varuint(file)
                except (EOFError, ValueError):
                    break

            if value_len < 0 or file.tell() + value_len > record_body_end:
                break

            raw = file.read(value_len)
            if len(raw) < value_len:
                break

            if field_num == 0 and len(raw) > 0:
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
            elif field_num == 6 and len(raw) > 0:
                out_dict['sample_status'] = self._decode_varuint_bytes(raw)
            elif field_num == 7 and len(raw) > 0:
                out_dict['sample_cnt'] = int.from_bytes(raw, 'little', signed=False)
            elif field_num == 8 and len(raw) > 0:
                out_dict['shade_avail'] = int(raw[0])
            elif field_num == 9 and len(raw) > 0:
                out_dict['scposn_lat'] = int.from_bytes(raw, 'little', signed=True)
            elif field_num == 10 and len(raw) > 0:
                out_dict['scposn_lon'] = int.from_bytes(raw, 'little', signed=True)
            elif field_num == 11 and len(raw) == 4:
                out_dict['water_temp'] = float(np.frombuffer(raw, dtype='<f4')[0])
            elif field_num == 12 and len(raw) > 0:
                out_dict['beam'] = self._decode_varuint_bytes(raw)
            elif field_num == 13:
                out_dict.update(self._parseBeamInfoPayload(raw))
            elif field_num == 14 and len(raw) > 0:
                out_dict['interrogation_id'] = self._decode_varuint_bytes(raw)

            


        # Next ping header is from current position + ping_cnt
        # next_ping = file.tell() + out_dict['packet_size']
        next_ping = i + pingHeaderLen + out_dict['data_size'] + 12 #12 for trailer magic, chunk size & crc

        out_dict['index'] = i

        sample_cnt = out_dict.get('sample_cnt', 0)
        if pd.isna(sample_cnt):
            sample_cnt = 0
        out_dict['son_offset'] = (out_dict['data_size']) - (sample_cnt*2) + pingHeaderLen

        # out_dict['son_offset'] = record_body_header_len+1
 
        return out_dict, self._align_next_record(file, next_ping)

    # ======================================================================
    def _parseStateDataPayload(self, payload: bytes):
        out = {}
        fields = self._parse_var_struct_payload(payload)

        for field_num, raw in fields:
            if field_num == 0 and len(raw) > 0:
                out['state'] = self._decode_varuint_bytes(raw)
            elif field_num == 1:
                values = self._parse_varuint_varray(raw)
                if values:
                    out['channel_id'] = values[0]

        return out

    # ======================================================================
    def _parse_varuint_varray(self, payload: bytes):
        """Parse Garmin varray[VarUInt32] values."""
        stream = BytesIO(payload)
        values = []

        try:
            count = self._read_varuint(stream)
            total_len = self._read_varuint(stream)
        except (EOFError, ValueError):
            return values

        end_pos = min(len(payload), stream.tell() + total_len)
        for _ in range(count):
            if stream.tell() >= end_pos:
                break
            try:
                values.append(self._read_varuint(stream))
            except (EOFError, ValueError):
                break

        return values

    # ======================================================================
    def _find_next_record(self, file, start_pos: int):
        magic = self.magicNum.to_bytes(4, 'little')
        file.seek(start_pos)
        data = file.read(max(0, self.file_len - start_pos))
        pos = data.find(magic)
        if pos < 0:
            return self.file_len
        return start_pos + pos - 2

    # ======================================================================
    def _align_next_record(self, file, expected_pos: int):
        if expected_pos >= self.file_len:
            return self.file_len

        file.seek(expected_pos)
        b = file.read(6)
        if len(b) >= 6 and b[0] == 6 and int.from_bytes(b[2:6], 'little') == self.magicNum:
            return expected_pos

        return self._find_next_record(file, max(self.headBytes, expected_pos - 16))

    # ======================================================================
    def extract_raw_sample_arrays(self, df: pd.DataFrame=None):
        """Return raw uint16 sonar samples per ping grouped by channel_id.

        Each read uses the decoded record index, data_size, sample_cnt, and
        son_offset fields. Pings with incomplete or inconsistent sample blocks
        are skipped instead of raising.
        """
        if df is None:
            df = self.header_dat

        samples_by_channel = {}

        with open(self.sonFile, 'rb') as file:
            for _, row in df.iterrows():
                try:
                    channel_id = row['channel_id']
                    index = int(row['index'])
                    data_size = int(row['data_size'])
                    sample_cnt = int(row['sample_cnt'] if 'sample_cnt' in row else row['ping_cnt'])
                    son_offset = int(row['son_offset'])
                    ping_header_len = int(row.get('ping_header_len', self.pingHeaderLen))
                except (KeyError, TypeError, ValueError):
                    continue

                sample_bytes = sample_cnt * 2
                if pd.isna(channel_id) or index < 0 or data_size <= 0 or sample_cnt <= 0:
                    continue
                if son_offset < ping_header_len:
                    continue
                if son_offset + sample_bytes > ping_header_len + data_size:
                    continue

                file.seek(index + son_offset)
                raw = file.read(sample_bytes)
                if len(raw) != sample_bytes:
                    continue

                arr = np.frombuffer(raw, dtype='<u2').copy()
                samples_by_channel.setdefault(int(channel_id), []).append(arr)

        return samples_by_channel

    # ======================================================================
    def write_channel_waterfall_pngs(self, out_dir: str, df: pd.DataFrame=None,
                                     prefix: str=None, width: int=None):
        """Write one Garmin-style waterfall PNG per channel and return paths."""
        from PIL import Image

        os.makedirs(out_dir, exist_ok=True)
        samples_by_channel = self.extract_raw_sample_arrays(df)
        out_paths = {}

        for channel_id, pings in samples_by_channel.items():
            if not pings:
                continue

            max_len = width or max(len(p) for p in pings)
            img_arr = np.zeros((len(pings), max_len), dtype=np.uint16)

            for row_idx, ping in enumerate(pings):
                n = min(len(ping), max_len)
                img_arr[row_idx, :n] = ping[:n]

            scaled = self._scale_samples_for_waterfall(img_arr)
            image = Image.fromarray(scaled, mode='P')
            image.putpalette(self._garmin_waterfall_palette())

            stem = prefix or os.path.splitext(os.path.basename(self.sonFile))[0]
            safe_stem = ''.join(c if c.isalnum() or c in ('-', '_') else '_' for c in stem)
            out_path = os.path.join(out_dir, f'{safe_stem}_channel_{channel_id}.png')
            image.save(out_path)
            out_paths[channel_id] = out_path

        return out_paths

    # ======================================================================
    def write_sonar_data_player_project(self, out_dir: str, include_pngs: bool=True,
                                        prefix: str=None):
        """Write a SonarDataPlayer processed project for this Garmin RSD file.

        The project contains ping telemetry CSV, synchronized frame metadata,
        a raw uint16 little-endian sample blob, and optional per-channel PNG
        previews. The parser state is initialized automatically when needed.
        """
        os.makedirs(out_dir, exist_ok=True)
        meta_dir = os.path.join(out_dir, 'meta')
        channel_dir = os.path.join(out_dir, 'channels')
        os.makedirs(meta_dir, exist_ok=True)
        os.makedirs(channel_dir, exist_ok=True)

        # Keep PINGverter's metadata side effect inside the project folder.
        self.metaDir = meta_dir
        self._ensure_garmin_metadata()

        pings_csv = os.path.join(out_dir, 'pings.csv')
        self.header_dat.to_csv(pings_csv, index=False)

        samples_path = os.path.join(out_dir, 'samples.u16le')
        frames_path = os.path.join(out_dir, 'frames.jsonl')
        frame_count = self.write_sonar_data_player_frames(samples_path, frames_path)

        waterfall_paths = {}
        if include_pngs:
            waterfall_paths = self.write_channel_waterfall_pngs(channel_dir, prefix=prefix)

        channels = []
        df = self.header_dat
        channel_ids = sorted(int(c) for c in df['channel_id'].dropna().unique())
        for channel_id in channel_ids:
            group = df[df['channel_id'] == channel_id]
            channel_info = self._channel_info_for_id(channel_id)
            channel_desc = self.describe_channel(channel_id, group, channel_info)
            sample_col = 'ping_cnt' if 'ping_cnt' in group.columns else 'sample_cnt'
            max_samples = int(group[sample_col].max()) if len(group) and sample_col in group.columns else 0

            channel = {
                'channelId': channel_id,
                'label': channel_desc['label'],
                'mode': channel_desc['mode'],
                'orientation': channel_desc['orientation'],
                'beam': channel_desc['beam'],
                'startFrequencyHz': channel_desc['startFrequencyHz'],
                'endFrequencyHz': channel_desc['endFrequencyHz'],
                'rows': int(len(group)),
                'maxSamples': max_samples,
                'timeStart': self._none_if_nan(group['time_s'].min()) if 'time_s' in group else None,
                'timeEnd': self._none_if_nan(group['time_s'].max()) if 'time_s' in group else None,
            }

            if channel_id in waterfall_paths:
                channel['waterfall'] = self._relpath(waterfall_paths[channel_id], out_dir)

            channels.append(channel)

        manifest = {
            'formatVersion': 2,
            'source': os.path.abspath(self.sonFile),
            'telemetry': self._relpath(pings_csv, out_dir),
            'frames': self._relpath(frames_path, out_dir),
            'samples': {
                'path': self._relpath(samples_path, out_dir),
                'encoding': 'uint16-le',
            },
            'frameCount': frame_count,
            'channels': channels,
        }

        manifest_path = os.path.join(out_dir, 'manifest.json')
        with open(manifest_path, 'w', encoding='utf-8') as file:
            json.dump(manifest, file, indent=2)

        return manifest_path

    # ======================================================================
    def write_sonar_data_player_frames(self, samples_path: str, frames_path: str,
                                       df: pd.DataFrame=None):
        """Write synchronized frame metadata and raw uint16 sonar samples."""
        if df is None:
            df = self.header_dat

        sample_col = 'ping_cnt' if 'ping_cnt' in df.columns else 'sample_cnt'
        frame_count = 0
        offset = 0

        with open(self.sonFile, 'rb') as rsd, open(samples_path, 'wb') as samples, open(frames_path, 'w', encoding='utf-8') as frames:
            for sequence_count, group in df.groupby('sequence_cnt', sort=True):
                channels = []

                for _, row in group.sort_values('channel_id').iterrows():
                    try:
                        sample_count = int(row[sample_col])
                        byte_count = sample_count * 2
                        data_size = int(row['data_size'])
                        ping_header_len = int(row.get('ping_header_len', self.pingHeaderLen))
                        son_offset = int(row['son_offset'])
                        record_index = int(row['index'])
                    except (KeyError, TypeError, ValueError):
                        continue

                    if sample_count <= 0:
                        continue
                    if son_offset < ping_header_len:
                        continue
                    if son_offset + byte_count > ping_header_len + data_size:
                        continue

                    rsd.seek(record_index + son_offset)
                    raw = rsd.read(byte_count)
                    if len(raw) != byte_count:
                        continue

                    samples.write(raw)
                    channels.append({
                        'channelId': int(row['channel_id']),
                        'sampleOffset': offset,
                        'sampleCount': sample_count,
                        'byteLength': byte_count,
                        'minRangeMeters': self._none_if_nan(row.get('min_range')),
                        'maxRangeMeters': self._none_if_nan(row.get('max_range')),
                        'bottomDepthMeters': self._none_if_nan(row.get('inst_dep_m')),
                    })
                    offset += byte_count

                if not channels:
                    continue

                frame = {
                    'frameIndex': frame_count,
                    'sequenceCount': int(sequence_count),
                    'timeSeconds': self._none_if_nan(group['time_s'].mean()) if 'time_s' in group else None,
                    'lat': self._none_if_nan(group['lat'].mean()) if 'lat' in group else None,
                    'lon': self._none_if_nan(group['lon'].mean()) if 'lon' in group else None,
                    'speedMetersPerSecond': self._none_if_nan(group['speed_ms'].mean()) if 'speed_ms' in group else None,
                    'trackDistanceMeters': self._none_if_nan(group['trk_dist'].mean()) if 'trk_dist' in group else None,
                    'headingDegrees': self._none_if_nan(group['instr_heading'].mean()) if 'instr_heading' in group else None,
                    'temperatureCelsius': self._none_if_nan(group['tempC'].mean()) if 'tempC' in group else None,
                    'channels': channels,
                }
                frames.write(json.dumps(frame, separators=(',', ':')) + '\n')
                frame_count += 1

        return frame_count

    # ======================================================================
    def describe_channel(self, channel_id: int, group: pd.DataFrame=None,
                         channel_info: dict=None):
        """Return display metadata for a Garmin channel."""
        if group is None:
            group = self.header_dat[self.header_dat['channel_id'] == channel_id]
        if channel_info is None:
            channel_info = self._channel_info_for_id(channel_id)

        beam = None
        if group is not None and 'beam' in group.columns and len(group['beam'].dropna()) > 0:
            beam = int(group['beam'].dropna().mode().iloc[0])

        start_hz = int(channel_info.get('start_freq_hz', 0) or 0)
        end_hz = int(channel_info.get('end_freq_hz', 0) or 0)

        orientation = None
        mode = 'Unknown'

        if beam == 1:
            mode = 'Traditional CHIRP'
        elif beam == 4:
            mode = 'Down Imaging'
        elif beam in (2, 3):
            mode = 'SideVu'
            port_star = group['port_star_id'].median() if group is not None and 'port_star_id' in group else None
            if port_star is not None and port_star < 0:
                orientation = 'Port'
            elif port_star is not None and port_star > 0:
                orientation = 'Starboard'
            elif beam == 2:
                orientation = 'Port'
            elif beam == 3:
                orientation = 'Starboard'

        freq = self._format_frequency_range(start_hz, end_hz)
        label_parts = [mode]
        if orientation:
            label_parts.append(orientation)
        if freq:
            label_parts.append(freq)

        return {
            'label': ' '.join(label_parts) if label_parts else 'Channel {}'.format(channel_id),
            'mode': mode,
            'orientation': orientation,
            'beam': beam,
            'startFrequencyHz': start_hz or None,
            'endFrequencyHz': end_hz or None,
        }

    # ======================================================================
    def _ensure_garmin_metadata(self):
        if not hasattr(self, 'file_len'):
            self._getFileLen()
        if not hasattr(self, 'file_header'):
            self._parseFileHeader()
        if not hasattr(self, 'header_dat'):
            if not hasattr(self, 'metaDir'):
                self.metaDir = os.path.dirname(os.path.abspath(self.sonFile))
            self._parsePingHeader()
            self._recalcRecordNum()

    # ======================================================================
    def _channel_info_for_id(self, channel_id: int):
        for item in getattr(self, 'channel_info', []):
            try:
                if int(item.get('channel_id', -1)) == int(channel_id):
                    return item
            except (TypeError, ValueError):
                continue
        return {}

    # ======================================================================
    def _format_frequency_range(self, start_hz: int, end_hz: int):
        if start_hz <= 0 or end_hz <= 0:
            return None

        return '{:.0f}-{:.0f} kHz'.format(start_hz / 1000, end_hz / 1000)

    # ======================================================================
    def _none_if_nan(self, value):
        try:
            f = float(value)
        except (TypeError, ValueError):
            return None
        return None if math.isnan(f) else f

    # ======================================================================
    def _relpath(self, path: str, root: str):
        return os.path.relpath(os.path.abspath(path), os.path.abspath(root)).replace(os.sep, '/')

    # ======================================================================
    def _scale_samples_for_waterfall(self, samples: np.ndarray):
        """Compress raw uint16 sonar intensities to an 8-bit palette index."""
        arr = samples.astype(np.float32)
        positive = arr[arr > 0]
        if positive.size == 0:
            return np.zeros(arr.shape, dtype=np.uint8)

        arr = np.log1p(arr)
        positive = arr[samples > 0]
        lo, hi = np.percentile(positive, [1, 99.5])
        if hi <= lo:
            hi = float(positive.max())
            lo = float(positive.min())
        if hi <= lo:
            return np.zeros(arr.shape, dtype=np.uint8)

        scaled = (arr - lo) * (255.0 / (hi - lo))
        return np.clip(scaled, 0, 255).astype(np.uint8)

    # ======================================================================
    def _garmin_waterfall_palette(self):
        """Approximate Garmin sonar colors as a 256-entry PIL palette."""
        stops = [
            (0, (0, 0, 0)),
            (24, (0, 20, 72)),
            (64, (0, 92, 160)),
            (104, (0, 178, 196)),
            (144, (30, 185, 70)),
            (184, (230, 205, 45)),
            (222, (230, 78, 30)),
            (255, (255, 245, 210)),
        ]
        palette = []
        for idx in range(256):
            for stop_idx in range(len(stops) - 1):
                a_idx, a_col = stops[stop_idx]
                b_idx, b_col = stops[stop_idx + 1]
                if a_idx <= idx <= b_idx:
                    t = 0 if b_idx == a_idx else (idx - a_idx) / (b_idx - a_idx)
                    col = tuple(int(round(a_col[c] + (b_col[c] - a_col[c]) * t)) for c in range(3))
                    palette.extend(col)
                    break
        return palette
    
    
    ### Ping Header Conversions ###
    # ======================================================================
    def _calcSpeedTrkDist(self, df: pd.DataFrame, jump_thresh: float=1.0):

        if len(df) <= 1:
            df['dist'] = 0.0
            df['speed_ms'] = 0.0
            df['trk_dist'] = 0.0
            return df

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
    def _recomputeTrackSpeedFromWgs(self, df: pd.DataFrame):
        """Recompute speed and cumulative track distance from WGS84 positions."""
        required = {'sequence_cnt', 'time_s', 'lat', 'lon'}
        if not required.issubset(df.columns) or len(df) == 0:
            return df

        frame_track = (
            df.groupby('sequence_cnt', sort=True)
            .agg(time_s=('time_s', 'mean'), lat=('lat', 'mean'), lon=('lon', 'mean'))
            .reset_index()
        )
        if len(frame_track) == 0:
            return df

        distances = [0.0]
        speeds = [np.nan]
        cumulative = [0.0]

        for idx in range(1, len(frame_track)):
            prev = frame_track.iloc[idx - 1]
            cur = frame_track.iloc[idx]
            dist = self._haversineMeters(prev['lat'], prev['lon'], cur['lat'], cur['lon'])
            dt = float(cur['time_s'] - prev['time_s'])

            distances.append(dist)
            speeds.append(dist / dt if dt > 0 else np.nan)
            cumulative.append(cumulative[-1] + dist)

        frame_track['dist'] = distances
        frame_track['speed_ms'] = speeds
        frame_track['trk_dist'] = cumulative
        frame_track['speed_ms'] = frame_track['speed_ms'].interpolate().bfill().ffill().round(2)

        replacements = frame_track.set_index('sequence_cnt')[['dist', 'speed_ms', 'trk_dist']]
        for column in replacements.columns:
            df[column] = df['sequence_cnt'].map(replacements[column])

        return df

    # ======================================================================
    def _haversineMeters(self, lat1: float, lon1: float, lat2: float, lon2: float):
        radius_m = 6371000.0
        p1 = math.radians(float(lat1))
        p2 = math.radians(float(lat2))
        dp = math.radians(float(lat2) - float(lat1))
        dl = math.radians(float(lon2) - float(lon1))
        a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
        return 2 * radius_m * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
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



        self.humDat['wgs'] = "EPSG:4326"

        if pyproj is not None and len(df) > 0:
            # Determine epsg code
            self.humDat['epsg'] = "EPSG:"+str(int(float(self._convert_wgs_to_utm(df['lon'].iloc[0], df['lat'].iloc[0]))))

            # Configure re-projection function
            self.trans = pyproj.Proj(self.humDat['epsg'])

            # Reproject lat/lon to UTM zone
            e, n = self.trans(df['lon'], df['lat'])
            df['e'] = e
            df['n'] = n
        else:
            self.humDat['epsg'] = self.humDat['wgs']
            self.trans = None
            df['e'] = df['lon']
            df['n'] = df['lat']


        #########################
        # Calculate COG (heading)
        ## Garmin does not appear to store heading....
        heading = self._getCOG(df)
        # self._getBearing() returns n-1 values because last ping can't
        ## have a COG value.  We will duplicate the last COG value and use it for
        ## the last ping.
        last = heading[-1] if len(heading) > 0 else 0
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
