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

# RSD structur
rsdStruct = np.dtype([
    ("test", "<u4"),
])

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

        self.son_struct = rsdStruct

        self.humDat = {} # Store general sonar recording metadata

        self.son8bit = True

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
    
    
    # ======================================================================
    def _parseFileHeader(self):
        '''
        '''

        self.headBytes = 20480 # Hopefully a fixed value for all RSD files
        chanInfoLen = 1069 # It is not clear if there is helpful info in channel information...

        # Get the file header structure
        headStruct, firstHeadBytes = self._getFileHeaderStruct()

        print('\n\n\nheadStruct:')
        for v in headStruct:
            print(v)
        
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

        for k,v in out_dict.items():
            print(k, v)

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
            #     field_cnt_0 = self._fread_dat(file, 1, 'B')[0] # Decode the spacer byte

            #     chanidx = 0
            #     fidx_0 = 0

            #     # Iterate each channel
            #     while chanidx < chan_cnt:
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
            print('lastPos:', lastPos)        

        return headStruct, lastPos-1
    
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

        # Test file to see outputs
        out_test = os.path.join(self.metaDir, 'All-Garmin-Sonar-MetaData.csv')
        print(out_test)
        df.to_csv(out_test, index=False)



        return
    
    # ======================================================================
    def _getPingHeader(self, file, i: int):

        print('\n\n\n', i)

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
            10:[('SP0a', '<u1'), ('unknown_sp0a', '<u2')],
            11:[('SP0b', '<u1'), ('bottom_depth_unknown', '<u1'), ('bottom_depth', '<u2')], #0b bottom depth
            13:[('SP0d', '<u1'), ('unknown_sp0d', '<u4'), ('unknown_sp0d_1', '<u1')],
            18:[('SP12', '<u1'), ('unknown_sp12', '<u2')],
            19:[('SP13', '<u1'), ('drawn_bottom_depth_unknown', '<u1'), ('drawn_bottom_depth', '<u2')], #13 drawn bottom depth
            21:[('SP15', '<u1'), ('unknown_sp15', '<u4'), ('unknown_sp15_1', '<u1')],
            25:[('SP19', '<u1'), ('first_sample_depth', '<u1')], #19 first sample depth
            35:[('SP23', '<u1'), ('last_sample_depth_unknown', '<u1'), ('last_sample_depth', '<u2')], #23 last sample depth
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
                ('su2_f0', '<u1'), ('su2_f0_unknown', '<f4'),
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
        record_body_header_len = 1

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

            


        # Next ping header is from current position + ping_cnt
        # next_ping = file.tell() + out_dict['packet_size']
        next_ping = i + pingHeaderLen + out_dict['data_size'] + 12 #12 for magic number & crc

        out_dict['index'] = i

        out_dict['son_offset'] = out_dict['data_size'] - out_dict['sample_cnt']
 
        return out_dict, next_ping
    
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

        pingBodyHeaderToCheck = {
            -1:('record_body_fcnt', '<u1'),
            1:[('SP1', '<u1'), ('channel_id_1', '<u1')], #01 channel_id
            11:[('SP0b', '<u1'), ('bottom_depth_unknown', '<u1'), ('bottom_depth', '<u2')], #0b bottom depth
            13:[('SP0d', '<u1'), ('unknown_sp0d', '<u4'), ('unknown_sp0d_1', '<u1')],
            18:[('SP12', '<u1'), ('unknown_sp12', '<u2')],
            19:[('SP13', '<u1'), ('drawn_bottom_depth_unknown', '<u1'), ('drawn_bottom_depth', '<u2')], #13 drawn bottom depth
            21:[('SP15', '<u1'), ('unknown_sp15', '<u4'), ('unknown_sp15_1', '<u1')],
            25:[('SP19', '<u1'), ('first_sample_depth', '<u1')], #19 first sample depth
            35:[('SP23', '<u1'), ('last_sample_depth_unknown', '<u1'), ('last_sample_depth', '<u2')], #23 last sample depth
            41:[('SP29', '<u1'), ('gain', '<u1')], #29 gain
            49:[('SP31', '<u1'), ('sample_status', '<u1')], #31 sample status
            60:[('SP3c', '<u1'), ('sample_cnt', '<u4')], #3c sample count
            65:[('SP41', '<u1'), ('shade_avail', '<u1')], #41 shade available
            76:[('SP4c', '<u1'), ('scposn_lat', '<u4')], #4c latitude
            84:[('SP54', '<u1'), ('scposn_lon', '<u4')], #54 longitude
            92:[('SP5c', '<u1'), ('water_temp', '<f4')], #5c temperature
            97:[('SP61', '<u1'), ('beam', '<u1')], #61 beam
        }

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

                # # Add pingheader
                # for i in pingHeader:
                #     headStruct.append(i)

                # Get field count
                field_cnt = self._fread_dat(file, 1, 'B')[0]
                headStruct.append(pingBodyHeaderToCheck[-1])

                if field_cnt > 13: # Only 13 known fields. Some beams have up to 15
                    field_cnt = 13

                fidx = 0
                record_body_header_len = 1

                while fidx < field_cnt:

                    byte = self._fread_dat(file, 1, 'B')[0]

                    if byte in pingBodyHeaderToCheck:
                        elen = 0
                        for v in pingBodyHeaderToCheck[byte]:
                            headStruct.append(v)

                            # Get length of element
                            elen += (np.dtype(v[-1]).itemsize)

                        # Move forward elen amount
                        cpos = file.tell()
                        npos = cpos + elen - 1

                        record_body_header_len += elen

                        file.seek(npos)

                        fidx += 1

                    else:
                        print('{} not in sonar body. Terminating.'.format(byte))
                        print('Offset: {}'.format(file.tell()))
                        sys.exit()

                    foundChanInfo = True

        # self.son_header_struct = pingHeader
        # self.son_struct = headStruct

        return headStruct, pingHeader, record_body_header_len
            







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