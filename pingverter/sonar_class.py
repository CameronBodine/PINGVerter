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

import os
import pandas as pd
import numpy as np

class sonar(object):

    def __init__(self):

        print('sonar class poop')

        pass

    def _parsePingHeader(self, in_file: str, out_file: str=None):
        '''
        '''

        # Get file length
        file_len = self.file_len = os.path.getsize(in_file)

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

            print(i)

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
        if out_file:
            header_dat_all.to_csv(out_file, index=False)
        else:
            self.header_dat = header_dat_all

        return
    
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
        next_ping = int(file.tell() + header[0][-2])

        return out_dict, next_ping
