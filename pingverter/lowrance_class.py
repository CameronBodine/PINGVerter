 
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
import os
import numpy as np
import pandas as pd

try:
    import pyproj
except ImportError:
    pyproj = None

'''
Based on Sonarlight by Kenneth Thorø Martinsen
The package is inspired by and builds upon other 
tools and descriptions for processing Lowrance 
sonar data, e.g. SL3Reader which includes a usefull 
paper, python-sllib, sonaR, Navico_SLG_Format notes, 
older blog post.
'''

#dtype for '.sl2' files (144 bytes)
sl2Struct = np.dtype([
    ("frame_offset", "<u4"),
    ("prev_primary_offset", "<u4"),
    ("prev_secondary_offset", "<u4"),
    ("prev_downscan_offset", "<u4"),
    ("prev_left_sidescan_offset", "<u4"),
    ("prev_right_sidescan_offset", "<u4"),
    ("prev_sidescan_offset", "<u4"),
    ("frame_size", "<u2"),
    ("prev_frame_size", "<u2"),
    ("survey_type", "<u2"),
    ("packet_size", "<u2"),
    ("id", "<u4"),
    ("min_range", "<f4"),
    ("max_range", "<f4"),
    ("unknown48", "<u2"),
    ("unknown50", "<B"),
    ("unknown51", "<B"),
    ("unknown52", "<B"),
    ("frequency_type", "<B"), # Frequency
    ("unknown54", "<u2"),
    ("unknown56", "<u2"),
    ("unknown58", "<u2"),
    ("hardware_time", "<u4"),
    ("depth_ft", "<f4"),
    ("keel_depth_ft", "<f4"),
    ('unknown72', '<B'),
    ('unknown73', '<B'),
    ('unknown74', '<u2'),
    ('unknown76', '<B'),
    ('unknown77', '<B'),
    ('unknown78', '<u2'),
    ('unknown80', '<f4'),    
    ('unknown84', '<f4'),
    ('unknown88', '<f4'),
    ('unknown92', '<f4'),
    ('unknown96', '<B'),
    ('unknown97', '<B'),
    ('unknown98', '<B'),
    ('unknown99', '<B'),
    ("gps_speed", "<f4"), #[knots]
    ("water_temperature", "<f4"), #[C]
    ("utm_e", "<i4"), #Easting in mercator [meters]
    ("utm_n", "<i4"), #Northing in mercator [meters]
    ("water_speed", "<f4"), #Water speed through paddlewheel or GPS if not present [knots]
    ("track_cog", "<f4"), # Track (COG) [radians]
    ("altitude", "<f4"), # Above sea level [feet]
    ("heading", "<f4"), #[radians]
    ("flags", "<u2"),
    ('unknown134', '<u2'),
    ('unknown136', '<B'),
    ('unknown137', '<B'),
    ('unknown138', '<B'),
    ('unknown139', '<B'),
    ("time_s", "<u4") # Time since beginning of log [ms]
])

#dtype for '.sl3' files (168 bytes)
# sl3Struct = np.dtype([
#     ("frame_offset", "<u4"),
#     ("frame_version", '<u4'),
#     ("frame_size", "<u2"),
#     ("prev_frame_size", "<u2"),
#     ("survey_type", "<u2"),
#     ("unknown14", "<i2"),
#     ("id", "<u4"),
#     ("min_range", "<f4"),
#     ("max_range", "<f4"),
#     ("unknown28", "<f4"),
#     ("unknown32", "<f4"),
#     ("unknown36", "<f4"),
#     ("hardware_time", "<u4"),
#     ("packet_size", "<u2"),
#     ("unknown46", "<u2"),
#     ("depth_ft", "<f4"),
#     ("frequency_type", "<B"),
#     ("unknown53", "<B"),
#     ("unknown54", "<B"),
#     ("unknown55", "<B"),
#     ("unknown56", "<B"),
#     ("unknown57", "<B"),
#     ("unknown58", "<B"),
#     ("unknown59", "<B"),
#     ("unknown60", "<B"),
#     ("unknown61", "<B"),
#     ("unknown62", "<B"),
#     ("unknown63", "<B"),
#     ("unknown64", "<f4"),
#     ("unknown68", "<f4"),
#     ("unknown72", "<f4"),
#     ("unknown76", "<f4"),
#     ("unknown80", "<B"),
#     ("unknown81", "<B"),
#     ("unknown82", "<B"),
#     ("unknown83", "<B"),
#     ("gps_speed", "<f4"), #[knots]
#     ("water_temperature", "<f4"), #[C]
#     ("utm_e", "<i4"), #Easting in mercator [meters]
#     ("utm_n", "<i4"), #Northing in mercator [meters]
#     ("water_speed", "<u4"), #Water speed through paddlewheel or GPS if not present [knots]
#     ("track_cog", "<f4"), # Track (COG) [radians]
#     ("altitude", "<f4"), # Above sea level [feet]
#     ("heading", "<f4"), #[radians]
#     ("unknown116", 'i4'),
#     ("unknown120", "<B"),
#     ("unknown121", "<B"),
#     ("unknown122", "<B"),
#     ("unknown123", "<B"),
#     ("time_s", "<u4"), # Time since beginning of log [ms]
#     ("prev_primary_offset", "<u4"),
#     ("prev_secondary_offset", "<u4"),
#     ("prev_downscan_offset", "<u4"),
#     ("prev_left_sidescan_offset", "<u4"),
#     ("prev_right_sidescan_offset", "<u4"),
#     ("prev_sidescan_offset", "<u4"),
#     ("unknown152", "<u4"),
#     ("unknown156", "<u4"),
#     ("unknown160", "<u4"),
#     ("prev_3d_offseft", "<u4")
# ])

sl3Struct = np.dtype([
    ("frame_offset", "<u4"),
    ("frame_version", "<u4"),
    ("frame_size", "<u2"),
    ("prev_frame_size", "<u2"),
    ("survey_type", "<u2"),
    ("unknown14", "<i2"),
    ("id", "<u4"),
    ("min_range", "<f4"),
    ("max_range", "<f4"),
    ("unknown28", "<f4"),
    ("unknown32", "<f4"),
    ("unknown36", "<f4"),
    ("hardware_time", "<u4"),
    ("packet_size", "<u4"),
    ("depth_ft", "<f4"),
    ("frequency_type", "<u2"),
    ("unknown54", "<f4"),
    ("unknown58", "<f4"),
    ("unknown62", "<i2"),
    ("unknown64", "<f4"),
    ("unknown68", "<f4"),
    ("unknown72", "<f4"),
    ("unknown76", "<f4"),
    ("unknown80", "<f4"),
    ("gps_speed", "<f4"), #[knots]
    ("water_temperature", "<f4"), #[C]
    ("utm_e", "<i4"), #Easting in mercator [meters]
    ("utm_n", "<i4"), #Northing in mercator [meters]
    ("water_speed", "<f4"), #Water speed through paddlewheel or GPS if not present [knots]
    ("track_cog", "<f4"), # Track (COG) [radians]
    ("altitude", "<f4"), # Above sea level [feet]
    ("heading", "<f4"), #[radians]
    ("flags", "<u2"),
    ("unknown118", "<u2"),
    ("unknown120", "<u4"),
    ("time_s", "<u4"), # Time since beginning of log [ms]
    ("prev_primary_offset", "<u4"),
    ("prev_secondary_offset", "<u4"),
    ("prev_downscan_offset", "<u4"),
    ("prev_left_sidescan_offset", "<u4"),
    ("prev_right_sidescan_offset", "<u4"),
    ("prev_sidescan_offset", "<u4"),
    ("unknown152", "<u4"),
    ("unknown156", "<u4"),
    ("unknown160", "<u4"),
    ("prev_3d_offseft", "<u4")
])

# Map Lowrance ping attribute names to PING-Mapper (PM)
lowCols2PM = {
    'track_cog': 'instr_heading',
    'heading': 'heading_magnetic',
    'gps_speed': 'speed_ms',
    'depth_ft': 'inst_dep_m',
    'packet_size': 'ping_cnt',
    'frame_offset': 'index',
    'keel_depth_ft': 'keel_depth_m'

}

class low(object):

    def __init__(self, inFile: str, nchunk: int=0, exportUnknown: bool=False):

        '''
        '''

        self.humFile = None
        self.sonFile = inFile
        self.nchunk = nchunk
        self.exportUnknown = exportUnknown

        self.file_header_size = 8

        self.extension = os.path.basename(inFile).split('.')[-1]

        self.frame_header_size = 168 if "sl3" in self.extension else 144
        self.son_struct = sl3Struct if "sl3" in self.extension else sl2Struct

        self.lowCols2PM = lowCols2PM

        self.humDat = {} # Store general sonar recording metadata

        self.survey_dict = {0: 'primary', 1: 'secondary', 2: 'downscan',
                            3: 'left_sidescan', 4: 'right_sidescan', 5: 'sidescan',
                            9: '3D', 10: 'debug_digital', 11: 'debug_noise'}
        
        self.frequency_dict = {0: "200kHz", 1: "50kHz", 2: "83kHz",
                               3: "455kHz", 4: "800kHz", 5: "38kHz", 
                               6: "28kHz", 7: "130kHz_210kHz", 8: "90kHz_150kHz", 
                               9: "40kHz_60kHz", 10: "25kHz_45kHz"}
        
        self.son8bit = True
        
        return
    
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
    
    def _getFileLen(self):
        self.file_len = os.path.getsize(self.sonFile)

        return
    
    def _parseFileHeader(self):

        self.header = {0: [0, 0, 2, 'format', '<u2'],
                       2: [2, 0, 2, 'version', '<u2'],
                       4: [4, 0, 2, 'bytes_per_sounding', '<u2'],
                       6: [6, 0, 1, 'debug', 'B'],
                       7: [7, 0, 1, 'byte', 'B']}
        

        # Open sonar log
        f = open(self.sonFile, 'rb')

        # Iterate known file header items
        header = dict()
        for k, v in self.header.items():
            offset = v[0]
            length = v[2]
            name = v[3]
            type = v[4]
            f.seek(offset)

            v = self._fread_dat(f, length, type)
            header[name] = v.item()

        # Set class attribtutes
        self.file_header = header

        return
    
    def _parsePingHeader(self):
        '''
        '''

        # Get the file length
        file_len = self.file_len

        # Initialize offset after file header
        i = self.file_header_size

        # Open the file
        file = open(self.sonFile, 'rb')

        # Store contents in list
        header_dat_all = []

        # # counter for testing
        # test_cnt = 0

        # Decode ping header
        while i < file_len:

            # Get header data at offset i
            header_dat, cpos = self._getPingHeader(file, i)

            # Store the data
            header_dat_all.append(header_dat)

            # Update counter
            i = cpos


            # test_cnt += 1

            # if test_cnt == 50:
            #     break

        # Convert to dataframe
        df = pd.DataFrame.from_dict(header_dat_all)

        # Do unit conversions to PING-Mapper units
        df = self._doUnitConversion(df)

        # Do column conversions to PING-Mapper column names
        df.rename(columns=self.lowCols2PM, inplace=True)

        # Calculate along-track distance from 'time's and 'speed_ms'. Approximate distance estimate
        df = self._calcTrkDistTS(df)

        # Determine beams present
        df = self._convertBeam(df)

        # Convert Lowrance frequency
        df = self._convertLowFrequency(df)

        # Store sonar offset 
        df['son_offset'] = self.frame_header_size

        # Test file to see outputs
        out_test = os.path.join(self.metaDir, 'All-Lowrance-Sonar-MetaData.csv')
        df.to_csv(out_test, index=False)

        self.header_dat = df

        return
    
    def _getPingHeader(self, file, i: int):

        # Get necessary attributes
        head_struct = self.son_struct
        length = self.frame_header_size

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
        # next_ping = file.tell() + out_dict['packet_size']
        next_ping = i + out_dict['frame_size']

        return out_dict, next_ping
    
    def _doUnitConversion(self, df: pd.DataFrame):

        # Convert feet to meters
        if self.extension == "sl2":
            df[["depth_ft", "keel_depth_ft", "min_range", "max_range", "altitude"]] /= 3.2808399
        else:
            df[["depth_ft", "min_range", "max_range", "altitude"]] /= 3.2808399

        # convert time [ms] to s
        df['time_s'] /= 1000

        # Convert speed [knots] to m/s
        df['gps_speed'] *= 0.514444

        # Calculate caltime
        hardware_time_start = df["hardware_time"][0]
        df['caltime'] = pd.to_datetime(hardware_time_start + df['time_s'], unit='s')

        df['date'] = df['caltime'].dt.date
        df['time'] = df['caltime'].dt.time
        df = df.drop('caltime', axis=1)

        # Calculate latitude and longitude
        df['lat'] = (((2*np.arctan(np.exp(df['utm_n']/6356752.3142)))-(np.pi/2))*(180/np.pi))
        df['lon'] = (df['utm_e']/6356752.3142*(180/np.pi))

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

        # Convert radians to degrees
        df['track_cog'] = np.rad2deg(df['track_cog'])
        df['heading'] = np.rad2deg(df['heading'])

        # Store survey temperature
        df['tempC'] = self.tempC*10

        # Add transect number (for aoi processing)
        df['transect'] = 0

        # Calculate pixel size [m]  *** ....MAYBE.... ***
        df['pixM'] = (df['max_range'] - df['min_range']) / df['packet_size']

        # Calculate frequency and type of beam
        df["survey"] = [self.survey_dict.get(i, "unknown") for i in df["survey_type"]]
        df["frequency"] = [self.frequency_dict.get(i, "unknown") for i in df["frequency_type"]]
        
        

        return df
    
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
    
    def _convertBeam(self, df: pd.DataFrame):
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

        df['beam'] = [beam_xwalk.get(i, "unknown") for i in df['survey_type']]

        return df

    def _convertLowFrequency(self, df: pd.DataFrame):

        '''
        Crosswalk Lowrance frequency to PING-Mapper.
        PM has slots for frequency, min-frequency, max-frequency

        {lowrance-frequency: [PM Frequecy, min, max]}
        '''
        
        frequency_xwalk = {'200kHz': [200, 200, 200], '50kHz': [50, 50, 50],
                           '83kHz': [83, 83, 83], '455kHz': [455, 455, 455],
                           '800kHz': [800, 800, 800], '38kHz': [38, 38, 38],
                           '28kHz': [28, 28, 28], '130kHz_210kHz': [170, 130, 210],
                           '90kHz_150kHz': [120, 90, 150], '40kHz_60kHz': [50, 40, 60],
                           '25kHz_45kHz': [35, 25, 45]}
        
        frequency_min = {200: 200, 50: 50, 83: 83, 455: 455, 800: 800, 38: 38,
                         28: 28, 170: 130, 120:90, 50: 40, 35: 25}
        
        # df['f'] = [frequency_xwalk[i][0] for i in df['frequency']]
        df["f"] = [frequency_xwalk.get(i, -1) for i in df["frequency_type"]]
        
        df['f_min'] = [frequency_xwalk.get(i, -1) for i in df["frequency_type"]]
        df['f_max'] = [frequency_xwalk.get(i, -1) for i in df["frequency_type"]]

        return df

    def _removeUnknownBeams(self):

        df = self.header_dat

        # Drop unknown
        df = df[df['beam'] != 'unknown']

        self.header_dat = df
        return
    
    def _removeDownBeams(self):
        '''
        PING-Mapper expects low-frequency (83kHz) stored as beam 0
        and high-frequency(200kHz) stored as beam 2
        '''

        df = self.header_dat
        dfDown = df[df['beam'] < 2]

        dfRest = df[df['beam'] > 1]

        for beam, group in dfDown.groupby('beam'):
            if beam == 0:
                f = group['f'].iloc[0]
                if -1 < f < 200:
                    dfRest = pd.concat([dfRest, group])

            elif beam == 1:
                f = group['f'].iloc[0]
                if f >= 200:
                    dfRest = pd.concat([dfRest, group])

        self.header_dat = dfRest

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

        # set min_range to 0
        port['min_range'] = 0
        star['min_range'] = 0

        # Concatenate df's
        dfAll = pd.concat([dfAll, port, star], ignore_index=True)

        dfAll.sort_values(by=['time_s', 'beam'], inplace=True)

        self.header_dat = dfAll

        return

    def extract_raw_sample_arrays(self, df: pd.DataFrame=None, expand_to_uint16: bool=False):
        """Return raw Lowrance samples per ping grouped by beam.

        Lowrance SL2/SL3 sample payloads are one byte per sounding. When
        expand_to_uint16 is True, samples are scaled into the full uint16 range
        for consumers that share Garmin's uint16 project format.
        """
        if df is None:
            df = self.header_dat

        samples_by_channel = {}

        with open(self.sonFile, 'rb') as file:
            for _, row in df.iterrows():
                try:
                    channel_id = int(row['beam'])
                    record_index = int(row['index'])
                    sample_count = int(row['ping_cnt'])
                    son_offset = int(row['son_offset'])
                    frame_size = int(row['frame_size'])
                except (KeyError, TypeError, ValueError):
                    continue

                if sample_count <= 0 or son_offset < self.frame_header_size:
                    continue
                if son_offset + sample_count > frame_size:
                    continue

                file.seek(record_index + son_offset)
                raw = file.read(sample_count)
                if len(raw) != sample_count:
                    continue

                arr = np.frombuffer(raw, dtype=np.uint8)
                if expand_to_uint16:
                    arr = arr.astype('<u2') * 257
                else:
                    arr = arr.copy()

                samples_by_channel.setdefault(channel_id, []).append(arr)

        return samples_by_channel

    def write_channel_waterfall_pngs(self, out_dir: str, df: pd.DataFrame=None,
                                     prefix: str=None, width: int=None):
        """Write one Lowrance waterfall PNG per beam and return paths."""
        from PIL import Image

        os.makedirs(out_dir, exist_ok=True)
        samples_by_channel = self.extract_raw_sample_arrays(df)
        out_paths = {}

        for channel_id, pings in samples_by_channel.items():
            if not pings:
                continue

            max_len = width or max(len(p) for p in pings)
            img_arr = np.zeros((len(pings), max_len), dtype=np.uint8)

            for row_idx, ping in enumerate(pings):
                n = min(len(ping), max_len)
                img_arr[row_idx, :n] = ping[:n]

            image = Image.fromarray(img_arr, mode='P')
            image.putpalette(self._lowrance_waterfall_palette())

            stem = prefix or os.path.splitext(os.path.basename(self.sonFile))[0]
            safe_stem = ''.join(c if c.isalnum() or c in ('-', '_') else '_' for c in stem)
            out_path = os.path.join(out_dir, '{}_channel_{}.png'.format(safe_stem, channel_id))
            image.save(out_path)
            out_paths[channel_id] = out_path

        return out_paths

    def write_sonar_data_player_project(self, out_dir: str, include_pngs: bool=True,
                                        include_unknown: bool=False, prefix: str=None):
        """Write a SonarDataPlayer processed project for this Lowrance file."""
        os.makedirs(out_dir, exist_ok=True)
        meta_dir = os.path.join(out_dir, 'meta')
        channel_dir = os.path.join(out_dir, 'channels')
        os.makedirs(meta_dir, exist_ok=True)
        os.makedirs(channel_dir, exist_ok=True)

        self.metaDir = meta_dir
        self._ensure_lowrance_metadata(include_unknown=include_unknown)

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
        channel_ids = sorted(int(c) for c in df['beam'].dropna().unique())
        for channel_id in channel_ids:
            group = df[df['beam'] == channel_id]
            channel_desc = self.describe_channel(channel_id, group)
            max_samples = int(group['ping_cnt'].max()) if len(group) else 0
            channel = {
                'channelId': channel_id,
                'label': channel_desc['label'],
                'mode': channel_desc['mode'],
                'orientation': channel_desc['orientation'],
                'beam': channel_id,
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
                'sourceEncoding': 'uint8-expanded',
            },
            'frameCount': frame_count,
            'channels': channels,
        }

        manifest_path = os.path.join(out_dir, 'manifest.json')
        with open(manifest_path, 'w', encoding='utf-8') as file:
            json.dump(manifest, file, indent=2)

        return manifest_path

    def write_sonar_data_player_frames(self, samples_path: str, frames_path: str,
                                       df: pd.DataFrame=None):
        """Write synchronized frame metadata and uint16-expanded Lowrance samples."""
        if df is None:
            df = self.header_dat

        frame_count = 0
        offset = 0

        channel_groups = {
            int(channel_id): group.sort_values('time_s').reset_index(drop=True)
            for channel_id, group in df.groupby('beam', sort=True)
        }
        max_frames = max((len(group) for group in channel_groups.values()), default=0)

        with open(self.sonFile, 'rb') as source, open(samples_path, 'wb') as samples, open(frames_path, 'w', encoding='utf-8') as frames:
            for frame_idx in range(max_frames):
                channels = []
                rows = []

                for channel_id, group in channel_groups.items():
                    if frame_idx >= len(group):
                        continue

                    row = group.iloc[frame_idx]
                    try:
                        sample_count = int(row['ping_cnt'])
                        frame_size = int(row['frame_size'])
                        son_offset = int(row['son_offset'])
                        record_index = int(row['index'])
                    except (TypeError, ValueError):
                        continue

                    if sample_count <= 0:
                        continue
                    if son_offset < self.frame_header_size:
                        continue
                    if son_offset + sample_count > frame_size:
                        continue

                    source.seek(record_index + son_offset)
                    raw = source.read(sample_count)
                    if len(raw) != sample_count:
                        continue

                    expanded = (np.frombuffer(raw, dtype=np.uint8).astype('<u2') * 257).tobytes()
                    byte_count = len(expanded)
                    samples.write(expanded)

                    channels.append({
                        'channelId': channel_id,
                        'sampleOffset': offset,
                        'sampleCount': sample_count,
                        'byteLength': byte_count,
                        'minRangeMeters': self._none_if_nan(row.get('min_range')),
                        'maxRangeMeters': self._none_if_nan(row.get('max_range')),
                        'bottomDepthMeters': self._none_if_nan(row.get('inst_dep_m')),
                    })
                    offset += byte_count
                    rows.append(row)

                if not channels:
                    continue

                frame_df = pd.DataFrame(rows)
                frame = {
                    'frameIndex': frame_count,
                    'sequenceCount': frame_idx,
                    'timeSeconds': self._none_if_nan(frame_df['time_s'].mean()) if 'time_s' in frame_df else None,
                    'lat': self._none_if_nan(frame_df['lat'].mean()) if 'lat' in frame_df else None,
                    'lon': self._none_if_nan(frame_df['lon'].mean()) if 'lon' in frame_df else None,
                    'speedMetersPerSecond': self._none_if_nan(frame_df['speed_ms'].mean()) if 'speed_ms' in frame_df else None,
                    'trackDistanceMeters': self._none_if_nan(frame_df['trk_dist'].mean()) if 'trk_dist' in frame_df else None,
                    'headingDegrees': self._none_if_nan(frame_df['instr_heading'].mean()) if 'instr_heading' in frame_df else None,
                    'temperatureCelsius': self._none_if_nan(frame_df['tempC'].mean()) if 'tempC' in frame_df else None,
                    'channels': channels,
                }
                frames.write(json.dumps(frame, separators=(',', ':')) + '\n')
                frame_count += 1

        return frame_count

    def describe_channel(self, channel_id: int, group: pd.DataFrame=None):
        """Return display metadata for a Lowrance beam."""
        mode_by_beam = {
            0: ('Primary', None),
            1: ('Secondary', None),
            2: ('SideScan', 'Port'),
            3: ('SideScan', 'Starboard'),
            4: ('DownScan', None),
        }
        mode, orientation = mode_by_beam.get(channel_id, ('Unknown', None))

        frequency = None
        if group is not None and 'frequency' in group and len(group['frequency'].dropna()) > 0:
            frequency = str(group['frequency'].dropna().mode().iloc[0])

        start_hz, end_hz = self._frequency_range_hz(frequency)
        label_parts = [mode]
        if orientation:
            label_parts.append(orientation)
        if frequency and frequency != 'unknown':
            label_parts.append(frequency.replace('_', '-'))

        return {
            'label': ' '.join(label_parts),
            'mode': mode,
            'orientation': orientation,
            'startFrequencyHz': start_hz,
            'endFrequencyHz': end_hz,
        }

    def _ensure_lowrance_metadata(self, include_unknown: bool=False):
        if not hasattr(self, 'tempC'):
            self.tempC = 1.0
        if not hasattr(self, 'file_len'):
            self._getFileLen()
        if not hasattr(self, 'file_header'):
            self._parseFileHeader()
        if not hasattr(self, 'header_dat'):
            self._parsePingHeader()

        if not include_unknown:
            self._removeUnknownBeams()

        if 5 in set(self.header_dat['beam'].dropna()):
            self._splitLowSS()

        self._recalcRecordNum()

    def _frequency_range_hz(self, frequency: str):
        if not frequency or frequency == 'unknown':
            return None, None
        text = frequency.replace('kHz', '')
        if '_' in text:
            parts = text.split('_')
            try:
                return int(float(parts[0]) * 1000), int(float(parts[-1]) * 1000)
            except (TypeError, ValueError):
                return None, None
        try:
            value = int(float(text) * 1000)
            return value, value
        except (TypeError, ValueError):
            return None, None

    def _lowrance_waterfall_palette(self):
        """Approximate Lowrance/sonar display colors as a 256-entry PIL palette."""
        stops = [
            (0, (0, 0, 0)),
            (32, (0, 30, 90)),
            (72, (0, 105, 180)),
            (116, (0, 190, 210)),
            (156, (40, 190, 75)),
            (196, (235, 210, 45)),
            (226, (235, 88, 30)),
            (255, (255, 248, 218)),
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

    def _none_if_nan(self, value):
        try:
            f = float(value)
        except (TypeError, ValueError):
            return None
        return None if math.isnan(f) else f

    def _relpath(self, path: str, root: str):
        return os.path.relpath(os.path.abspath(path), os.path.abspath(root)).replace(os.sep, '/')
    
    def _recalcRecordNum(self):

        df = self.header_dat

        # Reset index and recalculate record num
        ## Record num is unique for each ping across all sonar beams
        df = df.reset_index(drop=True)
        df['record_num'] = df.index

        self.header_dat = df
        return
    
    def _splitBeamsToCSV(self):

        '''
        '''

        # Dictionary to store necessary attributes for PING-Mapper
        self.beamMeta = beamMeta = {}

        # Get df
        df = self.header_dat

        # Iterate each beam
        for beam, group in df.groupby('beam'):
            meta = {}

            # Set pixM based on side scan
            if beam == 2 or beam == 3:
                self.pixM = group['pixM'].iloc[0]
            

            # Determine beam name
            beam = 'B00'+str(beam)
            meta['beamName'] = self._getBeamName(beam)

            # Store sonFile
            meta['sonFile'] = self.sonFile

            # Drop columns
            group.drop(columns=['survey_type', 'frequency_type', 'survey', 'frequency'], inplace=True)

            # Add chunk_id
            group = self._getChunkID(group)

            # Save csv
            outCSV = '{}_{}_meta.csv'.format(beam, meta['beamName'])
            outCSV = os.path.join(self.metaDir, outCSV)
            group.to_csv(outCSV, index=False)

            meta['metaCSV'] = outCSV

            # Store the beams metadata
            beamMeta[beam] = meta


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

   

