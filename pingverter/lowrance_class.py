
import os
import numpy as np
import pandas as pd
import pyproj

'''
Based on Sonarlight by Kenneth Thorø Martinsen
The package is inspired by and builds upon other tools and descriptions for processing Lowrance sonar data, e.g. SL3Reader which includes a usefull paper, python-sllib, sonaR, Navico_SLG_Format notes, older blog post.
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

        # Determine epsg code
        self.humDat['epsg'] = "epsg:"+str(int(float(self._convert_wgs_to_utm(df['lon'][0], df['lat'][0]))))
        self.humDat['wgs'] = "epsg:4326"

        # Configure re-projection function
        self.trans = pyproj.Proj(self.humDat['epsg'])

        # Reproject lat/lon to UTM zone
        e, n = self.trans(df['lon'], df['lat'])
        df['e'] = e
        df['n'] = n

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
        
        print(df['frequency'])
        
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

   

