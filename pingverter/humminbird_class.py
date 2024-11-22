
import os, sys
import numpy as np
import pandas as pd
import pyproj

from .lowrance_class import low

class hum(object):

    #===========================================================================
    def __init__(self,humFile: str):
        
        self.humFile = humFile
        self.sonFile = humFile.split('.')[0]

        self.head_start_val = 3235818273
        self.head_end_val = 33

        print(self.sonFile)

        return
    
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






        




        