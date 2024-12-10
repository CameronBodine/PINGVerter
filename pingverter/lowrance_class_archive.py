
import os
import numpy as np
import pandas as pd

'''
Based on Sonarlight by Kenneth Thorø Martinsen
The package is inspired by and builds upon other tools and descriptions for processing Lowrance sonar data, e.g. SL3Reader which includes a usefull paper, python-sllib, sonaR, Navico_SLG_Format notes, older blog post.
'''

# #dtype for '.sl2' files (144 bytes)
# sl2Struct = np.dtype([
#     ("first_byte", "<u4"),
#     ("frame_version", "<u4"),
#     ("unknown8", "<f4"),
#     ("unknown12", "<f4"),
#     ("unknown16", "<f4"),
#     ("unknown20", "<f4"),
#     ("unknown24", "<f4"),
#     ("frame_size", "<u2"),
#     ("prev_frame_size", "<u2"),
#     ("survey_type", "<u2"),
#     ("packet_size", "<u2"),
#     ("id", "<u4"),
#     ("min_range", "<f4"),
#     ("max_range", "<f4"),
#     ("unknown48", "<f4"),
#     ("unknown52", "<u1"),
#     ("frequency_type", "<u2"),
#     ("unknown55", "<u1"),
#     ("unknown56", "<f4"),
#     ("hardware_time", "<u4"),
#     ("water_depth", "<f4"),
#     ("unknown68", "<f4"),
#     ("unknown72", "<f4"),
#     ("unknown76", "<f4"),
#     ("unknown80", "<f4"),
#     ("unknown84", "<f4"),
#     ("unknown88", "<f4"),
#     ("unknown92", "<f4"),
#     ("unknown96", "<f4"),
#     ("gps_speed", "<f4"),
#     ("water_temperature", "<f4"),
#     ("x", "<i4"),
#     ("y", "<i4"),
#     ("water_speed", "<f4"),
#     ("gps_heading", "<f4"),
#     ("gps_altitude", "<f4"),
#     ("magnetic_heading", "<f4"),
#     ("flags", "<u2"),
#     ("unknown132", "<u2"),
#     ("unknown136", "<f4"),
#     ("seconds", "<u4")
# ])

# #dtype for '.sl3' files (168 bytes)
# sl3Struct = np.dtype([
#     ("first_byte", "<u4"),
#     ("frame_version", "<u4"),
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
#     ("echo_size", "<u4"),
#     ("water_depth", "<f4"),
#     ("frequency_type", "<u2"),
#     ("unknown54", "<f4"),
#     ("unknown58", "<f4"),
#     ("unknown62", "<i2"),
#     ("unknown64", "<f4"),
#     ("unknown68", "<f4"),
#     ("unknown72", "<f4"),
#     ("unknown76", "<f4"),
#     ("unknown80", "<f4"),
#     ("gps_speed", "<f4"),
#     ("water_temperature", "<f4"),
#     ("x", "<i4"),
#     ("y", "<i4"),
#     ("water_speed", "<f4"),
#     ("gps_heading", "<f4"),
#     ("gps_altitude", "<f4"),
#     ("magnetic_heading", "<f4"),
#     ("flags", "<u2"),
#     ("unknown118", "<u2"),
#     ("unknown120", "<u4"),
#     ("seconds", "<u4"), #milliseconds
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


# For Lowrance to Humminbird proof of concept...
header = {0: [0, 0, 2, 'format', '<u2'],
          2: [2, 0, 2, 'version', '<u2'],
          4: [4, 0, 2, 'bytes_per_sounding', '<u2'],
          6: [6, 0, 1, 'debug', 'B'],
          7: [7, 0, 1, 'byte', 'B']}

# For Lowrance to Humminbird proof of concept...
# {offset: [offset, offset from current position (always 0 for Lowrance), length, name]}
sl2Struct_forHum = {0: [0, 0, 4, 'frame_offset', '<u4'], 
             4: [4, 0, 4, 'prev_primary_offset', '<u4'],
             8: [8, 0, 4, 'prev_secondary_offset', '<u4'],
             12: [12, 0, 4, 'prev_downscan_offset', '<u4'], 
             16: [16, 0, 4, 'prev_left_sidescan_offset', '<u4'],
             20: [20, 0, 4, 'prev_right_sidescan_offset', '<u4'],
             24: [24, 0, 4, 'prev_sidescan_offset', '<u4'],
             28: [28, 0, 2, 'frame_size', '<u2'],
             30: [30, 0, 2, 'prev_frame_size', '<u2'],
             32: [32, 0, 2, 'channel_type', '<u2'],
             34: [34, 0, 2, 'packet_size', '<u2'], # ping count
             36: [36, 0, 4, 'frame_index', '<u4'],
             40: [40, 0, 4, 'min_range', '<f4'],
             44: [44, 0, 4, 'max_range', '<f4'],
             48: [44, 0, 2, 'unknown48', '<u2'],
             50: [50, 0, 1, 'unknown50', 'B'],
             51: [51, 0, 1, 'unknown51', 'B'],
             52: [52, 0, 1, 'unknown52', 'B'],
             53: [53, 0, 1, 'frequency', 'B'], # Frequency
             54: [54, 0, 2, 'unknown54', '<u2'],
             56: [56, 0, 2, 'unknown56', '<u2'],
             58: [58, 0, 2, 'unknown58', '<u2'],
             60: [60, 0, 4, 'creation_date_time', '<u4'],
             64: [64, 0, 4, 'depth_ft', '<f4'],
             68: [68, 0, 4, 'keel_depth_ft', '<f4'],
             72: [72, 0, 1, 'unkown72', 'B'],
             73: [73, 0, 1, 'unknown73', 'B'],
             74: [74, 0, 2, 'unknown74', '<u2'],
             76: [76, 0, 1, 'unknown76', 'B'],
             77: [77, 0, 1, 'unknown77', 'B'],
             78: [78, 0, 2, 'unknown78', '<u2'],
             80: [80, 0, 4, 'unknown80', '<f4'],
             84: [84, 0, 4, 'unknown84', '<f4'],
             88: [88, 0, 4, 'unknown88', '<f4'],
             92: [92, 0, 4, 'unknown92', '<f4'],
             96: [96, 0, 1, 'unknown96', 'B'],
             97: [97, 0, 1, 'unknown97', 'B'],
             98: [98, 0, 1, 'unknown98', 'B'],
             99: [99, 0, 1, 'unknown99', 'B'],
             100: [100, 0, 4, 'gps_speed', '<f4'], #[knots]
             104: [104, 0, 4, 'water_temperature', '<f4'], #[C]
             108: [108, 0, 4, 'utm_e', '<i4'], #Easting in mercator [meters]
             112: [112, 0, 4, 'utm_n', '<i4'], #Northing in mercator [meters]
             116: [116, 0, 4, 'water_speed', '<f4'], #Water speed through paddlewheel or GPS if not present [knots]
             120: [120, 0, 4, 'track_cog', '<f4'], # Track (COG) [radians]
             124: [124, 0, 4, 'altitude', '<f4'], # Above sea level [feet]
             128: [128, 0, 4, 'heading', '<f4'], #[radians]
             132: [132, 0, 2, 'flags', '<u2'],
             134: [134, 0, 2, 'unknown134', '<u2'],
             136: [136, 0, 1, 'unknown136', 'B'],
             137: [137, 0, 1, 'unknown137', 'B'],
             138: [138, 0, 1, 'unknown138', 'B'],
             139: [139, 0, 1, 'unknown139', 'B'],
             140: [140, 0, 4, 'time', '<u4'], # Time since beginning of log [ms]
             } # sonar returns begin at 144 

class low(object):


    # Below used for Lowrance to Humminbird proof of concept...
    def __init__(self, path: str):
        self.path = path
        self.file_header_size = 8
        self.extension = os.path.basename(path).split('.')[-1]
        self.header = header

        self.frame_header_size = 168 if "sl3" in self.extension else 144
        self.son_struct = sl3Struct if "sl3" in self.extension else sl2Struct

        self.supported_channels = ["primary", "secondary", "downscan", "sidescan"]
        self.valid_channels = []
        self.valid_channels_records = []
        
        self.survey_dict = {0: 'primary', 1: 'secondary', 2: 'downscan',
                            3: 'left_sidescan', 4: 'right_sidescan', 5: 'sidescan',
                            9: '3D', 10: 'debug_digital', 11: 'debug_noise'}
        
        self.frequency_dict = {0: "200kHz", 1: "50kHz", 2: "83kHz",
                               3: "455kHz", 4: "800kHz", 5: "38kHz", 
                               6: "28kHz", 7: "130kHz_210kHz", 8: "90kHz_150kHz", 
                               9: "40kHz_60kHz", 10: "25kHz_45kHz"}
        
        self.vars_to_keep = ["id", "survey", "datetime",
                             "x", "y", "longitude", "latitude", 
                             "min_range", "max_range", "water_depth", 
                             "gps_speed", "gps_heading", "gps_altitude", 
                             "bottom_index", "frames"]

        return
    
    def _fread(self,
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
        # dat = arr(typ)
        # dat.fromfile(infile, num)
        # return(list(dat))

        buffer = infile.read(num)
        data = np.frombuffer(buffer, dtype=typ)

        return data
    
    def _parseFileHeader(self):
        # Open sonar log
        f = open(self.path, 'rb')

        # Iterate known file header items
        header = dict()
        for k, v in self.header.items():
            offset = v[0]
            length = v[2]
            name = v[3]
            type = v[4]
            f.seek(offset)

            v = self._fread(f, length, type)
            header[name] = v.item()

        # Set class attribtutes
        self.file_header = header

        return

    def _getFileLen(self):
        self.file_len = os.path.getsize(self.path)

    def _parsePingHeader(self):

        # Initialize offset after file header
        i = self.file_header_size

        # Decode ping header
        while i < self.file_len:

            # Get header data at offset i
            header_dat = self._getPingHeader(i)

            # Set offset for next frame
            i = header_dat['frame_offset'][0] + header_dat['frame_size'][0]

            # Store in global dictionary
            if 'header_dat_all' not in locals():
                header_dat_all = header_dat
            else:
                for k, v in header_dat.items():
                    header_dat_all[k] += v

            del header_dat

        # Store in class attribute as dataframe
        self.header_dat = pd.DataFrame.from_dict(header_dat_all)

        # self.header_dat.to_csv('lowrance_test.csv')

        return
    
    def _getPingHeader(self, i):

        # Get necessary attributes
        head_struct = self.son_struct

        # Open sonar file
        file = open(self.path, 'rb')

        # For storing header contents
        son_head = dict()

        for k, v in head_struct.items():
            byte_index = v[0] # Offset within header
            offset_from_byte = v[1] # Additional offset from byte_index (Lowrance always 0)
            length = v[2] # Length of data to be decoded
            struct_name = v[3] # Name of attribute
            struct_type = v[4] # Data type to be docoded

            # Calculate index of global offset
            index = i + byte_index + offset_from_byte

            # Move to location in file
            file.seek(index)

            # Get the data
            byte = self._fread(file, length, struct_type)

            # Store the data
            son_head[struct_name] = [byte.item()]

        return son_head
    
    def _convertPingAttributes(self):
        '''
        Convert ping attributes (headers) to known units
        '''

        df = self.header_dat



        df[["depth_ft", "min_range", "max_range", "altitude"]] /= 3.2808399 #feet to meter
        df["gps_speed"] *=  0.5144 #knots to m/s
        df["survey"] = [self.survey_dict.get(i, "unknown") for i in df["channel_type"]]
        df["frequency"] = [self.frequency_dict.get(i, "unknown") for i in df["frequency"]]
        df["time"] /= 1000 #milliseconds to seconds
        hardware_time_start = df["creation_date_time"][0]
        df["datetime"] = pd.to_datetime(hardware_time_start+df["time"], unit='s')

        self.header_dat = df
        return

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
