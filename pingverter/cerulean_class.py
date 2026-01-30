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

Copyright (c) 2025 Cameron S. Bodine

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

'''
Notes:

I think sonar time is based on the os_uptime and;
mavlink (navigation) is using process_uptime. Annoying!
'''

import os, sys
import numpy as np
import json
import pandas as pd
from datetime import datetime
import pyproj

# Structure is of Blue Robotis Ping Protocol: https://github.com/bluerobotics/ping-protocol
# Documented at Cerulean: https://docs.ceruleansonar.com/c/cerulean-ping-protocol

packetHeadStruct = np.dtype([
    ("B", "<u1"),
    ("R", "<u1"),
    ("packet_len", "<u2"),
    ("packet_id", "<u2"),
    ("SP6", "<u1"),
    ("SP7", "<u1")
])

svlogStruct = np.dtype([
    ("ping_number", "<u4"),
    ("start_mm", "<u4"),
    ("length_mm", "<u4"),
    ("timestamp_ms", "<u4"),
    ("ping_hz", "<u4"),
    ("gain_index", "<u2"),
    ("num_results", "<u2"),
    ("sos_dmps", "<u2"),
    ("channel_number", "<u1"),
    ("SP27", "<u1"),
    ("pulse_duration_sec", "<f4"),
    ("analog_gain", "<f4"),
    ("max_pwr_db", "<f4"),
    ("min_pwr_db", "<f4"),
    ("transducer_heading_deg", "<f4"),
    ("vehicle_heading_deg", "<f4"),
])

cerulCols2PM = {
    'hdg': 'instr_heading',
    'alt': 'altitude',
    'ping_number': 'record_num',
    'num_results': 'ping_cnt',

}

class cerul(object):

    #===========================================================================
    def __init__(self, svlog: str, nchunk: int=0, exportUnknown: bool=False, port=0, star=1):
        '''
        '''

        self.humFile = None
        self.isOnix = 0
        self.sonFile = svlog
        self.nchunk = nchunk
        self.exportUnknown = exportUnknown

        self.packet_header_size = 8
        self.checksum_len = 2

        self.packetHeadStruct = packetHeadStruct
        self.son_struct = svlogStruct
        self.headBytes = 52

        self.humDat = {} # Store general sonar recording metadata

        self.cerulCols2PM = cerulCols2PM

        self.port = port
        self.star = star

        self.son8bit = False

        return
    
    #===========================================================================
    def _getFileLen(self):
        self.file_len = os.path.getsize(self.sonFile)

        return
    
    #===========================================================================
    def _parseFileHeader(self):
        '''
        '''
        # Get necessary attributes
        packet_struct = self.packetHeadStruct
        length = self.packet_header_size
        checksum_len = self.checksum_len

        # Open sonar log
        file = open(self.sonFile, 'rb')

        # Get packet header
        packet_head, _ = self._getPacketHeader(file, 0)

        # Set the file header
        self.file_header_size = packet_head['packet_len'] + length + checksum_len

        # If json, do conversion
        if packet_head['packet_id'] == 10:

            packet = self._getJSONdat(file, location=length, length=packet_head['packet_len'])

            for k,v in packet.items():

                print('\n\n', k, v)

            # Store time variables
            self.hardware_time_start = datetime.fromisoformat(packet['timestamp']).timestamp()

            # GPS sensor appears to use the process_uptime
            ## Substrace from all gps readings then add to timestamp to get absolute date/time
            self.nav_time_init = packet['process_uptime'] * 1000 # sec to ms

            # Sonar sensor appears to use the process_uptime
            ## Substrace from all sonar readings then add to timestamp to get absolute date/time
            try:
                self.sonar_time_init = packet['os_uptime'] * 1000 # sec to ms
            except:
                self.sonar_time_init = packet['session_uptime'] * 1000 # sec to ms
                


            # Set class attributes
            self.file_header = packet

        return

    
    #===========================================================================
    def _getJSONdat(self, file, location: int, length: int):
        '''
        '''
        # Move to file location
        file.seek(location)

        # Get file header contents
        buffer = file.read(length)

        # Convert buffer to string
        string_val = buffer.decode("utf-8")

        json_val = json.loads(string_val)

        return json_val
    
    #===========================================================================
    def _locatePacketsRaw(self):
        '''
        '''

        nav_time_name = 'time_boot_ms'
        son_time_name = 'timestamp_ms'

        # Get the file length
        file_len = self.file_len

        # 
        headBytes = self.headBytes
        son_struct = self.son_struct

        # Initialize offset after file header
        i = self.file_header_size

        # Open the file
        file = open(self.sonFile, 'rb')

        # Store contents in list
        header_dat_all = []

        # Store data from current time offset
        cur_nav_dat = {}

        # counter for testing
        test_cnt = 0

        while i < file_len:

            # Get packet data at offset i
            header_dat, cpos = self._getPacketHeader(file, i)

            # If json, do conversion
            if header_dat['packet_id'] == 150:
                packet_dat = {}
                packet = self._getJSONdat(file, location=cpos, length=header_dat['packet_len'])

                for k,v in packet['header'].items():
                    packet_dat[k] = v
                for k,v in packet['message'].items():
                    packet_dat[k] = v

                header_dat_all.append(packet_dat)

            # If sonar data, get ping attributes
            if header_dat['packet_id'] == 2198:
                
                # Move to offset
                file.seek(cpos)
                
                # Get the data
                buffer = file.read(headBytes)

                # Read the data
                header = np.frombuffer(buffer, dtype=son_struct)

                # Populate dictionary
                packet_dat = {}
                for name, typ in header.dtype.fields.items():
                    packet_dat[name] = header[name][0].item()

                header_dat_all.append(packet_dat)


            i = cpos + header_dat['packet_len'] + self.checksum_len



            test_cnt += 1
            # if test_cnt == 100:
            #     break

        # Convert to dataframe
        df = pd.DataFrame.from_dict(header_dat_all)

        # Save raw data. Does not include anything that didn't have a time reported.
        outCSV = 'All-Cerulean-Sonar-MetaData-RAW.csv'
        outCSV = os.path.join(self.metaDir, outCSV)
        df.to_csv(outCSV, index=False)

        return
    
    #===========================================================================
    def _locatePackets(self):
        '''
        '''

        nav_time_name = 'time_boot_ms'
        son_time_name = 'timestamp_ms'

        # Get the file length
        file_len = self.file_len

        # 
        headBytes = self.headBytes
        son_struct = self.son_struct

        # Initialize offset after file header
        i = self.file_header_size

        # Open the file
        file = open(self.sonFile, 'rb')

        # Store contents in list
        header_dat_all = []

        # Store data from current time offset
        cur_nav_dat = {}

        # counter for testing
        test_cnt = 0

        while i < file_len:

            # Get packet data at offset i
            header_dat, cpos = self._getPacketHeader(file, i)

            # If json, do conversion
            if header_dat['packet_id'] == 150:
                packet_dat = {}
                packet = self._getJSONdat(file, location=cpos, length=header_dat['packet_len'])

                found_time = False

                for k,v in packet['header'].items():
                    packet_dat[k] = v
                for k,v in packet['message'].items():
                    packet_dat[k] = v
                    if nav_time_name in k:
                        found_time = True

                if found_time:
                    if len(cur_nav_dat) == 0:
                        cur_nav_dat = packet_dat
                    else:
                        cur_time = cur_nav_dat[nav_time_name]
                        nex_time = packet_dat[nav_time_name]

                        if cur_time == nex_time:
                            for k,v in packet_dat.items():
                                cur_nav_dat[k] = v

                        else:
                            header_dat_all.append(cur_nav_dat)
                            cur_nav_dat = packet_dat


                if found_time:
                    # Calculate time offset
                    packet_dat['time_s'] = (packet_dat[nav_time_name] - self.nav_time_init) / 1000
                    packet_dat['index'] = i

            # If sonar data, get ping attributes
            if header_dat['packet_id'] == 2198:

                found_time = False
                
                # Move to offset
                file.seek(cpos)
                
                # Get the data
                buffer = file.read(headBytes)

                # Read the data
                header = np.frombuffer(buffer, dtype=son_struct)

                # Populate dictionary
                packet_dat = {}
                for name, typ in header.dtype.fields.items():
                    packet_dat[name] = header[name][0].item()
                    if son_time_name in name:
                        found_time = True

                if found_time:
                    # Calculate time offset
                    packet_dat['time_s'] = (packet_dat[son_time_name] - self.sonar_time_init) / 1000
                    packet_dat['index'] = i

                    header_dat_all.append(packet_dat)

            # if found_time:
            #     header_dat_all.append(packet_dat)


            i = cpos + header_dat['packet_len'] + self.checksum_len



            test_cnt += 1
            # if test_cnt == 100:
            #     break

        # Convert to dataframe
        df = pd.DataFrame.from_dict(header_dat_all)

        # Do interpolation of position/imu information for each ping
        df = self._doPosInterp(df)

        # Drop nan's
        df = df.dropna(subset=['ping_number','lat', 'lon']).reset_index(drop=True)

        # Do unit conversion
        df = self._doUnitConversion(df)

        # Do column conversions to PING-Mapper column names
        df.rename(columns=self.cerulCols2PM, inplace=True)

        # Drop non-sonar columns
        df = df[df['ping_cnt'].notna()]

        # Calculate vessel speed
        df = self._calcSpeedDist(df)

        # # Caclculate cog
        # df = self._calcCOG(df)

        # Test file to see outputs
        out_test = os.path.join(self.metaDir, 'All-Cerulean-Sonar-MetaData.csv')
        df.to_csv(out_test, index=False)

        self.header_dat = df

        return
    
    # ======================================================================
    def _getPacketHeader(self, file, i: int):

        # Get necessary attributes
        packet_struct = self.packetHeadStruct
        length = self.packet_header_size

        # Move to offset
        file.seek(i)

        # Get the data
        buffer = file.read(length)

        # Read the data
        header = np.frombuffer(buffer, dtype=packet_struct)

        out_dict = {}
        for name, typ in header.dtype.fields.items():
            out_dict[name] = header[name][0].item()

        return out_dict, file.tell()


    # ======================================================================
    def _doPosInterp(self, df: pd.DataFrame):
        '''
        '''

        field2Interp = ['pitch', 'pitchspeed', 'roll', 'rollspeed', 'time_boot_ms',
                        'yaw', 'yawspeed', 'alt', 'hdg', 'lat', 'lon', 'relative_alt',
                        'vx', 'vy', 'vz', 'x', 'y', 'z',]
        
        for f in field2Interp:
            try:
                df[f] = df[f].interpolate(method='linear')
            except:
                pass

        return df

    # ======================================================================
    def _doUnitConversion(self, df: pd.DataFrame):
        '''
        '''

        # Calculate speed (m/s)

        # Convert time to timestamp
        df['caltime'] = pd.to_datetime(self.hardware_time_start + df['time_s'], unit='s')
        df['date'] = df['caltime'].dt.date
        df['time'] = df['caltime'].dt.time

        df.drop(columns=['caltime'], inplace=True)

        # Convert lat lon to decimal
        df['lat'] = df['lat'] * 1e-7
        df['lon'] = df['lon'] * 1e-7

        # Determine epsg code 
        self.humDat['epsg'] = "epsg:"+str(int(float(self._convert_wgs_to_utm(df['lon'][0], df['lat'][0]))))
        self.humDat['wgs'] = "epsg:4326"

        # Configure re-projection function
        self.trans = pyproj.Proj(self.humDat['epsg'])

        # Reproject lat/lon to UTM zone
        e, n = self.trans(df['lon'], df['lat'])
        df['e'] = e
        df['n'] = n

        # Convert heading
        df['hdg'] = df['hdg'] * 1e-2

        # Convert altitude (?)
        df['alt'] = df['alt'] * 1e-3
        df['relative_alt'] = df['relative_alt'] * 1e-3

        # Store survey temperature
        df['tempC'] = self.tempC*10

        # Add transect number (for aoi processing)
        df['transect'] = 0

        # Calculate min/max range
        df['min_range'] = df['start_mm'] * 1e-3
        df['max_range'] = df['length_mm'] * 1e-3

        # Calculate pixel size [m]
        df['pixM'] = (df['max_range'] - df['min_range']) / (df['num_results'])

        # Calculate frequency
        df['frequency'] = df['ping_hz'] * 1e-3

        # Calculate offset to sonar
        df['son_offset'] = self.headBytes

        return df
    
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

    # ======================================================================
    def _calcSpeedDist(self, df: pd.DataFrame):
        '''
        '''

        e = 'e'
        n = 'n'
        tim = 'time_s'

        # Initialize an empty array to store the speed values
        speed_values = np.zeros(len(df))

        for name, group in df.groupby(['channel_number']):

            # Prepare pntA values [0:n-1]
            lonA = group[e].to_numpy() # Store longitude coordinates in numpy array
            latA = group[n].to_numpy() # Store longitude coordinates in numpy array
            lonA = lonA[:-1] # Omit last coordinate
            latA = latA[:-1] # Omit last coordinate
            pntA = [lonA,latA] # Store in array of arrays

            # Prepare pntB values [0+1:n]
            lonB = group[e].to_numpy() # Store longitude coordinates in numpy array
            latB = group[n].to_numpy() # Store longitude coordinates in numpy array
            lonB = lonB[1:] # Omit first coordinate
            latB = latB[1:] # Omit first coordinate
            pntB = [lonB,latB] # Store in array of arrays

            # Calculate time difference
            timeA = group[tim].to_numpy()
            timeA = timeA[:-1]

            timeB = group[tim].to_numpy()
            timeB = timeB[1:]
        
            timeDif = timeB - timeA

            # Calculate distance
            dist = np.sqrt( (pntA[0] - pntB[0])**2 + (pntA[1] - pntB[1])**2 )

            # Calculate meters per second
            mps = dist/timeDif
            last = mps[-1]
            mps = np.append(mps, last)

            # Update the speed values in the original DataFrame
            # speed_values[group.index] = mps
            df.loc[group.index, 'speed_ms'] = mps

            # Add distance
            last = dist[-1]
            dist = np.append(dist, last)
            df.loc[group.index, 'dist_m'] = dist

            # Calculate cumulative distance
            cum_dist = np.cumsum(dist)
            df.loc[group.index, 'trk_dist'] = cum_dist

        return df
    
    # ======================================================================
    # def _calcCOG(self, df: pd.DataFrame):
    #     '''
    #     '''

    #     lat = 'lat'
    #     lon = 'lon'

    #     for name, group in df.groupby(['channel_number']):

    #         # Prepare pntA values [0:n-1]
    #         lonA = group[lon].to_numpy() # Store longitude coordinates in numpy array
    #         latA = group[lat].to_numpy() # Store longitude coordinates in numpy array
    #         lonA = lonA[:-1] # Omit last coordinate
    #         latA = latA[:-1] # Omit last coordinate
    #         pntA = [lonA,latA] # Store in array of arrays

    #         # Prepare pntB values [0+1:n]
    #         lonB = group[lon].to_numpy() # Store longitude coordinates in numpy array
    #         latB = group[lat].to_numpy() # Store longitude coordinates in numpy array
    #         lonB = lonB[1:] # Omit first coordinate
    #         latB = latB[1:] # Omit first coordinate
    #         pntB = [lonB,latB] # Store in array of arrays

    #         # Convert latitude values into radians
    #         lat1 = np.deg2rad(pntA[1])
    #         lat2 = np.deg2rad(pntB[1])

    #         diffLong = np.deg2rad(pntB[0] - pntA[0]) # Calculate difference in longitude then convert to degrees
    #         bearing = np.arctan2(np.sin(diffLong) * np.cos(lat2), np.cos(lat1) * np.sin(lat2) - (np.sin(lat1) * np.cos(lat2) * np.cos(diffLong))) # Calculate bearing in radians

    #         db = np.degrees(bearing) # Convert radians to degrees
    #         db = (db + 360) % 360 # Ensure degrees in range 0-360

    #         last = db[-1]
    #         db = np.append(db, last)
            
    #         df.loc[group.index, 'cog'] = db

    #     return df

    # ======================================================================
    def _convertBeam(self):
        '''
        '''
        df = self.header_dat

        beam_xwalk = {self.port: 2, self.star: 3}

        df['beam'] = [beam_xwalk.get(i, "unknown") for i in df['channel_number']]

        self.header_dat = df

        return 
    
    # ======================================================================
    def _convertFrequency(self):

        '''
        Only one frequency known: 450
        '''
        df = self.header_dat

        df['f'] = 450
        df['f_min'] = 450
        df['f_max'] = 450

        self.header_dat = df

        return

    # ======================================================================
    def _recalcRecordNum(self):

        '''
        '''

        df = self.header_dat

        # Reset index and recalculate record num
        ## Record num is unique for each ping across all sonar beams
        df = df.reset_index(drop=True)
        df['record_num'] = df.index

        self.header_dat = df

        return

    # ======================================================================
    def _splitBeamsToCSV(self):

        '''
        '''

        df = self.header_dat

        # Dictionary to store necessary attributes for PING-Mapper
        self.beamMeta = beamMeta = {}

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
            group.drop(columns=['frequency', 'component_id', 'sequence', 'system_id', 'type'], inplace=True)

            # Add chunk_id
            group = self._getChunkID(group)

            # Save csv
            outCSV = '{}_{}_meta.csv'.format(beam, meta['beamName'])
            outCSV = os.path.join(self.metaDir, outCSV)
            group.to_csv(outCSV, index=False)

            meta['metaCSV'] = outCSV

            # Store the beams metadata
            beamMeta[beam] = meta

        self.header_dat = df

        return
    
    # ======================================================================
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
        output = "Cerulean Class Contents"
        output += '\n\t'
        output += self.__repr__()
        temp = vars(self)
        for item in temp:
            output += '\n\t'
            output += "{} : {}".format(item, temp[item])
        return output
        