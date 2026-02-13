import os
import struct

import numpy as np
import pandas as pd
import pyproj


class jsf(object):

    def __init__(self, inFile: str, nchunk: int = 0, exportUnknown: bool = False):
        self.humFile = None
        self.isOnix = 0
        self.sonFile = inFile
        self.nchunk = nchunk
        self.exportUnknown = exportUnknown

        self.file_header_size = 0
        self.msg_header_size = 16
        self.msg80_header_size = 240

        self.humDat = {}
        self.headBytes = self.msg_header_size + self.msg80_header_size
        self.son8bit = False
        self.sample_dtype = '<u2'

    def _getFileLen(self):
        self.file_len = os.path.getsize(self.sonFile)
        return

    def _parseFileHeader(self):
        self.file_header = {
            'format': 'JSF',
            'message_header_size': self.msg_header_size,
            'message80_header_size': self.msg80_header_size,
            'endianness': 'little',
        }

        out_file = os.path.join(self.metaDir, 'DAT_meta.csv')
        pd.DataFrame.from_dict(self.file_header, orient='index').T.to_csv(out_file, index=False)
        self.datMetaFile = out_file
        return

    def _parsePingHeader(self):
        rows = []

        with open(self.sonFile, 'rb') as file:
            file_len = self.file_len
            i = self.file_header_size

            while i + self.msg_header_size <= file_len:
                file.seek(i)
                msg_head = file.read(self.msg_header_size)
                if len(msg_head) < self.msg_header_size:
                    break

                marker = struct.unpack_from('<H', msg_head, 0)[0]
                if marker != 0x1601:
                    i += 1
                    continue

                protocol_version = msg_head[2]
                message_type = struct.unpack_from('<H', msg_head, 4)[0]
                subsystem_number = msg_head[7]
                channel = msg_head[8]
                msg_size = struct.unpack_from('<i', msg_head, 12)[0]

                if msg_size <= 0:
                    i += 1
                    continue

                next_i = i + self.msg_header_size + msg_size
                if next_i > file_len:
                    break

                if message_type == 80 and msg_size >= self.msg80_header_size:
                    file.seek(i + self.msg_header_size)
                    msg80 = file.read(self.msg80_header_size)
                    rows.append(self._decode_msg80(i, protocol_version, subsystem_number, channel, msg80))

                i = next_i

        df = pd.DataFrame.from_dict(rows)

        if len(df) == 0:
            raise ValueError('No JSF Message Type 80 records were parsed.')

        df = self._doUnitConversion(df)
        df.sort_values(by=['time_s', 'beam'], inplace=True)
        df.reset_index(drop=True, inplace=True)
        df = self._calcTrkDistTS(df)
        df['record_num'] = np.arange(len(df), dtype=np.int64)

        out_test = os.path.join(self.metaDir, 'All-JSF-Sonar-MetaData.csv')
        df.to_csv(out_test, index=False)

        self.header_dat = df
        return

    def _decode_msg80(self, record_start, protocol_version, subsystem_number, channel, msg80):
        time_since_1970 = struct.unpack_from('<i', msg80, 0)[0]
        ping_number = struct.unpack_from('<I', msg80, 8)[0]
        msb1 = struct.unpack_from('<H', msg80, 16)[0]
        lsb1 = struct.unpack_from('<H', msg80, 18)[0]
        lsb2 = struct.unpack_from('<H', msg80, 20)[0]
        validity_flag = struct.unpack_from('<H', msg80, 30)[0]
        data_format = struct.unpack_from('<h', msg80, 34)[0]

        longitude_raw = struct.unpack_from('<i', msg80, 80)[0]
        latitude_raw = struct.unpack_from('<i', msg80, 84)[0]
        coord_units = struct.unpack_from('<h', msg80, 88)[0]

        samples_lsb = struct.unpack_from('<H', msg80, 114)[0]
        samples_msb = (msb1 >> 8) & 0x0F
        samples = int(samples_lsb + (samples_msb << 16))

        sample_interval_ns = struct.unpack_from('<I', msg80, 116)[0]

        start_freq_dahz = struct.unpack_from('<H', msg80, 126)[0]
        end_freq_dahz = struct.unpack_from('<H', msg80, 128)[0]

        depth_mm = struct.unpack_from('<i', msg80, 136)[0]
        altitude_mm = struct.unpack_from('<i', msg80, 144)[0]
        sound_speed = struct.unpack_from('<f', msg80, 148)[0]

        weighting_factor = struct.unpack_from('<h', msg80, 168)[0]
        compass_heading = struct.unpack_from('<H', msg80, 172)[0]
        pitch_raw = struct.unpack_from('<h', msg80, 174)[0]
        roll_raw = struct.unpack_from('<h', msg80, 176)[0]

        course_tenths = struct.unpack_from('<h', msg80, 192)[0]
        speed_tenths_knots = struct.unpack_from('<h', msg80, 194)[0]

        milli_seconds_today = struct.unpack_from('<I', msg80, 200)[0]
        water_temp_tenths = struct.unpack_from('<h', msg80, 226)[0]
        layback_m = struct.unpack_from('<f', msg80, 228)[0]
        cable_out_dm = struct.unpack_from('<H', msg80, 236)[0]

        ints_per_sample = 2 if data_format in [1, 9] else 1
        ping_cnt = int(samples * ints_per_sample)

        son_offset = self.msg_header_size + self.msg80_header_size

        beam = self._map_beam(subsystem_number, channel)

        lat, lon, e, n = self._decode_position(latitude_raw, longitude_raw, coord_units)

        heading = float(compass_heading) / 100.0
        pitch = float(pitch_raw) * 180.0 / 32768.0
        roll = float(roll_raw) * 180.0 / 32768.0
        course_frac = float((lsb1 >> 8) & 0xFF) / 100.0
        course = float(course_tenths) / 10.0 + course_frac

        speed_frac = float(lsb2 & 0x0F) / 100.0
        speed_kn = float(speed_tenths_knots) / 10.0 + speed_frac

        if sound_speed <= 0 or not np.isfinite(sound_speed):
            sound_speed = 1500.0

        pix_m = np.nan
        if sample_interval_ns > 0:
            pix_m = (float(sound_speed) * (sample_interval_ns / 1e9)) / 2.0

        start_freq_khz = start_freq_dahz / 100.0 if start_freq_dahz > 0 else np.nan
        end_freq_khz = end_freq_dahz / 100.0 if end_freq_dahz > 0 else np.nan
        if np.isfinite(start_freq_khz) and np.isfinite(end_freq_khz):
            f = (start_freq_khz + end_freq_khz) / 2.0
        elif np.isfinite(start_freq_khz):
            f = start_freq_khz
        elif np.isfinite(end_freq_khz):
            f = end_freq_khz
        else:
            f = np.nan

        dep_m = np.nan
        if depth_mm > 0 and (validity_flag & (1 << 9)):
            dep_m = depth_mm / 1000.0

        altitude_m = np.nan
        if altitude_mm > 0 and (validity_flag & (1 << 6)):
            altitude_m = altitude_mm / 1000.0

        time_s = float(time_since_1970)
        if milli_seconds_today > 0:
            frac = (milli_seconds_today % 1000) / 1000.0
            time_s = float(time_since_1970) + frac

        row = {
            'index': int(record_start),
            'son_offset': int(son_offset),
            'protocol_version': int(protocol_version),
            'message_type': 80,
            'subsystem_number': int(subsystem_number),
            'channel': int(channel),
            'beam': int(beam),
            'ping_number': int(ping_number),
            'time_s': float(time_s),
            'ping_cnt': int(ping_cnt),
            'data_format': int(data_format),
            'weighting_factor': int(weighting_factor),
            'f': float(f) if np.isfinite(f) else np.nan,
            'f_min': float(start_freq_khz) if np.isfinite(start_freq_khz) else np.nan,
            'f_max': float(end_freq_khz) if np.isfinite(end_freq_khz) else np.nan,
            'pixM': float(pix_m) if np.isfinite(pix_m) else np.nan,
            'instr_heading': float(heading),
            'pitch': float(pitch),
            'roll': float(roll),
            'yaw': np.nan,
            'course': float(course),
            'speed_ms': float(speed_kn) * 0.514444,
            'inst_dep_m': float(dep_m) if np.isfinite(dep_m) else np.nan,
            'dep_m': float(dep_m) if np.isfinite(dep_m) else np.nan,
            'altitude': float(altitude_m) if np.isfinite(altitude_m) else np.nan,
            'lat': lat,
            'lon': lon,
            'e': e,
            'n': n,
            'validity_flag': int(validity_flag),
            'water_temp_raw': int(water_temp_tenths),
            'layback_m': float(layback_m),
            'cable_out_m': float(cable_out_dm) / 10.0,
            'transect': 0,
        }

        return row

    def _decode_position(self, latitude_raw: int, longitude_raw: int, coord_units: int):
        if coord_units == 2:
            lat = latitude_raw / (10000.0 * 60.0)
            lon = longitude_raw / (10000.0 * 60.0)
            if np.isfinite(lat) and np.isfinite(lon) and abs(lat) <= 90 and abs(lon) <= 180:
                epsg = self._convert_wgs_to_utm(lon, lat)
                self.humDat['epsg'] = f'EPSG:{epsg}'
                self.humDat['wgs'] = 'EPSG:4326'
                self.trans = pyproj.Proj(self.humDat['epsg'])
                e, n = self.trans(lon, lat)
                return float(lat), float(lon), float(e), float(n)

        if coord_units == 1:
            x_m = longitude_raw / 1000.0
            y_m = latitude_raw / 1000.0
        elif coord_units == 3:
            x_m = longitude_raw / 10.0
            y_m = latitude_raw / 10.0
        elif coord_units == 4:
            x_m = longitude_raw / 100.0
            y_m = latitude_raw / 100.0
        else:
            x_m = np.nan
            y_m = np.nan

        self.humDat['epsg'] = 'UNKNOWN'
        self.humDat['wgs'] = 'EPSG:4326'
        self.trans = lambda lon, lat: (lon, lat)
        return np.nan, np.nan, float(x_m) if np.isfinite(x_m) else np.nan, float(y_m) if np.isfinite(y_m) else np.nan

    def _convert_wgs_to_utm(self, lon: float, lat: float):
        utm_band = str(int((np.floor((lon + 180) / 6) % 60) + 1))
        if len(utm_band) == 1:
            utm_band = '0' + utm_band
        if lat >= 0:
            return '326' + utm_band
        return '327' + utm_band

    def _map_beam(self, subsystem_number: int, channel: int):
        if channel == 0:
            return 2
        if channel == 1:
            return 3

        if subsystem_number == 0:
            return 1
        return 4

    def _doUnitConversion(self, df: pd.DataFrame):
        if 'inst_dep_m' not in df.columns:
            if 'dep_m' in df.columns:
                df['inst_dep_m'] = df['dep_m']
            else:
                df['inst_dep_m'] = np.nan

        if 'dep_m' not in df.columns:
            df['dep_m'] = df['inst_dep_m']

        df['tempC'] = np.float32(self.tempC * 10)

        if 'water_temp_raw' in df.columns:
            is_valid_temp = (df['validity_flag'] & (1 << 8)) > 0
            df.loc[is_valid_temp, 'tempC'] = df.loc[is_valid_temp, 'water_temp_raw'] / 10.0

        return df

    def _calcTrkDistTS(self, df: pd.DataFrame):
        ts = df['time_s'].to_numpy(dtype=float)
        ss = df['speed_ms'].fillna(0).to_numpy(dtype=float)

        if len(ts) == 0:
            df['trk_dist'] = []
            return df

        ds = np.zeros((len(ts),), dtype=float)
        if len(ts) > 1:
            d = np.maximum(0, (ts[1:] - ts[:-1]) * ss[1:])
            ds[1:] = d
            ds = np.cumsum(ds)

        df['trk_dist'] = ds
        return df

    def _recalcRecordNum(self):
        df = self.header_dat.reset_index(drop=True)
        df['record_num'] = df.index
        self.header_dat = df
        return

    def _splitBeamsToCSV(self):
        self.beamMeta = beamMeta = {}
        df = self.header_dat

        for beam, group in df.groupby('beam'):
            meta = {}

            if beam in [2, 3] and 'pixM' in group.columns and len(group) > 0:
                self.pixM = group['pixM'].iloc[0]

            beam_name = f'B00{int(beam)}'
            meta['beamName'] = self._getBeamName(beam_name)
            meta['sonFile'] = self.sonFile

            group = self._getChunkID(group.copy())

            out_csv = f'{beam_name}_{meta["beamName"]}_meta.csv'
            out_csv = os.path.join(self.metaDir, out_csv)
            group.to_csv(out_csv, index=False)

            meta['metaCSV'] = out_csv
            beamMeta[beam_name] = meta

        return

    def _getBeamName(self, beam: str):
        if beam == 'B000':
            return 'ds_lowfreq'
        if beam == 'B001':
            return 'ds_highfreq'
        if beam == 'B002':
            return 'ss_port'
        if beam == 'B003':
            return 'ss_star'
        if beam == 'B004':
            return 'ds_vhighfreq'
        return 'unknown'

    def _getChunkID(self, df: pd.DataFrame):
        df.reset_index(drop=True, inplace=True)
        df['chunk_id'] = int(-1)

        chunk = 0
        start_idx = 0
        end_idx = self.nchunk

        while start_idx < len(df):
            df.iloc[start_idx:end_idx, df.columns.get_loc('chunk_id')] = int(chunk)
            chunk += 1
            start_idx = end_idx
            end_idx += self.nchunk

        last_chunk = df[df['chunk_id'] == chunk]
        if len(last_chunk) <= self.nchunk / 2:
            df.loc[df['chunk_id'] == chunk, 'chunk_id'] = chunk - 1

        return df
