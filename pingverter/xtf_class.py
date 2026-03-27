import os
import struct
from datetime import datetime

import numpy as np
import pandas as pd
import pyproj


class xtf(object):

    def __init__(self, inFile: str, nchunk: int = 0, exportUnknown: bool = False):
        self.humFile = None
        self.isOnix = 0
        self.sonFile = inFile
        self.nchunk = nchunk
        self.exportUnknown = exportUnknown

        self.file_header_size = 1024
        self.ping_header_size = 256
        self.ping_chan_header_size = 64

        self.humDat = {}
        self.headBytes = 0
        self.son8bit = False
        self.sample_dtype = '<u2'
        self.nav_units = None

    def _getFileLen(self):
        self.file_len = os.path.getsize(self.sonFile)
        return

    def _parseFileHeader(self):
        with open(self.sonFile, 'rb') as file:
            base = file.read(1024)

        if len(base) < 1024:
            raise ValueError('Invalid XTF file: file header shorter than 1024 bytes.')

        nav_units = struct.unpack_from('<H', base, 164)[0]
        self.nav_units = int(nav_units)
        sonar_channels = struct.unpack_from('<H', base, 166)[0]

        if sonar_channels > 6:
            extra = sonar_channels - 6
            extra_blocks = int(np.ceil(extra / 8.0))
            self.file_header_size = 1024 + (extra_blocks * 1024)
        else:
            self.file_header_size = 1024

        with open(self.sonFile, 'rb') as file:
            header = file.read(self.file_header_size)

        chaninfo = self._parse_chaninfo(header, sonar_channels)

        self.file_header = {
            'nav_units': nav_units,
            'number_of_sonar_channels': sonar_channels,
            'file_header_size': self.file_header_size,
        }
        self.chaninfo = chaninfo

        out_file = os.path.join(self.metaDir, 'DAT_meta.csv')
        pd.DataFrame.from_dict(self.file_header, orient='index').T.to_csv(out_file, index=False)
        self.datMetaFile = out_file

        return

    def _parse_chaninfo(self, header_bytes: bytes, sonar_channels: int):
        chaninfo = {}
        base_offset = 256
        chan_size = 128

        for idx in range(sonar_channels):
            start = base_offset + idx * chan_size
            end = start + chan_size
            if end > len(header_bytes):
                break

            block = header_bytes[start:end]
            bytes_per_sample = struct.unpack_from('<H', block, 6)[0]
            frequency = struct.unpack_from('<f', block, 32)[0]
            sample_format = block[74]
            type_of_channel = block[0]
            sub_channel_number = block[1]

            chaninfo[idx] = {
                'type_of_channel': type_of_channel,
                'sub_channel_number': sub_channel_number,
                'bytes_per_sample': bytes_per_sample,
                'frequency': frequency,
                'sample_format': sample_format,
            }

        return chaninfo

    def _parsePingHeader(self):
        with open(self.sonFile, 'rb') as file:
            file_len = self.file_len
            i = self.file_header_size
            rows = []

            while i + self.ping_header_size <= file_len:
                file.seek(i)
                header = file.read(self.ping_header_size)

                if len(header) < self.ping_header_size:
                    break

                magic = struct.unpack_from('<H', header, 0)[0]
                if magic != 0xFACE:
                    i += 1
                    continue

                header_type = header[2]
                num_chans_to_follow = struct.unpack_from('<H', header, 4)[0]
                num_bytes_this_record = struct.unpack_from('<I', header, 10)[0]

                if num_bytes_this_record <= 0:
                    i += 1
                    continue

                if i + num_bytes_this_record > file_len:
                    break

                if header_type == 0 and num_chans_to_follow > 0:
                    rows.extend(self._parse_sonar_record(i, header, num_chans_to_follow, num_bytes_this_record, file))

                i += num_bytes_this_record

        df = pd.DataFrame.from_dict(rows)

        if len(df) == 0:
            raise ValueError('No XTF sonar ping packets (HeaderType 0) were parsed.')

        df = self._split_combined_sidescan(df)
        df = self._sync_port_star_metadata(df)

        df.sort_values(by=['time_s', 'beam'], inplace=True)
        df.reset_index(drop=True, inplace=True)

        df = self._doUnitConversion(df)
        df = self._calcTrkDistTS(df)

        df['record_num'] = np.arange(len(df), dtype=np.int64)

        out_test = os.path.join(self.metaDir, 'All-XTF-Sonar-MetaData.csv')
        df.to_csv(out_test, index=False)

        self.header_dat = df
        return

    def _split_combined_sidescan(self, df: pd.DataFrame):
        return self._split_combined_sidescan_group(df)

    def _split_combined_sidescan_group(self, df: pd.DataFrame):
        beams = set(df['beam'].dropna().astype(int).unique().tolist())

        if 2 in beams and 3 in beams:
            return df

        if beams == {2} or beams == {3}:
            source_beam = 2 if 2 in beams else 3
            src = df[df['beam'] == source_beam].copy()

            if len(src) == 0:
                return df

            port = src.copy()
            star = src.copy()

            port['beam'] = 2
            star['beam'] = 3

            half_cnt = (src['ping_cnt'] / 2).astype(int)
            port['ping_cnt'] = half_cnt
            star['ping_cnt'] = half_cnt

            bytes_per_sample = src['bytes_per_sample'].fillna(1).astype(int)
            star['son_offset'] = src['son_offset'] + (half_cnt * bytes_per_sample)

            if 'channel_number' in port.columns:
                port['channel_number'] = 0
            if 'channel_number' in star.columns:
                star['channel_number'] = 1

            combined = pd.concat([port, star], ignore_index=True)
            combined.sort_values(by=['time_s', 'beam'], inplace=True)
            combined.reset_index(drop=True, inplace=True)
            return combined

        return df

    def _sync_port_star_metadata(self, df: pd.DataFrame):
        if 'beam' not in df.columns:
            return df

        if 'freq_band' not in df.columns:
            return self._sync_port_star_metadata_group(df)

        grouped = []
        for _, group in df.groupby('freq_band', dropna=False):
            grouped.append(self._sync_port_star_metadata_group(group))

        out = pd.concat(grouped, ignore_index=True)
        out.sort_values(by=['time_s', 'beam'], inplace=True)
        out.reset_index(drop=True, inplace=True)
        return out

    def _invalid_geom_mask(self, df: pd.DataFrame):
        pix = pd.to_numeric(df['pixM'], errors='coerce') if 'pixM' in df.columns else pd.Series(np.nan, index=df.index)
        spp = pd.to_numeric(df['seconds_per_ping'], errors='coerce') if 'seconds_per_ping' in df.columns else pd.Series(np.nan, index=df.index)
        ping = pd.to_numeric(df['ping_cnt'], errors='coerce') if 'ping_cnt' in df.columns else pd.Series(np.nan, index=df.index)

        return (
            pix.isna() |
            (pix <= 1e-4) |
            spp.isna() |
            (spp <= 1e-6) |
            ping.isna() |
            (ping <= 0) |
            (ping > 1_000_000)
        )

    def _repair_invalid_metadata_rows(self, df: pd.DataFrame):
        if len(df) == 0:
            return df

        out = df.copy()
        invalid = self._invalid_geom_mask(out)
        if not invalid.any():
            return out

        key = 'ping_number' if 'ping_number' in out.columns else 'time_s'
        out.sort_values(by=[key], inplace=True)
        invalid = self._invalid_geom_mask(out)

        repair_fields = ['pixM', 'seconds_per_ping', 'f', 'f_min', 'f_max', 'bytes_per_sample', 'ping_cnt']
        for field in repair_fields:
            if field not in out.columns:
                continue

            series = pd.to_numeric(out[field], errors='coerce')
            series = series.where(~invalid, np.nan)
            series = series.ffill().bfill()

            if field in ['bytes_per_sample', 'ping_cnt']:
                out[field] = np.round(series).astype(int)
            else:
                out[field] = series

        out.reset_index(drop=True, inplace=True)
        return out

    def _sync_port_star_metadata_group(self, df: pd.DataFrame):
        beams = set(df['beam'].dropna().astype(int).unique().tolist())
        if not ({2, 3}.issubset(beams)):
            return self._repair_invalid_metadata_rows(df)

        port = df[df['beam'] == 2].copy()
        star = df[df['beam'] == 3].copy()

        if len(port) == 0 or len(star) == 0:
            return df

        key = 'ping_number' if 'ping_number' in df.columns else 'index'

        port_lookup = port.sort_values(by=[key]).drop_duplicates(subset=[key], keep='first').set_index(key)
        star_lookup = star.sort_values(by=[key]).drop_duplicates(subset=[key], keep='first').set_index(key)

        star_needs_geom = self._invalid_geom_mask(star)
        if 'channel_number' in star.columns:
            star_chan = pd.to_numeric(star['channel_number'], errors='coerce')
            star_needs_geom |= (~star_chan.isin([0, 1, 2, 3, 4, 5]))

        src_ping_for_star = pd.to_numeric(star[key].map(port_lookup['ping_cnt']), errors='coerce')
        dst_ping_star = pd.to_numeric(star['ping_cnt'], errors='coerce')
        ping_mismatch_star = np.isfinite(src_ping_for_star) & np.isfinite(dst_ping_star) & (
            (dst_ping_star > (src_ping_for_star * 2.0)) |
            (dst_ping_star < (src_ping_for_star * 0.5))
        )

        src_spp_for_star = pd.to_numeric(star[key].map(port_lookup['seconds_per_ping']), errors='coerce')
        dst_spp_star = pd.to_numeric(star['seconds_per_ping'], errors='coerce')
        spp_mismatch_star = np.isfinite(src_spp_for_star) & np.isfinite(dst_spp_star) & (
            (dst_spp_star > max(5.0, np.nanmax(src_spp_for_star) * 5.0 if np.isfinite(np.nanmax(src_spp_for_star)) else 5.0)) |
            (dst_spp_star > (src_spp_for_star * 5.0)) |
            (dst_spp_star < (src_spp_for_star * 0.2))
        )

        star_needs_geom |= ping_mismatch_star | spp_mismatch_star

        for idx in star[star_needs_geom].index:
            k = star.loc[idx, key]
            if k not in port_lookup.index:
                continue
            src = port_lookup.loc[k]

            for field in ['pixM', 'seconds_per_ping', 'f', 'f_min', 'f_max', 'bytes_per_sample', 'ping_cnt']:
                if field in star.columns and field in src.index:
                    star.at[idx, field] = src[field]

            if 'channel_number' in star.columns:
                star.at[idx, 'channel_number'] = 1

        port_needs_geom = self._invalid_geom_mask(port)
        if 'channel_number' in port.columns:
            port_chan = pd.to_numeric(port['channel_number'], errors='coerce')
            port_needs_geom |= (~port_chan.isin([0, 1, 2, 3, 4, 5]))

        src_ping_for_port = pd.to_numeric(port[key].map(star_lookup['ping_cnt']), errors='coerce')
        dst_ping_port = pd.to_numeric(port['ping_cnt'], errors='coerce')
        ping_mismatch_port = np.isfinite(src_ping_for_port) & np.isfinite(dst_ping_port) & (
            (dst_ping_port > (src_ping_for_port * 2.0)) |
            (dst_ping_port < (src_ping_for_port * 0.5))
        )

        src_spp_for_port = pd.to_numeric(port[key].map(star_lookup['seconds_per_ping']), errors='coerce')
        dst_spp_port = pd.to_numeric(port['seconds_per_ping'], errors='coerce')
        spp_mismatch_port = np.isfinite(src_spp_for_port) & np.isfinite(dst_spp_port) & (
            (dst_spp_port > max(5.0, np.nanmax(src_spp_for_port) * 5.0 if np.isfinite(np.nanmax(src_spp_for_port)) else 5.0)) |
            (dst_spp_port > (src_spp_for_port * 5.0)) |
            (dst_spp_port < (src_spp_for_port * 0.2))
        )

        port_needs_geom |= ping_mismatch_port | spp_mismatch_port

        for idx in port[port_needs_geom].index:
            k = port.loc[idx, key]
            if k not in star_lookup.index:
                continue
            src = star_lookup.loc[k]
            for field in ['pixM', 'seconds_per_ping', 'f', 'f_min', 'f_max', 'bytes_per_sample', 'ping_cnt']:
                if field in port.columns and field in src.index:
                    port.at[idx, field] = src[field]

            if 'channel_number' in port.columns:
                port.at[idx, 'channel_number'] = 0

        port = self._repair_invalid_metadata_rows(port)
        star = self._repair_invalid_metadata_rows(star)

        out = pd.concat([port, star], ignore_index=True)
        out.sort_values(by=['time_s', 'beam'], inplace=True)
        out.reset_index(drop=True, inplace=True)
        return out

    def _normalize_channel_number(self, channel_number: int, chan_idx: int):
        try:
            channel_number = int(channel_number)
        except Exception:
            channel_number = -1

        if channel_number in (0, 1, 2, 3, 4, 5):
            return channel_number

        # Some Klein files store corrupted/non-canonical channel IDs in packet
        # channel headers while channel order (chan_idx) remains stable.
        if 0 <= int(chan_idx) <= 5:
            return int(chan_idx)

        return None

    def _parse_sonar_record(self, record_start: int, header: bytes, num_chans_to_follow: int, record_bytes: int, file):
        year = struct.unpack_from('<H', header, 14)[0]
        month = header[16]
        day = header[17]
        hour = header[18]
        minute = header[19]
        second = header[20]
        hsecond = header[21]

        ping_number = struct.unpack_from('<I', header, 28)[0]
        sensor_speed_kn = struct.unpack_from('<f', header, 152)[0]
        sensor_y = struct.unpack_from('<d', header, 160)[0]
        sensor_x = struct.unpack_from('<d', header, 168)[0]
        sensor_depth = struct.unpack_from('<f', header, 192)[0]
        sensor_altitude = struct.unpack_from('<f', header, 196)[0]
        sensor_pitch = struct.unpack_from('<f', header, 204)[0]
        sensor_roll = struct.unpack_from('<f', header, 208)[0]
        sensor_heading = struct.unpack_from('<f', header, 212)[0]
        sensor_yaw = struct.unpack_from('<f', header, 220)[0]

        rows = []

        sample_offset = self.ping_header_size + (num_chans_to_follow * self.ping_chan_header_size)

        try:
            dt = datetime(year, max(month, 1), max(day, 1), hour, minute, second)
            time_s = dt.timestamp() + (hsecond / 100.0)
        except Exception:
            time_s = float(ping_number)

        lat, lon, e, n = self._decode_position(sensor_x, sensor_y, self.nav_units)

        ref_ping_cnt = None

        for chan_idx in range(num_chans_to_follow):
            if int(sample_offset) >= int(record_bytes):
                break

            chan_off = self.ping_header_size + (chan_idx * self.ping_chan_header_size)
            file.seek(record_start + chan_off)
            chan_header = file.read(self.ping_chan_header_size)

            if len(chan_header) < self.ping_chan_header_size:
                break

            channel_number = struct.unpack_from('<H', chan_header, 0)[0]
            slant_range = struct.unpack_from('<f', chan_header, 4)[0]
            seconds_per_ping = struct.unpack_from('<f', chan_header, 20)[0]
            frequency = struct.unpack_from('<H', chan_header, 26)[0]
            num_samples = struct.unpack_from('<I', chan_header, 42)[0]
            try:
                num_samples = int(num_samples)
            except Exception:
                num_samples = 0

            chan_cfg = self.chaninfo.get(chan_idx, {})
            type_of_channel = chan_cfg.get('type_of_channel', None)
            bytes_per_sample = max(int(chan_cfg.get('bytes_per_sample', 1)), 1)

            remaining_payload_bytes = max(int(record_bytes) - int(sample_offset), 0)
            max_samples = remaining_payload_bytes // bytes_per_sample

            if type_of_channel is not None and type_of_channel not in [1, 2]:
                skip_samples = max(0, min(int(num_samples), int(max_samples)))
                sample_offset += int(skip_samples) * int(bytes_per_sample)
                continue

            channel_number = self._normalize_channel_number(channel_number, chan_idx)
            if channel_number is None:
                skip_samples = max(0, min(int(num_samples), int(max_samples)))
                sample_offset += int(skip_samples) * int(bytes_per_sample)
                continue

            chan_freq = float(chan_cfg.get('frequency', frequency)) if chan_cfg.get('frequency', frequency) is not None else np.nan

            beam = self._map_beam(channel_number, type_of_channel, chan_idx)
            freq_band = self._map_freq_band(channel_number, chan_idx)

            ping_cnt = int(num_samples)
            if ping_cnt < 0:
                ping_cnt = 0
            if max_samples <= 0:
                ping_cnt = 0
            elif ping_cnt > max_samples:
                ping_cnt = int(max_samples)

            if ref_ping_cnt is None and ping_cnt > 0:
                ref_ping_cnt = int(ping_cnt)

            if (
                ref_ping_cnt is not None and
                ref_ping_cnt > 0 and
                (num_samples <= 0 or int(num_samples) > int(max_samples) or ping_cnt > (ref_ping_cnt * 2))
            ):
                ping_cnt = int(min(ref_ping_cnt, max_samples))

            if ping_cnt <= 0:
                continue

            pix_m = np.nan
            if slant_range > 0:
                pix_m = slant_range / ping_cnt

            sample_format = int(chan_cfg.get('sample_format', 0))
            if bytes_per_sample == 1:
                self.son8bit = True
                self.sample_dtype = '>u1'
            elif bytes_per_sample == 2:
                self.son8bit = False
                self.sample_dtype = '<u2'
            elif bytes_per_sample == 4:
                self.son8bit = False
                # XTF sample format 5 is commonly float32 amplitudes (e.g., Klein).
                if sample_format == 5:
                    self.sample_dtype = '<f4'
                else:
                    self.sample_dtype = '<u4'
            else:
                self.son8bit = False
                self.sample_dtype = '<u2'

            sensor_depth_m = float(sensor_depth) if np.isfinite(sensor_depth) and float(sensor_depth) > 0 else np.nan
            sensor_altitude_m = float(sensor_altitude) if np.isfinite(sensor_altitude) and float(sensor_altitude) > 0 else np.nan

            dep_m = sensor_depth_m
            if not np.isfinite(dep_m) and np.isfinite(sensor_altitude_m):
                dep_m = sensor_altitude_m

            row = {
                'index': int(record_start),
                'son_offset': int(sample_offset),
                'ping_number': int(ping_number),
                'time_s': float(time_s),
                'beam': int(beam),
                'channel_number': int(channel_number),
                'ping_cnt': ping_cnt,
                'bytes_per_sample': int(bytes_per_sample),
                'f': float(chan_freq) if np.isfinite(chan_freq) else np.nan,
                'f_min': float(chan_freq) if np.isfinite(chan_freq) else np.nan,
                'f_max': float(chan_freq) if np.isfinite(chan_freq) else np.nan,
                'freq_band': freq_band,
                'pixM': float(pix_m) if np.isfinite(pix_m) else np.nan,
                'speed_ms': float(sensor_speed_kn) * 0.514444,
                'inst_dep_m': sensor_depth_m,
                'instr_heading': float(sensor_heading),
                'pitch': float(sensor_pitch) if np.isfinite(sensor_pitch) else np.nan,
                'roll': float(sensor_roll) if np.isfinite(sensor_roll) else np.nan,
                'yaw': float(sensor_yaw) if np.isfinite(sensor_yaw) else np.nan,
                'dep_m': dep_m,
                'altitude': sensor_altitude_m,
                'lat': lat,
                'lon': lon,
                'e': e,
                'n': n,
                'seconds_per_ping': float(seconds_per_ping) if np.isfinite(seconds_per_ping) else np.nan,
                'transect': 0,
            }
            rows.append(row)

            sample_offset += int(ping_cnt) * int(bytes_per_sample)

        return rows

    def _decode_position(self, x: float, y: float, nav_units=None):
        if nav_units == 3 and np.isfinite(y) and np.isfinite(x) and abs(y) <= 90 and abs(x) <= 180:
            lat = float(y)
            lon = float(x)
            epsg = self._convert_wgs_to_utm(lon, lat)
            self.humDat['epsg'] = f'EPSG:{epsg}'
            self.humDat['wgs'] = 'EPSG:4326'
            self.trans = pyproj.Proj(self.humDat['epsg'])
            e, n = self.trans(lon, lat)
            return lat, lon, float(e), float(n)

        if nav_units is None and np.isfinite(y) and np.isfinite(x) and abs(y) <= 90 and abs(x) <= 180:
            lat = float(y)
            lon = float(x)
            epsg = self._convert_wgs_to_utm(lon, lat)
            self.humDat['epsg'] = f'EPSG:{epsg}'
            self.humDat['wgs'] = 'EPSG:4326'
            self.trans = pyproj.Proj(self.humDat['epsg'])
            e, n = self.trans(lon, lat)
            return lat, lon, float(e), float(n)

        self.humDat['epsg'] = 'UNKNOWN'
        self.humDat['wgs'] = 'EPSG:4326'
        self.trans = lambda lon, lat: (lon, lat)
        if np.isfinite(x) and np.isfinite(y):
            return np.nan, np.nan, float(x), float(y)
        return np.nan, np.nan, np.nan, np.nan

    def _convert_wgs_to_utm(self, lon: float, lat: float):
        utm_band = str(int((np.floor((lon + 180) / 6) % 60) + 1))
        if len(utm_band) == 1:
            utm_band = '0' + utm_band
        if lat >= 0:
            return '326' + utm_band
        return '327' + utm_band

    def _map_beam(self, channel_number: int, type_of_channel, chan_idx: int):
        if type_of_channel == 1:
            return 2
        if type_of_channel == 2:
            return 3

        if chan_idx % 2 == 0:
            return 2
        return 3

    def _map_freq_band(self, channel_number: int, chan_idx: int = None):
        if channel_number in (0, 1):
            return 'low'
        if channel_number in (2, 3):
            return 'high'
        if channel_number in (4, 5):
            return 'vhigh'

        if chan_idx is not None:
            pair_idx = int(chan_idx) // 2
            if pair_idx == 0:
                return 'low'
            if pair_idx == 1:
                return 'high'
            if pair_idx == 2:
                return 'vhigh'

        return None

    def _doUnitConversion(self, df: pd.DataFrame):
        if 'inst_dep_m' not in df.columns:
            if 'dep_m' in df.columns:
                df['inst_dep_m'] = df['dep_m']
            else:
                df['inst_dep_m'] = np.nan

        if 'dep_m' not in df.columns:
            df['dep_m'] = df['inst_dep_m']

        if 'dep_m' not in df.columns:
            df['dep_m'] = np.nan
        if 'speed_ms' not in df.columns:
            df['speed_ms'] = np.nan
        if 'instr_heading' not in df.columns:
            df['instr_heading'] = np.nan

        if 'f' in df.columns:
            invalid_f = ~np.isfinite(df['f'])
            df.loc[invalid_f, 'f'] = np.nan
            df.loc[invalid_f, 'f_min'] = np.nan
            df.loc[invalid_f, 'f_max'] = np.nan

        df['tempC'] = np.float32(self.tempC * 10)
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

        group_cols = ['beam']
        if 'freq_band' in df.columns:
            group_cols.append('freq_band')

        for group_key, group in df.groupby(group_cols):
            if isinstance(group_key, tuple):
                beam, freq_band = group_key
            else:
                beam, freq_band = group_key, None
            meta = {}

            if beam in [2, 3] and 'pixM' in group.columns and len(group) > 0:
                self.pixM = group['pixM'].iloc[0]

            beam_name = f'B00{int(beam)}'
            base_name = self._getBeamName(beam_name)
            if freq_band:
                meta['beamName'] = f'{base_name}_{freq_band}'
            else:
                meta['beamName'] = base_name
            meta['sonFile'] = self.sonFile

            group = self._getChunkID(group.copy())

            out_csv = f'{beam_name}_{meta["beamName"]}_meta.csv'
            out_csv = os.path.join(self.metaDir, out_csv)
            group.to_csv(out_csv, index=False)

            meta['metaCSV'] = out_csv
            key = beam_name if not freq_band else f'{beam_name}_{freq_band}'
            beamMeta[key] = meta

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
