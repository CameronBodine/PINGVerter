"""Klein SonarPro SDF support for PINGVerter / PINGMapper.

Strict parser for the observed MAX VIEW 600 variant:
1) Parse outer chain to locate inner section start.
2) Detect valid inner pages via config signature.
3) Extract only the primary side-scan stream for PINGMapper.
"""

import os
import struct
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import pyproj


RECORD_HDR = struct.Struct('<HHI')

_OUTER_NAV_PREFIX_BYTES = 40
_FALLBACK_INNER_START = 158024

_INNER_CFG_DL = 1536
_INNER_CFG_FIRST_U32 = 124
_INNER_CFG_SIZE = 8 + _INNER_CFG_DL
_INNER_SONAR_CH = 2
_INNER_SONAR_DL_MIN = 80000
_INNER_SONAR_DL_MAX = 120000
_INNER_MIN_PAGE_GAP = 2000
_INNER_SPEED_OFF = 568
_INNER_LAT_OFF  = 580  # float64 radians; lon follows at +8

_CFG_SIG = struct.pack('<HHI', 0x0000, 0x0000, _INNER_CFG_DL)

_CFG_HEADING_OFF = 76
_CFG_ROLL_OFF = 80
_CFG_PITCH_OFF = 84


class sdf(object):
    """Klein SonarPro SDF file reader compatible with PINGMapper."""

    def __init__(self, inFile: str, nchunk: int = 0, exportUnknown: bool = False):
        self.humFile = None
        self.isOnix = 0
        self.sonFile = inFile
        self.nchunk = nchunk
        self.exportUnknown = exportUnknown

        self.file_header_size = 0
        self.headBytes = 0
        self.humDat = {}
        self.son8bit = False
        self.sample_dtype = '<u2'
        self.flip_port = False

        self.frequency_khz = 150.0
        self.sound_speed_ms = 1500.0
        self._global_timestamp = 0.0
        self._inner_start = _FALLBACK_INNER_START

    def _getFileLen(self):
        self.file_len = os.path.getsize(self.sonFile)

    def _parseFileHeader(self):
        read_len = min(200000, self.file_len)
        with open(self.sonFile, 'rb') as fh:
            raw = fh.read(read_len)

        inner_start = _FALLBACK_INNER_START
        pos = 0
        while pos + 8 <= read_len:
            rt, ch, dl = RECORD_HDR.unpack_from(raw, pos)
            next_pos = pos + 8 + dl

            if rt == 0xFFFF and pos == 0:
                pb = pos + 8
                if pb + 36 <= read_len:
                    freq_hz = struct.unpack_from('<I', raw, pb + 32)[0]
                    if 10000 < freq_hz < 500000:
                        self.frequency_khz = freq_hz / 1000.0
                if pb + 96 <= read_len:
                    try:
                        y, mo, d, h, mi, s = struct.unpack_from('<6I', raw, pb + 64)
                        if 2000 <= y <= 2100 and 1 <= mo <= 12 and 1 <= d <= 31:
                            self._global_timestamp = datetime(y, mo, d, h, mi, s, tzinfo=timezone.utc).timestamp()
                    except Exception:
                        pass

            if next_pos > self.file_len or dl > 5_000_000:
                inner_start = pos + 8 + _OUTER_NAV_PREFIX_BYTES
                break

            pos = next_pos

        self._inner_start = inner_start
        self.file_header_size = inner_start

        self.file_header = {
            'format': 'SDF',
            'frequency_khz': self.frequency_khz,
            'sound_speed_ms': self.sound_speed_ms,
            'inner_start': inner_start,
            'timestamp_utc': self._global_timestamp,
        }

        out_file = os.path.join(self.metaDir, 'DAT_meta.csv')
        pd.DataFrame.from_dict(self.file_header, orient='index').T.to_csv(out_file, index=False)
        self.datMetaFile = out_file

    def _parsePingHeader(self):
        """Parse inner SDF pages extracting only the primary side-scan beam."""
        with open(self.sonFile, 'rb') as fh:
            raw = fh.read()

        inner_start = int(getattr(self, '_inner_start', _FALLBACK_INNER_START))
        pages = self._find_inner_pages(raw, inner_start)
        if not pages:
            raise ValueError(f'No strict SDF inner pages found from offset {inner_start}.')

        conv_path = os.path.join(self.metaDir, 'sdf_converted.bin')
        rows = []
        write_off = 0

        with open(conv_path, 'wb') as out_fh:
            for ping_idx, page_pos in enumerate(pages):
                cfg_payload_start = page_pos + 8

                nav = {'heading': 0.0, 'roll': 0.0, 'pitch': 0.0}
                nav = self._parse_config_nav(raw, cfg_payload_start, _INNER_CFG_DL, nav)
                speed_ms = self._parse_config_speed(raw, cfg_payload_start, _INNER_CFG_DL)
                lat_d, lon_d = self._parse_config_position(raw, cfg_payload_start, _INNER_CFG_DL)

                sonar_records = self._extract_page_sonar_records(raw, page_pos)
                if not sonar_records:
                    continue

                if self._global_timestamp > 0:
                    time_s = float(self._global_timestamp + ping_idx)
                else:
                    time_s = float(ping_idx)

                for beam_id, samps, rt in sonar_records:
                    if len(samps) < 32:
                        continue

                    out_fh.write(samps.tobytes())

                    rows.append({
                        'index': int(write_off),
                        'son_offset': 0,
                        'bytes_per_sample': 2,
                        'beam': int(beam_id),
                        'chan_id': int(beam_id),
                        'ping_number': int(ping_idx),
                        'time_s': time_s,
                        'ping_cnt': int(len(samps)),
                        'f': float(self.frequency_khz),
                        'f_min': float(self.frequency_khz),
                        'f_max': float(self.frequency_khz),
                        'pixM': np.nan,
                        'speed_ms': float(speed_ms),
                        'instr_heading': float(nav.get('heading', 0.0)),
                        'pitch': float(nav.get('pitch', 0.0)),
                        'roll': float(nav.get('roll', 0.0)),
                        'yaw': np.nan,
                        'inst_dep_m': np.nan,
                        'dep_m': np.nan,
                        'altitude': np.nan,
                        'lat': lat_d,
                        'lon': lon_d,
                        'e': np.nan,
                        'n': np.nan,
                        'layback_m': 0.0,
                        'flags': int(rt),
                        'transect': 0,
                    })
                    write_off += int(len(samps) * 2)

        if not rows:
            raise ValueError('SDF strict parser found pages, but no valid sonar samples were decoded.')

        self.sonFile = conv_path

        df = pd.DataFrame(rows)

        # Compute UTM e/n from decoded GPS lat/lon
        valid_mask = df['lat'].notna() & df['lon'].notna()
        if valid_mask.any():
            first_valid = df[valid_mask].iloc[0]
            epsg = self._convert_wgs_to_utm(first_valid['lon'], first_valid['lat'])
            self.humDat['epsg'] = f'EPSG:{epsg}'
            self.humDat['wgs'] = 'EPSG:4326'
            proj = pyproj.Proj(self.humDat['epsg'])
            es, ns = proj(df.loc[valid_mask, 'lon'].to_numpy(),
                          df.loc[valid_mask, 'lat'].to_numpy())
            df.loc[valid_mask, 'e'] = es
            df.loc[valid_mask, 'n'] = ns
        else:
            self.humDat['epsg'] = 'UNKNOWN'
            self.humDat['wgs'] = 'EPSG:4326'

        df = self._doUnitConversion(df)
        df.sort_values(by=['time_s', 'beam'], inplace=True)
        df.reset_index(drop=True, inplace=True)
        df = self._calcTrkDistTS(df)
        df['record_num'] = np.arange(len(df), dtype=np.int64)

        out_all = os.path.join(self.metaDir, 'All-SDF-Sonar-MetaData.csv')
        df.to_csv(out_all, index=False)
        self.header_dat = df

    def _extract_page_sonar_records(self, raw: bytes, page_pos: int):
        """Extract primary side-scan record and split dual-channel payload into
        port (beam 2) and starboard (beam 3).

        The main sonar record (ch=2, dl≈98 KB) contains sub-channels separated
        by u32 words whose upper 16 bits are small integers (≤30) acting as
        channel markers.  CH2 — the largest segment — is a dual swath with the
        port half stored near→far, a nadir blank, then the starboard half also
        near→far.  Both halves are returned as separate beams.
        """
        scan_pos = page_pos + _INNER_CFG_SIZE
        page_end = min(scan_pos + 200000, len(raw))

        while scan_pos + 8 < page_end:
            try:
                rt, ch, dl = RECORD_HDR.unpack_from(raw, scan_pos)
            except Exception:
                break

            if rt == 0x0000 and ch == 0 and dl == _INNER_CFG_DL:
                break

            if dl < 1000 or dl > 300000:
                scan_pos += 8 + max(dl, 0)
                continue

            data_end = scan_pos + 8 + dl
            if data_end > len(raw):
                break

            if ch == _INNER_SONAR_CH and _INNER_SONAR_DL_MIN <= dl <= _INNER_SONAR_DL_MAX:
                try:
                    u32 = np.frombuffer(raw[scan_pos + 8:data_end], dtype='<u4')
                    ch2_data = self._extract_ch2(u32)
                except Exception:
                    scan_pos += 8 + dl
                    continue

                if ch2_data is not None and len(ch2_data) > 200:
                    port, star = self._split_port_star(ch2_data)
                    results = []
                    if len(port) > 32:
                        results.append((2, port, rt))
                    if len(star) > 32:
                        results.append((3, star, rt))
                    if results:
                        return results

            scan_pos += 8 + dl

        return []

    def _extract_ch2(self, u32: np.ndarray) -> np.ndarray:
        """Return the lower-16-bit sample values for the dual-channel (CH2) segment.

        CH2 is the segment between the first and second true channel-marker
        groups.  True markers have small upper-16-bit IDs (≤30); the large
        nav-data packet embedded later in the payload has IEEE-float upper bytes
        (≥100) and is excluded from channel detection.
        """
        upper16 = (u32 >> 16).astype(np.uint32)
        marker_pos = np.where(upper16 > 0)[0]

        if len(marker_pos) == 0:
            # No sub-channel structure detected — return all samples as-is
            return (u32 & 0xFFFF).astype(np.float32)

        # Collect consecutive marker runs; keep only those where every marker
        # in the run has a small upper-16 value (true channel separator, not nav data).
        groups = []
        i = 0
        while i < len(marker_pos):
            j = i
            while j + 1 < len(marker_pos) and marker_pos[j + 1] - marker_pos[j] <= 2:
                j += 1
            grp_upper = upper16[marker_pos[i:j + 1]]
            if np.all(grp_upper <= 30):
                groups.append((int(marker_pos[i]), int(marker_pos[j])))
            i = j + 1

        if len(groups) < 1:
            return (u32 & 0xFFFF).astype(np.float32)

        # CH2 occupies the range between the end of group[0] and the start of group[1]
        ch2_start = groups[0][1] + 1
        ch2_end = groups[1][0] if len(groups) > 1 else len(u32)
        return (u32[ch2_start:ch2_end] & 0xFFFF).astype(np.float32)

    def _split_port_star(self, ch2_data: np.ndarray, window: int = 100,
                         blank_thresh: float = 150.0,
                         star_wc_thresh: float = 3000.0):
        """Split a dual near→far swath into port and starboard halves.

        Returns (port_u16, star_u16) where each array is uint16, near→far.

        Strategy
        --------
        1. Locate blank_start: first sample where the rolling mean falls below
           *blank_thresh* AFTER the port returns have peaked (roll > 2000 at
           some prior point).  Using 2000 avoids false triggers on the weak
           secondary echo at ~40% into CH2.
        2. Locate star_start: walk back from the star water-column peak
           (roll > *star_wc_thresh*, found in the second portion) to the last
           sample still below star_wc_thresh/4 — that is the pre-star silence.
        """
        n = len(ch2_data)
        roll = np.convolve(ch2_data.astype(np.float64), np.ones(window) / window,
                           mode='valid')
        nr = len(roll)

        # ── Step 1: blank_start ──────────────────────────────────────────────
        high_idx = np.where(roll[:nr // 2] > 2000.0)[0]
        search_after = int(high_idx[-1]) + window if len(high_idx) > 0 else nr // 4

        low_cands = np.where(roll[search_after: search_after + nr // 2] < blank_thresh)[0]
        blank_start = (search_after + int(low_cands[0])) if len(low_cands) > 0 else nr // 3

        # ── Step 2: star_start ───────────────────────────────────────────────
        # Look for star water-column peak at least 200 samples beyond blank_start
        # but skip over the first quarter of the remaining range to avoid echoes.
        remaining = nr - blank_start
        search_for_star = blank_start + max(200, remaining // 4)

        high_star = np.where(roll[search_for_star:] > star_wc_thresh)[0]
        if len(high_star) > 0:
            star_wc_idx = search_for_star + int(high_star[0])
            # Walk backwards: find the last quiet sample (< star_wc_thresh/4)
            # between blank_start and star_wc_idx — that marks where star begins.
            before_peak = roll[blank_start: star_wc_idx]
            quiet = np.where(before_peak < star_wc_thresh / 4.0)[0]
            if len(quiet) > 0:
                star_start = blank_start + int(quiet[-1]) + window
            else:
                star_start = star_wc_idx
        else:
            # No clear star peak found; fall back to symmetric split
            star_start = blank_start + remaining // 2

        port = ch2_data[:blank_start].astype(np.uint16)
        star = ch2_data[max(0, star_start):].astype(np.uint16)
        return port, star

    def _find_inner_pages(self, raw: bytes, inner_start: int):
        pages = []
        search_pos = max(0, int(inner_start))
        last_page = -10**12
        raw_len = len(raw)

        while True:
            pos = raw.find(_CFG_SIG, search_pos)
            if pos < 0:
                break
            search_pos = pos + 1

            if pos - last_page < _INNER_MIN_PAGE_GAP:
                continue
            if pos + _INNER_CFG_SIZE + 8 > raw_len:
                continue

            first_u32 = struct.unpack_from('<I', raw, pos + 8)[0]
            if first_u32 != _INNER_CFG_FIRST_U32:
                continue

            sonar_hdr_pos = pos + _INNER_CFG_SIZE
            rt2, ch2, dl2 = RECORD_HDR.unpack_from(raw, sonar_hdr_pos)
            if ch2 != _INNER_SONAR_CH:
                continue
            if not (_INNER_SONAR_DL_MIN <= dl2 <= _INNER_SONAR_DL_MAX):
                continue
            if sonar_hdr_pos + 8 + dl2 > raw_len:
                continue

            pages.append(pos)
            last_page = pos

        return pages

    def _parse_config_position(self, raw: bytes, payload_start: int, payload_len: int):
        """Decode lat/lon from config payload stored as float64 radians."""
        off = payload_start + _INNER_LAT_OFF
        if off + 16 > len(raw):
            return np.nan, np.nan
        try:
            lat_r, lon_r = struct.unpack_from('<dd', raw, off)
            if not (np.isfinite(lat_r) and np.isfinite(lon_r)):
                return np.nan, np.nan
            lat_d = float(lat_r * 180.0 / np.pi)
            lon_d = float(lon_r * 180.0 / np.pi)
            if abs(lat_d) > 90 or abs(lon_d) > 180 or (lat_d == 0.0 and lon_d == 0.0):
                return np.nan, np.nan
            return lat_d, lon_d
        except Exception:
            return np.nan, np.nan

    def _parse_config_speed(self, raw: bytes, payload_start: int, payload_len: int):
        off = payload_start + _INNER_SPEED_OFF
        if off + 4 > len(raw):
            return 0.0
        try:
            v = struct.unpack_from('<f', raw, off)[0]
        except Exception:
            return 0.0

        if np.isfinite(v) and 0.0 <= v <= 20.0:
            return float(v)
        return 0.0

    def _parse_config_nav(self, raw: bytes, payload_start: int, payload_len: int, prev_nav: dict) -> dict:
        nav = dict(prev_nav)
        try:
            if payload_start + _CFG_PITCH_OFF + 4 <= len(raw):
                h = struct.unpack_from('<f', raw, payload_start + _CFG_HEADING_OFF)[0]
                r = struct.unpack_from('<f', raw, payload_start + _CFG_ROLL_OFF)[0]
                p = struct.unpack_from('<f', raw, payload_start + _CFG_PITCH_OFF)[0]
                if np.isfinite(h) and -720.0 <= h <= 720.0:
                    nav['heading'] = float(h) % 360.0
                if np.isfinite(r) and -180.0 <= r <= 180.0:
                    nav['roll'] = float(r)
                if np.isfinite(p) and -90.0 <= p <= 90.0:
                    nav['pitch'] = float(p)
        except Exception:
            pass
        return nav

    def _recalcRecordNum(self):
        df = self.header_dat.reset_index(drop=True)
        df['record_num'] = df.index
        self.header_dat = df

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

    def _getBeamName(self, beam: str):
        if beam == 'B002':
            return 'ss_port'
        if beam == 'B003':
            return 'ss_star'
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

    def _doUnitConversion(self, df: pd.DataFrame):
        if 'inst_dep_m' not in df.columns:
            if 'dep_m' in df.columns:
                df['inst_dep_m'] = df['dep_m']
            else:
                df['inst_dep_m'] = np.nan

        if 'dep_m' not in df.columns:
            df['dep_m'] = df['inst_dep_m']

        if 'altitude' in df.columns:
            alt = pd.to_numeric(df['altitude'], errors='coerce')
            inst = pd.to_numeric(df['inst_dep_m'], errors='coerce')
            dep = pd.to_numeric(df['dep_m'], errors='coerce')

            inst_fb = (~np.isfinite(inst)) & np.isfinite(alt)
            dep_fb = (~np.isfinite(dep)) & np.isfinite(alt)

            if inst_fb.any():
                df.loc[inst_fb, 'inst_dep_m'] = alt.loc[inst_fb]
            if dep_fb.any():
                df.loc[dep_fb, 'dep_m'] = alt.loc[dep_fb]

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

    def _safe_timestamp(self, timestamp: float):
        if not np.isfinite(timestamp):
            return 0.0
        ts = float(timestamp)
        if ts <= 0 or ts > 4102444800.0:
            return 0.0
        try:
            _ = datetime.fromtimestamp(ts, tz=timezone.utc)
            return ts
        except Exception:
            return 0.0

    def _decode_position(self, latitude: float, longitude: float):
        if (np.isfinite(latitude) and np.isfinite(longitude)
                and abs(latitude) <= 90 and abs(longitude) <= 180
                and not (latitude == 0.0 and longitude == 0.0)):
            lat = float(latitude)
            lon = float(longitude)
            epsg = self._convert_wgs_to_utm(lon, lat)
            self.humDat['epsg'] = f'EPSG:{epsg}'
            self.humDat['wgs'] = 'EPSG:4326'
            self.trans = pyproj.Proj(self.humDat['epsg'])
            e, n = self.trans(lon, lat)
            return lat, lon, float(e), float(n)

        self.humDat['epsg'] = 'UNKNOWN'
        self.humDat['wgs'] = 'EPSG:4326'
        self.trans = lambda lon, lat: (lon, lat)
        return np.nan, np.nan, np.nan, np.nan

    def _convert_wgs_to_utm(self, lon: float, lat: float):
        utm_band = str(int((np.floor((lon + 180) / 6) % 60) + 1))
        if len(utm_band) == 1:
            utm_band = '0' + utm_band
        if lat >= 0:
            return '326' + utm_band
        return '327' + utm_band
