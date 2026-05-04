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
import os, sys
from pingverter import hum, low, cerul, gar, jsf, xtf
import time
import numpy as np
import pandas as pd
from glob import glob

from joblib import Parallel, delayed, cpu_count


def _identity_transform(lon, lat):
    return lon, lat

# =========================================================
# Humminbird to PINGMapper
# =========================================================

def hum2pingmapper(input: str, out_dir: str, nchunk: int=500, tempC: float=10, exportUnknown: bool=False):
    '''
    
    '''
    # Make sure input exists
    assert os.path.isfile(input), "{} does not exist.".format(input)

    # Create the class
    humminbird = hum(humFile=input, nchunk=nchunk, exportUnknown=exportUnknown)

    # Store temperatue
    humminbird.tempC = float(tempC)/10

    #################
    # Decode DAT File
    #################

    start_time = time.time()
    print("\nGetting DAT Metadata...")
    print(input)
    humminbird._getHumDatStruct()

    # Read in the humdat data
    if humminbird.isOnix == 0:
        humminbird._getHumdat()
    else:
        humminbird._decodeOnix()

    # Create 'meta' directory if it doesn't exist
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    metaDir = os.path.join(out_dir, 'meta')
    try:
        os.mkdir(metaDir)
    except:
        pass
    humminbird.metaDir = metaDir #Store metadata directory in sonObj

    print("\nDone!")
    print("Time (s):", round(time.time() - start_time, ndigits=1))

    #####################################
    # Generate sonar file meta/attributes
    #####################################

    # Find available SON files
    sonFiles = glob(os.path.join(humminbird.sonFile, '*.SON'))
    
    # Determine which sonar beams are present (B000,B001,..)
    chanAvail = {}
    for s in sonFiles:
        beam = os.path.split(s)[-1].split('.')[0] #Get beam number (B000,B001,..)
        chanAvail[beam] = s

    humminbird.chanAvail = chanAvail

    # Store sonar file meta in humminbird object
    ## Dictionary of dictionaries
    beamMeta = humminbird.beamMeta = {}
    for chan, file in chanAvail.items():
        
        beamMeta[chan] = {}

        # Get beam name
        beamMeta[chan]['beamName'] = humminbird._getBeamName(chan)

        beamMeta[chan]['beam'] = chan
        beamMeta[chan]['sonFile'] = file

        # Output csv name
        csv = '{}_{}_{}'.format(chan, beamMeta[chan]['beamName'], 'meta.csv')
        beamMeta[chan]['metaCSV'] = os.path.join(metaDir, csv)


    ##############################################################
    # Determine ping header structure (varies by Humminbird model)
    ##############################################################
    start_time = time.time()
    print("\nGetting Header Structure...")
    

    gotHeader = False # Flag indicating if length of header is found

    for beam, meta in beamMeta.items():
        
        # Get SON file
        son = meta['sonFile']

        # Count headbytes
        headbytes = humminbird._cntHead(son)

        if headbytes > 0: # Header length found
            print("Header Length: {}".format(headbytes))
            gotHeader = True

            # Add headbytes to humminbird object
            humminbird.frame_header_size = headbytes

            auto_decode = False

            break

        else:
            auto_decode = True

            gotHeader = False

    
    # Consider adding decodeHeadStruct() function back into PINGVerter....
    if not gotHeader:
        # sys.exit("\n#####\nERROR: Out of SON files... \n"+
        #         "Unable to determine header length.")
        print("\n#####\nERROR: Out of SON files... \n\n"+
                "Trying to automatically determine header length...")
        for beam, meta in beamMeta.items():
        
            # Get SON file
            son = meta['sonFile']

            # Autodecode
            headbytes = humminbird._decodeHeadStruct(son)

            if headbytes > 0:
                print("\n######################\nSLAMMA-JAMMA-DING-DONG \n"+
                      "Header Length Determined: {}\n\n".format(headbytes))
                print("As you were....\n\n")
                gotHeader = True

                # Add headbytes to humminbird object
                humminbird.frame_header_size = headbytes
                break
    if not gotHeader:
        sys.exit("\n#####\nERROR: Out of SON files... \n"+
                "Unable to automatically decode sonar header.")
        
    
    #############################################
    # Get the SON header structure and attributes
    #############################################
    if not auto_decode:
        humminbird._getHeadStruct()

    ##################
    # Parse son header
    ##################

    r = Parallel(n_jobs = len(beamMeta), verbose=10 )(delayed(humminbird._parsePingHeader)(meta['sonFile'], meta['metaCSV']) for beam, meta in beamMeta.items())

    # Store spatial transformation
    for (trans, humdat) in r:
        humminbird.trans = trans
        humminbird.humDat = humdat
        break

    # Save DAT metadata to file (csv)
    outFile = os.path.join(metaDir, 'DAT_meta.csv') # Specify file directory & name
    pd.DataFrame.from_dict(humminbird.humDat, orient='index').T.to_csv(outFile, index=False) # Export DAT df to csv
    humminbird.datMetaFile = outFile # Store metadata file path in sonObj
    del outFile

    print("\nDone!")
    print("Time (s):", round(time.time() - start_time, ndigits=1))
    
    return humminbird


# =========================================================
# Lowrance to PINGMapper
# =========================================================

def low2pingmapper(input: str, out_dir: str, nchunk: int=500, tempC: float=10, exportUnknown: bool=False):

    # Make sure input exists
    assert os.path.isfile(input), "{} does not exist.".format(input)

    # Create the class
    lowrance = low(inFile=input, nchunk=nchunk, exportUnknown=exportUnknown)

    # Store temperature
    lowrance.tempC = float(tempC)/10

    ######################
    # Decode Lowrance File
    ######################

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Create 'meta' directory if it doesn't exist
    metaDir = os.path.join(out_dir, 'meta')
    try:
        os.mkdir(metaDir)
    except:
        pass
    lowrance.metaDir = metaDir # Store metadata directory

    # Get Lowrance file length
    lowrance._getFileLen()

    # Parse file header ***Probably not needed***
    lowrance._parseFileHeader()

    # Parse ping headers (attributes) and do conversions
    lowrance._parsePingHeader()

    # Remove unknown beams
    lowrance._removeUnknownBeams()

    # Drop Beam 0 (83kHz) or 1 (200kHz) if necessary
    lowrance._removeDownBeams()

    # Split sidescan, if necessary
    beams = lowrance.header_dat['beam'].unique()
    if 5 in beams:
        lowrance._splitLowSS()
        flip_port = True
    else:
        flip_port = False

    # Recalculate record number
    lowrance._recalcRecordNum()

    # Drop unknown
    if not exportUnknown:
        cols = lowrance.header_dat.columns
        cols = [c for c in cols if 'unknown' in c]
        
        lowrance.header_dat.drop(columns=cols, inplace=True)

    # Save ping metadata to csv based on beam
    lowrance._splitBeamsToCSV()

    # Store headBytes
    lowrance.headBytes = lowrance.frame_header_size

    # Not Humminbird Onix
    lowrance.isOnix = 0
    
    return lowrance

# =========================================================
# Garmin to PINGMapper
# =========================================================

def gar2pingmapper(input: str, out_dir: str, nchunk: int=500, tempC: float=10, exportUnknown: bool=False):

    # Make sure input exists
    assert os.path.isfile(input), "{} does not exist.".format(input)

    # Create the class
    garmin = gar(inFile=input, nchunk=nchunk, exportUnknown=exportUnknown)
    
    # Store temperature
    garmin.tempC = float(tempC)/10

    ######################
    # Decode Lowrance File
    ######################

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Create 'meta' directory if it doesn't exist
    metaDir = os.path.join(out_dir, 'meta')
    try:
        os.mkdir(metaDir)
    except:
        pass
    garmin.metaDir = metaDir # Store metadata directory

    # Get Garmin file length
    garmin._getFileLen()

    # Parse file header
    garmin._parseFileHeader()

    # Create 'meta' directory if it doesn't exist
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    metaDir = os.path.join(out_dir, 'meta')
    try:
        os.mkdir(metaDir)
    except:
        pass
    garmin.metaDir = metaDir #Store metadata directory in sonObj

    # Save DAT metadata to file (csv)
    outFile = os.path.join(metaDir, 'DAT_meta.csv') # Specify file directory & name
    pd.DataFrame.from_dict(garmin.file_header, orient='index').T.to_csv(outFile, index=False) # Export DAT df to csv
    garmin.datMetaFile = outFile # Store metadata file path in sonObj
    del outFile

    # Parse ping headers (attributes) and do conversions
    garmin._parsePingHeader()

    # Drop unknown
    if not exportUnknown:
        cols = garmin.header_dat.columns
        cols = [c for c in cols if 'unknown' in c]

        garmin.header_dat.drop(columns=cols, inplace=True)

    # Recalculate record num
    garmin._recalcRecordNum()

    # Split and re-label beams to PING-Mapper convention
    garmin._splitBeamsToCSV()

    # Not Humminbird Onix
    garmin.isOnix = 0
    
    return garmin


    


# =========================================================
# Cerulean to PINGMapper
# =========================================================

def cerul2pingmapper(input: str, out_dir: str, nchunk: int=500, tempC: float=10, exportUnknown: bool=False):
    '''
    '''
    # Make sure input exists
    assert os.path.isfile(input), "{} does not exist.".format(input)

    # Create the class
    cerulean = cerul(svlog = input, nchunk=nchunk, exportUnknown=exportUnknown)

    # Store Temperature
    cerulean.tempC = float(tempC)/10

    ######################
    # Decode Cerulean File
    ######################

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Create 'meta' directory if it doesn't exist
    metaDir = os.path.join(out_dir, 'meta')
    try:
        os.mkdir(metaDir)
    except:
        pass
    cerulean.metaDir = metaDir # Store metadata directory

    # Get Cerulean file length
    cerulean._getFileLen()

    # Parse the file header
    cerulean._parseFileHeader()

    # Parse all packet headers
    if exportUnknown:
        cerulean._locatePacketsRaw()

    # Locate Packet Headers
    cerulean._locatePackets()

    # Set beam
    cerulean._convertBeam()

    # Set frequency
    cerulean._convertFrequency()
    
    # Recalculate record num
    cerulean._recalcRecordNum()

    # Save to file
    cerulean._splitBeamsToCSV()

    # print(cerulean)

    return cerulean


# =========================================================
# JSF to PINGMapper
# =========================================================

def jsf2pingmapper(input: str, out_dir: str, nchunk: int=500, tempC: float=10, exportUnknown: bool=False):
    assert os.path.isfile(input), "{} does not exist.".format(input)

    jsf_obj = jsf(inFile=input, nchunk=nchunk, exportUnknown=exportUnknown)
    jsf_obj.tempC = float(tempC)/10

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    metaDir = os.path.join(out_dir, 'meta')
    try:
        os.mkdir(metaDir)
    except:
        pass
    jsf_obj.metaDir = metaDir

    jsf_obj._getFileLen()
    jsf_obj._parseFileHeader()
    jsf_obj._parsePingHeader()
    jsf_obj._recalcRecordNum()
    jsf_obj._splitBeamsToCSV()

    if not hasattr(jsf_obj, 'trans'):
        jsf_obj.trans = lambda lon, lat: (lon, lat)

    return jsf_obj


# =========================================================
# XTF to PINGMapper
# =========================================================

def xtf2pingmapper(input: str, out_dir: str, nchunk: int=500, tempC: float=10, exportUnknown: bool=False):
    assert os.path.isfile(input), "{} does not exist.".format(input)

    xtf_obj = xtf(inFile=input, nchunk=nchunk, exportUnknown=exportUnknown)
    xtf_obj.tempC = float(tempC)/10

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    metaDir = os.path.join(out_dir, 'meta')
    try:
        os.mkdir(metaDir)
    except:
        pass
    xtf_obj.metaDir = metaDir

    xtf_obj._getFileLen()
    xtf_obj._parseFileHeader()
    xtf_obj._parsePingHeader()
    xtf_obj._recalcRecordNum()
    xtf_obj._splitBeamsToCSV()

    if not hasattr(xtf_obj, 'trans'):
        xtf_obj.trans = lambda lon, lat: (lon, lat)

    return xtf_obj


def _row_value(row, names, default=np.nan):
    for name in names:
        if name in row.index:
            value = row[name]
            if pd.notna(value):
                return value
    return default


def _safe_int(value, default=0):
    try:
        if pd.isna(value):
            return default
        return int(value)
    except Exception:
        return default


def _safe_float(value, default=np.nan):
    try:
        if pd.isna(value):
            return default
        return float(value)
    except Exception:
        return default


def _none_if_nan(value):
    if pd.isna(value):
        return None
    try:
        f = float(value)
    except Exception:
        return None
    return None if not np.isfinite(f) else f


def _detect_format(input_path: str, source_format: str=None):
    ext = (source_format or os.path.splitext(input_path)[1]).strip().lower()
    if ext not in SUPPORTED_SONAR_EXTENSIONS:
        raise ValueError(
            "Unsupported recording extension {}. Supported: {}".format(
                ext,
                ', '.join(SUPPORTED_SONAR_EXTENSIONS)
            )
        )
    return ext


def _build_sonar_object(input_path: str, work_dir: str, ext: str, nchunk: int, tempC: float, exportUnknown: bool):
    if ext == '.rsd':
        return gar2pingmapper(input_path, work_dir, nchunk=nchunk, tempC=tempC, exportUnknown=exportUnknown)
    if ext in ('.sl2', '.sl3'):
        return low2pingmapper(input_path, work_dir, nchunk=nchunk, tempC=tempC, exportUnknown=exportUnknown)
    if ext == '.svlog':
        return cerul2pingmapper(input_path, work_dir, nchunk=nchunk, tempC=tempC, exportUnknown=exportUnknown)
    if ext == '.jsf':
        return jsf2pingmapper(input_path, work_dir, nchunk=nchunk, tempC=tempC, exportUnknown=exportUnknown)
    if ext == '.xtf':
        return xtf2pingmapper(input_path, work_dir, nchunk=nchunk, tempC=tempC, exportUnknown=exportUnknown)
    if ext == '.dat':
        return hum2pingmapper(input_path, work_dir, nchunk=nchunk, tempC=tempC, exportUnknown=exportUnknown)
    raise ValueError("Unsupported recording extension {}".format(ext))


def _frame_key_series(df: pd.DataFrame):
    if 'sequence_cnt' in df.columns:
        return pd.to_numeric(df['sequence_cnt'], errors='coerce')
    if 'ping_number' in df.columns:
        return pd.to_numeric(df['ping_number'], errors='coerce')
    if 'id' in df.columns:
        return pd.to_numeric(df['id'], errors='coerce')
    if 'time_s' in df.columns:
        return pd.to_numeric(df['time_s'], errors='coerce').round(3)
    if 'record_num' in df.columns:
        return pd.to_numeric(df['record_num'], errors='coerce')
    return pd.Series(np.arange(len(df)), index=df.index, dtype='float64')


def _humminbird_sequence_from_time(df: pd.DataFrame):
    times = pd.to_numeric(df.get('time_s'), errors='coerce')
    if times is None or times.notna().sum() == 0:
        return pd.Series(np.arange(len(df)), index=df.index, dtype='int64')

    # Estimate ping period from within-channel deltas when available.
    period_candidates = []
    if 'channel_id' in df.columns:
        for _, g in df.groupby('channel_id'):
            t = pd.to_numeric(g.get('time_s'), errors='coerce').dropna().sort_values()
            if len(t) < 2:
                continue
            dt = np.diff(t.to_numpy())
            dt = dt[dt > 0]
            if len(dt) > 0:
                period_candidates.append(float(np.median(dt)))

    if period_candidates:
        period = float(np.median(period_candidates))
    else:
        all_t = np.sort(times.dropna().to_numpy())
        if len(all_t) < 2:
            return pd.Series(np.arange(len(df)), index=df.index, dtype='int64')
        dt_all = np.diff(all_t)
        dt_all = dt_all[dt_all > 0]
        if len(dt_all) == 0:
            return pd.Series(np.arange(len(df)), index=df.index, dtype='int64')
        period = float(np.percentile(dt_all, 75))

    # Keep threshold above timestamp jitter but well below ping period.
    threshold = max(0.003, period * 0.2)

    # Stable sort by time, then by record_num when present.
    work = pd.DataFrame(index=df.index)
    work['_time'] = times
    if 'record_num' in df.columns:
        work['_rn'] = pd.to_numeric(df['record_num'], errors='coerce')
    else:
        work['_rn'] = np.arange(len(df))
    order = work.sort_values(['_time', '_rn'], na_position='last', kind='mergesort').index

    seq = pd.Series(index=df.index, dtype='float64')
    current_seq = -1
    prev_time = None
    for idx in order:
        t = times.loc[idx]
        if not np.isfinite(t):
            continue
        if prev_time is None or (t - prev_time) > threshold:
            current_seq += 1
        seq.loc[idx] = current_seq
        prev_time = t

    # Fill any NaN sequence ids deterministically.
    if seq.notna().any():
        fill_value = int(seq.dropna().max()) + 1
    else:
        fill_value = 0
    seq = seq.fillna(fill_value).astype('int64')
    return seq


def _bytes_per_sample(sonar_obj, row):
    bps = _safe_int(_row_value(row, ['bytes_per_sample']), 0)
    if bps > 0:
        return bps

    if 'sample_dtype' in sonar_obj.__dict__:
        dtype = str(getattr(sonar_obj, 'sample_dtype', ''))
        if dtype.endswith('u1') or dtype.endswith('i1'):
            return 1
        if dtype.endswith('u2') or dtype.endswith('i2'):
            return 2
        if dtype.endswith('u4') or dtype.endswith('i4') or dtype.endswith('f4'):
            return 4

    son8 = getattr(sonar_obj, 'son8bit', False)
    return 1 if son8 else 2


def _normalize_to_u16(arr):
    if arr.size == 0:
        return np.array([], dtype=np.uint16)

    if arr.dtype == np.uint16:
        return arr

    values = np.asarray(arr)
    if np.issubdtype(values.dtype, np.integer):
        clipped = np.clip(values.astype(np.int64), 0, 65535)
        return clipped.astype(np.uint16)

    vals = values.astype(np.float64)
    finite = np.isfinite(vals)
    if not finite.any():
        return np.zeros(vals.shape, dtype=np.uint16)

    out = np.zeros(vals.shape, dtype=np.uint16)
    valid = vals[finite]
    vmin = np.nanpercentile(valid, 2)
    vmax = np.nanpercentile(valid, 98)
    if not np.isfinite(vmin):
        vmin = np.nanmin(valid)
    if not np.isfinite(vmax):
        vmax = np.nanmax(valid)

    if vmax <= vmin:
        out[finite] = 0
    else:
        scaled = (valid - vmin) / (vmax - vmin)
        scaled = np.clip(scaled, 0.0, 1.0)
        out_vals = np.rint(scaled * 65535.0).astype(np.uint16)
        out[finite] = out_vals

    return out


def _decode_raw_to_u16(raw: bytes, bytes_per_sample: int, sonar_obj):
    if bytes_per_sample == 1:
        arr8 = np.frombuffer(raw, dtype=np.uint8)
        return (arr8.astype(np.uint16) * 257)

    sample_dtype = str(getattr(sonar_obj, 'sample_dtype', ''))

    if bytes_per_sample == 2:
        dtype = np.dtype(sample_dtype) if sample_dtype in ('<u2', '>u2', '<i2', '>i2') else np.dtype('<u2')
        arr = np.frombuffer(raw, dtype=dtype)
        return _normalize_to_u16(arr)

    if bytes_per_sample == 4:
        if sample_dtype in ('<f4', '>f4'):
            arr = np.frombuffer(raw, dtype=np.dtype(sample_dtype))
            return _normalize_to_u16(arr)
        dtype = np.dtype(sample_dtype) if sample_dtype in ('<u4', '>u4', '<i4', '>i4') else np.dtype('<u4')
        arr = np.frombuffer(raw, dtype=dtype)
        return _normalize_to_u16(arr)

    return np.array([], dtype=np.uint16)


def _channel_metadata(group: pd.DataFrame, channel_id: int):
    beam = None
    if 'beam' in group.columns and len(group['beam'].dropna()) > 0:
        beam = _safe_int(group['beam'].dropna().mode().iloc[0], 0)

    start_freq = _safe_int(_row_value(group.iloc[0], ['f_min', 'start_freq_hz', 'startFrequencyHz']), 0)
    end_freq = _safe_int(_row_value(group.iloc[0], ['f_max', 'end_freq_hz', 'endFrequencyHz']), 0)

    orientation = None
    mode = 'Unknown'
    if beam == 1:
        mode = 'Down'
    elif beam == 2:
        mode = 'SideScan'
        orientation = 'Port'
    elif beam == 3:
        mode = 'SideScan'
        orientation = 'Starboard'
    elif beam == 4:
        mode = 'Down Imaging'

    label = "Channel {}".format(channel_id)
    if mode != 'Unknown':
        label = mode if orientation is None else "{} {}".format(mode, orientation)

    return {
        'label': label,
        'mode': mode,
        'orientation': orientation,
        'beam': beam,
        'startFrequencyHz': start_freq if start_freq > 0 else None,
        'endFrequencyHz': end_freq if end_freq > 0 else None,
    }


def _write_generic_sonar_data_player_project(sonar_obj, input_path: str, out_dir: str, file_map: dict = None):
    if not hasattr(sonar_obj, 'header_dat') or sonar_obj.header_dat is None:
        raise ValueError("PINGverter parser did not produce header_dat.")

    df = sonar_obj.header_dat.copy()
    if len(df) == 0:
        raise ValueError("No ping metadata rows were decoded.")

    if 'channel_id' not in df.columns:
        if 'beam' in df.columns:
            df['channel_id'] = pd.to_numeric(df['beam'], errors='coerce').fillna(0).astype(int)
        else:
            df['channel_id'] = 0

    if 'record_num' not in df.columns:
        df = df.reset_index(drop=True)
        df['record_num'] = df.index

    os.makedirs(out_dir, exist_ok=True)
    pings_csv = os.path.join(out_dir, 'pings.csv')
    samples_path = os.path.join(out_dir, 'samples.u16le')
    frames_path = os.path.join(out_dir, 'frames.jsonl')

    df.to_csv(pings_csv, index=False)

    frame_key = _frame_key_series(df)
    df['_frame_key'] = frame_key

    frame_count = 0
    sample_offset = 0

    # Open per-channel source file handles when file_map is provided (e.g. Humminbird
    # stores each beam in a separate .SON file), otherwise use a single input_path handle.
    if file_map:
        source_handles = {cid: open(path, 'rb') for cid, path in file_map.items()}
        default_source = None
    else:
        default_source = open(input_path, 'rb')
        source_handles = {}

    try:
        with open(samples_path, 'wb') as sample_file, open(frames_path, 'w', encoding='utf-8') as frame_file:
            grouped = df.groupby('_frame_key', dropna=False, sort=True)
            for key, group in grouped:
                channels = []
                sorted_group = group.sort_values('channel_id') if 'channel_id' in group.columns else group

                for _, row in sorted_group.iterrows():
                    index = _safe_int(_row_value(row, ['index', 'frame_offset']), -1)
                    son_offset = _safe_int(_row_value(row, ['son_offset']), -1)
                    sample_count = _safe_int(_row_value(row, ['sample_cnt', 'ping_cnt', 'num_results']), 0)
                    bytes_per_sample = _bytes_per_sample(sonar_obj, row)

                    if index < 0 or son_offset < 0 or sample_count <= 0 or bytes_per_sample <= 0:
                        continue

                    # Route read to the correct source file
                    row_channel_id = _safe_int(_row_value(row, ['channel_id', 'beam']), 0)
                    source = source_handles.get(row_channel_id, default_source)
                    if source is None:
                        continue

                    raw_len = sample_count * bytes_per_sample
                    source.seek(index + son_offset)
                    raw = source.read(raw_len)
                    if len(raw) != raw_len:
                        continue

                    values = _decode_raw_to_u16(raw, bytes_per_sample, sonar_obj)
                    if values.size == 0:
                        continue

                    data = values.tobytes(order='C')
                    sample_file.write(data)

                    channels.append({
                        'channelId': _safe_int(_row_value(row, ['channel_id', 'beam']), 0),
                        'sampleOffset': sample_offset,
                        'sampleCount': int(values.size),
                        'byteLength': len(data),
                        'minRangeMeters': _none_if_nan(_row_value(row, ['min_range', 'first_sample_depth'])),
                        'maxRangeMeters': _none_if_nan(_row_value(row, ['max_range', 'last_sample_depth'])),
                        'bottomDepthMeters': _none_if_nan(_row_value(row, ['inst_dep_m', 'bottom_depth'])),
                    })
                    sample_offset += len(data)

                if not channels:
                    continue

                frame = {
                    'frameIndex': frame_count,
                    'sequenceCount': _safe_int(key, frame_count),
                    'timeSeconds': _none_if_nan(group['time_s'].mean()) if 'time_s' in group.columns else None,
                    'lat': _none_if_nan(group['lat'].mean()) if 'lat' in group.columns else None,
                    'lon': _none_if_nan(group['lon'].mean()) if 'lon' in group.columns else None,
                    'speedMetersPerSecond': _none_if_nan(group['speed_ms'].mean()) if 'speed_ms' in group.columns else None,
                    'trackDistanceMeters': _none_if_nan(group['trk_dist'].mean()) if 'trk_dist' in group.columns else None,
                    'headingDegrees': _none_if_nan(group['instr_heading'].mean()) if 'instr_heading' in group.columns else None,
                    'temperatureCelsius': _none_if_nan(group['tempC'].mean()) if 'tempC' in group.columns else None,
                    'channels': channels,
                }
                frame_file.write(json.dumps(frame, separators=(',', ':')) + '\n')
                frame_count += 1
    finally:
        if default_source:
            default_source.close()
        for h in source_handles.values():
            h.close()

    manifest_channels = []
    channel_ids = sorted(int(c) for c in pd.to_numeric(df['channel_id'], errors='coerce').dropna().unique())
    for channel_id in channel_ids:
        group = df[pd.to_numeric(df['channel_id'], errors='coerce') == channel_id]
        if len(group) == 0:
            continue

        meta = _channel_metadata(group, channel_id)
        sample_col = 'sample_cnt' if 'sample_cnt' in group.columns else 'ping_cnt'
        max_samples = _safe_int(group[sample_col].max(), 0) if sample_col in group.columns else 0

        manifest_channels.append({
            'channelId': channel_id,
            'label': meta['label'],
            'mode': meta['mode'],
            'orientation': meta['orientation'],
            'beam': meta['beam'],
            'startFrequencyHz': meta['startFrequencyHz'],
            'endFrequencyHz': meta['endFrequencyHz'],
            'rows': int(len(group)),
            'maxSamples': max_samples,
            'timeStart': _none_if_nan(group['time_s'].min()) if 'time_s' in group.columns else None,
            'timeEnd': _none_if_nan(group['time_s'].max()) if 'time_s' in group.columns else None,
        })

    manifest = {
        'formatVersion': 2,
        'source': os.path.abspath(input_path),
        'telemetry': 'pings.csv',
        'frames': 'frames.jsonl',
        'samples': {
            'path': 'samples.u16le',
            'encoding': 'uint16-le',
        },
        'frameCount': frame_count,
        'channels': manifest_channels,
    }

    manifest_path = os.path.join(out_dir, 'manifest.json')
    with open(manifest_path, 'w', encoding='utf-8') as file:
        json.dump(manifest, file, indent=2)

    return manifest_path


def export_sonar_data_player_project(input: str, out_dir: str, include_pngs: bool=True,
                                     nchunk: int=500, tempC: float=10, exportUnknown: bool=True,
                                     source_format: str=None):
    """Export any supported raw recording into a SonarDataPlayer project folder."""
    assert os.path.isfile(input), "{} does not exist.".format(input)

    ext = _detect_format(input, source_format=source_format)
    os.makedirs(out_dir, exist_ok=True)

    if ext in ('.sl2', '.sl3'):
        sonar_obj = low(inFile=input, nchunk=nchunk, exportUnknown=exportUnknown)
        sonar_obj.tempC = float(tempC) / 10
        return sonar_obj.write_sonar_data_player_project(out_dir, include_pngs=include_pngs)

    parser_work_dir = os.path.join(out_dir, 'meta')
    os.makedirs(parser_work_dir, exist_ok=True)

    sonar_obj = _build_sonar_object(
        input_path=input,
        work_dir=parser_work_dir,
        ext=ext,
        nchunk=nchunk,
        tempC=tempC,
        exportUnknown=exportUnknown,
    )

    if hasattr(sonar_obj, 'write_sonar_data_player_project'):
        return sonar_obj.write_sonar_data_player_project(out_dir, include_pngs=include_pngs)

    # Humminbird: hum2pingmapper writes per-beam CSVs to disk rather than storing
    # header_dat in memory.  Reconstruct header_dat here and build a channel→file map
    # so the generic writer reads samples from the correct .SON file for each beam.
    file_map = None
    if ext == '.dat' and hasattr(sonar_obj, 'beamMeta'):
        beam_dfs = []
        file_map = {}
        for beam_key, meta in sonar_obj.beamMeta.items():
            csv_path = meta.get('metaCSV')
            son_path = meta.get('sonFile')
            if csv_path and os.path.isfile(csv_path):
                beam_df = pd.read_csv(csv_path)
                try:
                    channel_id = int(beam_key[1:])  # 'B000' -> 0, 'B001' -> 1, …
                except (ValueError, IndexError):
                    channel_id = len(beam_dfs)
                beam_df['channel_id'] = channel_id
                beam_dfs.append(beam_df)
                if son_path and os.path.isfile(son_path):
                    file_map[channel_id] = son_path
        if beam_dfs:
            sonar_obj.header_dat = pd.concat(beam_dfs, ignore_index=True)
        if not beam_dfs or not file_map:
            raise ValueError("Humminbird: could not load per-beam metadata CSVs from beamMeta.")
        # Build a robust frame key from time clustering. This tolerates non-zero
        # record starts and dropped records while still grouping same-event beams.
        sonar_obj.header_dat['sequence_cnt'] = _humminbird_sequence_from_time(sonar_obj.header_dat)
        # Ensure son8bit is set so _bytes_per_sample returns 1
        if not hasattr(sonar_obj, 'son8bit'):
            sonar_obj.son8bit = True

    return _write_generic_sonar_data_player_project(sonar_obj, input, out_dir, file_map=file_map)


SUPPORTED_SONAR_EXTENSIONS = (
    '.dat',
    '.jsf',
    '.rsd',
    '.sl2',
    '.sl3',
    '.svlog',
    '.xtf',
)


# =========================================================
# Lowrance to Humminbird
# =========================================================

def low2hum(input: str, output: str):

    # Make sure input exists
    assert os.path.isfile(input), "{} does not exist.".format(input)

    # Create the classes
    lowrance = low(input)
    humminbird = hum(output)

    # Make output files
    humminbird._makeOutFiles()

    # Start the decode
    start_time = time.time()
    print('\n\nDecoding Lowrance File...')

    # Get Lowrance file length
    lowrance._getFileLen()

    # Parse file header
    lowrance._parseFileHeader()

    # Parse ping headers (attributes)
    lowrance._parsePingHeader()

    # Convert ping attributes to known units
    lowrance._convertPingAttributes()

    # lowrance.header_dat.to_csv('lowrance_test.csv')

    print("Time (s):", round(time.time() - start_time, ndigits=1))

    #######################
    # Convert to Humminbird
    #######################

    start_time = time.time()
    print('\n\nConverting to Humminbird...')

    # Convert to Humminbird attributes
    humminbird._convertLowHeader(lowrance)

    # Drop unknown beams
    humminbird._removeUnknownBeams()

    # Split sidescan, if necessary
    beams = humminbird.header_dat['beam'].unique()
    if 5 in beams:
        humminbird._splitLowSS()
        flip_port = True
    else:
        flip_port = False

    # Recalculate record number
    humminbird._recalcRecordNum()

    # humminbird.header_dat.to_csv('hum_converted.csv')

    # Get necessary data for DAT file
    humminbird._convertLowDAT(lowrance)

    # Write DAT to file
    humminbird._writeDAT()

    # Save b001 and b002
    beams = humminbird.header_dat['beam'].unique()

    if 0 in beams:
        humminbird._writeSonfromLow(0, lowrance.frame_header_size, lowrance.path)

    if 1 in beams:
        humminbird._writeSonfromLow(1, lowrance.frame_header_size, lowrance.path)

    if 2 in beams:
        humminbird._writeSonfromLow(2, lowrance.frame_header_size, lowrance.path, flip_port)

    if 3 in beams:
        humminbird._writeSonfromLow(3, lowrance.frame_header_size, lowrance.path)

    if 4 in beams:
        humminbird._writeSonfromLow(4, lowrance.frame_header_size, lowrance.path)

    # Split b005 (lowrance sidescan) into port (2) and star (3)


    humminbird.header_dat.to_csv('hum_converted.csv')

    print("Time (s):", round(time.time() - start_time, ndigits=1))

    return



