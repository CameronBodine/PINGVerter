
import sys, os
import numpy as np
import pandas as pd

def filterGPS(df: pd.DataFrame, 
              jump_thresh: float=1):

    ####################
    # Add filtering step
    # Convert degrees to radians
    lat_rad = np.radians(df['lat'])
    lon_rad = np.radians(df['lon'])

    # Earth's radius in meters
    R = 6371000

    # Compute deltas
    dlat = lat_rad.diff()
    dlon = lon_rad.diff()

    # Haversine formula
    a = np.sin(dlat / 2)**2 + np.cos(lat_rad.shift()) * np.cos(lat_rad) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    dist_m = R * c

    df['gps_jump_m'] = dist_m.fillna(0)

    
    df['gps_bad'] = df['gps_jump_m'] > jump_thresh

    df.loc[df['gps_bad'], ['lat', 'lon']] = np.nan

    # Also filter by percentile
    lat_low, lat_high = np.percentile(df['lat'].dropna(), [0.5, 99.5])
    lon_low, lon_high = np.percentile(df['lon'].dropna(), [0.5, 99.5])

    valid_lat = df['lat'].between(lat_low, lat_high)
    valid_lon = df['lon'].between(lon_low, lon_high)
    valid_coords = valid_lat & valid_lon

    df.loc[~valid_coords, ['lat', 'lon']] = np.nan

    df['lat'] = df['lat'].interpolate(method='linear', limit_direction='both')
    df['lon'] = df['lon'].interpolate(method='linear', limit_direction='both')

    # df['gps_interpolated'] = df['lat'].isna() | df['lon'].isna()

    df.drop(['gps_bad', 'gps_jump_m'], axis=1, inplace=True)


    return df