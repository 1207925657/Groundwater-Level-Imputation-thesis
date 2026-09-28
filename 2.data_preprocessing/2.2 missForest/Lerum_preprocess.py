import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


START = pd.Timestamp("1990-01-01")
SGI_THRESHOLD = -1.0
RADIUS_KM = 50
MIN_OVERLAP_DAYS = 365
MIN_HISTORY_YEARS = 3
SEASONAL_HALF_WINDOW = 15
MAX_NEARBY = 10
GWL_DUPLICATE_POLICY = "first"
TIME = ["day_of_year", "month", "Unix_time"]
WELL_LEVELS = [1, 2, 3, 5, 8, 10]

EXCLUDED = {
    "Malgomaj_41", "Malgomaj_44", "Malgomaj_45", "Sollefteå_2",
    "Ströms_Vattudal_22", "Vaxholm_6", "Brattforsheden_3", "Kolmården_3",
    "Lysekil_4", "Motala_31", "Ormsjön_11", "Ronneby_101", "Torpshammar_24",
    "Torpshammar_28", "Vellinge_8", "Vimmerby_54",
}


def distance_km(lat, lon, other_lat, other_lon):
    lat, lon, other_lat, other_lon = map(np.radians, (lat, lon, other_lat, other_lon))
    a = (np.sin((other_lat - lat) / 2) ** 2
         + np.cos(lat) * np.cos(other_lat) * np.sin((other_lon - lon) / 2) ** 2)
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def deduplicate_gwl(gwl):

    keys = ["station_name", "date"]
    duplicates = gwl.loc[gwl.duplicated(keys, keep=False)]
    report = []
    for station, rows in duplicates.groupby("station_name", sort=True):
        counts = rows.groupby("date")["gwl_(inf)"].nunique(dropna=False)
        report.append({"station_name": station, "duplicate_dates": len(counts),
                       "conflicting_dates": int(counts.gt(1).sum()),
                       "discarded_rows": len(rows) - len(counts),
                       "kept_source_occurrence": GWL_DUPLICATE_POLICY})
    cleaned = gwl.drop_duplicates(keys, keep=GWL_DUPLICATE_POLICY)
    return cleaned, report


def read_sources(gwl_path, sgi_path, target_prefix="Lerum_", only_targets=None):


    if not sgi_path.is_file():
        raise FileNotFoundError(f"Run 2.1 SGI/preprocess_sgi.py first. Required output: {sgi_path}")
    columns = ["station_name", "date", "Longitude", "Latitude", "sgi"]
    sgi = pd.read_excel(sgi_path, sheet_name="Location_SGI_Final", usecols=columns)
    sgi["station_name"] = sgi["station_name"].astype("string").str.strip()
    sgi["date"] = pd.to_datetime(sgi["date"], errors="raise")
    for col in ["Longitude", "Latitude", "sgi"]:
        sgi[col] = pd.to_numeric(sgi[col], errors="raise")
    if sgi[columns[:-1]].isna().any().any():
        raise ValueError("Missing SGI station name, date or coordinates.")
    sgi = sgi.loc[~sgi.station_name.isin(EXCLUDED)].copy()
    sgi["year_month"] = sgi.date.dt.to_period("M")
    if sgi.duplicated(["station_name", "year_month"]).any():
        raise ValueError("Duplicate SGI station/month keys.")
    meta = sgi[["station_name", "Latitude", "Longitude"]].drop_duplicates()
    if meta.station_name.duplicated().any():
        raise ValueError("Conflicting coordinates for a station.")
    meta = meta.set_index("station_name")
    targets = sorted(name for name in meta.index
                     if target_prefix is None or name.startswith(target_prefix))
    if only_targets is not None:
        unknown = set(only_targets) - set(targets)
        if unknown:
            raise ValueError(f"Requested target wells are unavailable/excluded: {sorted(unknown)}")
        targets = sorted(only_targets)
    if not targets:
        raise ValueError("No eligible target wells found in the SGI file.")
    needed = set(targets)
    for name in targets:
        lat, lon = meta.loc[name, ["Latitude", "Longitude"]]
        nearby = distance_km(lat, lon, meta.Latitude, meta.Longitude) <= RADIUS_KM
        needed.update(meta.index[nearby])

    frames = []
    for chunk in pd.read_csv(gwl_path, usecols=["station_name", "date", "gwl_(inf)"],
                             chunksize=300_000):
        chunk["station_name"] = chunk.station_name.astype("string").str.strip()
        frames.append(chunk.loc[chunk.station_name.isin(needed)])
    gwl = pd.concat(frames, ignore_index=True)
    gwl["date"] = pd.to_datetime(gwl.date, errors="raise")
    gwl["gwl_(inf)"] = pd.to_numeric(gwl["gwl_(inf)"], errors="raise")
    if gwl.date.isna().any():
        raise ValueError("Missing GWL dates.")
    gwl, duplicates = deduplicate_gwl(gwl)
    if duplicates:
        print(f"Duplicate source records: kept the {GWL_DUPLICATE_POLICY} occurrence "
              f"for {sum(row['duplicate_dates'] for row in duplicates)} station/date keys.", flush=True)
    gwl = gwl.loc[gwl.date >= START].sort_values(["station_name", "date"])
    wide = gwl.pivot(index="date", columns="station_name", values="gwl_(inf)")
    wide.attrs["source_duplicates"] = duplicates
    absent = set(targets) - set(wide.columns)
    if absent:
        raise ValueError(f"Missing GWL target wells: {sorted(absent)}")
    return sgi, meta, targets, wide


def mask_target(well, wide, sgi):

    end = wide[well].last_valid_index()
    if end is None:
        raise ValueError(f"No valid GWL observations: {well}")
    dates = pd.date_range(START, end, freq="D", name="date")
    monthly = sgi.loc[sgi.station_name.eq(well)].set_index("year_month").sgi
    table = pd.DataFrame(index=dates)
    table["target_gwl_observed"] = wide[well].reindex(dates)
    table["sgi"] = monthly.reindex(dates.to_period("M")).to_numpy()
    if (table.target_gwl_observed.notna() & table.sgi.isna()).any():
        raise ValueError(f"Observed GWL without matching monthly SGI: {well}")
    table["is_sgi_drought"] = table.sgi.lt(SGI_THRESHOLD)
    table["was_removed"] = table.is_sgi_drought & table.target_gwl_observed.notna()
    table["target_gwl_original"] = table.target_gwl_observed.mask(table.is_sgi_drought)
    return table


def rank_nearby(well, target, wide, meta):

    lat, lon = meta.loc[well, ["Latitude", "Longitude"]]
    distances = distance_km(lat, lon, meta.Latitude, meta.Longitude)
    records = []
    for candidate in meta.index[(distances <= RADIUS_KM) & (meta.index != well)]:
        if candidate not in wide:
            continue
        pairs = pd.concat([target.rename("target"), wide[candidate].rename("reference")],
                          axis=1, join="inner").dropna()
        if len(pairs) < MIN_OVERLAP_DAYS:
            continue
        if pairs.target.nunique() < 2 or pairs.reference.nunique() < 2:
            continue

        rho = pairs.target.corr(pairs.reference, method="pearson")
        records.append({"target_well": well, "reference_well": candidate,
                        "distance_km": distances.loc[candidate], "pearson": rho,
                        "abs_pearson": abs(rho), "overlap_days": len(pairs)})
    if not records:
        return pd.DataFrame(columns=["target_well", "reference_well", "distance_km",
                                     "pearson", "abs_pearson", "overlap_days", "rank"])
    ranking = pd.DataFrame(records).sort_values(
        ["abs_pearson", "distance_km", "reference_well"], ascending=[False, True, True]
    ).head(MAX_NEARBY).reset_index(drop=True)
    ranking["rank"] = np.arange(1, len(ranking) + 1)
    return ranking


def seasonal_features(series):


    dates = series.index
    seasonal_dates = pd.to_datetime({"year": np.full(len(dates), 2001),
                                    "month": dates.month,
                                    "day": np.where((dates.month == 2) & (dates.day == 29),
                                                    28, dates.day)})
    days = seasonal_dates.dt.dayofyear.to_numpy() - 1
    values = series.to_numpy(dtype=float)
    valid = np.isfinite(values)
    anchor = values[valid][0] if valid.any() else 0.0
    centered = values - anchor
    anomaly = np.full(len(series), np.nan)
    zscore = np.full(len(series), np.nan)
    count, total, squares = np.zeros((3, 365))
    history_years = 0
    width = SEASONAL_HALF_WINDOW
    kernel = np.ones(2 * width + 1)
    for year in sorted(set(dates.year)):
        current = dates.year == year
        if history_years >= MIN_HISTORY_YEARS:
            n, sums, sumsq = [np.convolve(np.r_[x[-width:], x, x[:width]], kernel,
                                         mode="valid") for x in [count, total, squares]]
            mean = np.divide(sums, n, out=np.full(365, np.nan), where=n > 0)
            variance = np.divide(sumsq - sums * mean, n - 1,
                                 out=np.full(365, np.nan), where=n > 1)
            sd = np.sqrt(np.maximum(variance, 0))
            delta = centered[current] - mean[days[current]]
            anomaly[current] = delta
            daily_sd = sd[days[current]]
            zscore[current] = np.divide(delta, daily_sd,
                                        out=np.full(delta.shape, np.nan), where=daily_sd > 0)
        usable = current & valid
        if usable.any():
            count += np.bincount(days[usable], minlength=365)
            total += np.bincount(days[usable], weights=centered[usable], minlength=365)
            squares += np.bincount(days[usable], weights=centered[usable] ** 2, minlength=365)
            history_years += 1
    return anomaly, zscore


def np_features(climate):

    t, rh, wind = climate.tg, climate.hu, climate.fg
    es = 0.6108 * np.exp(17.27 * t / (t + 237.3))
    ea = es * rh / 100
    slope = 4098 * es / (t + 237.3) ** 2
    gamma = 0.665e-3 * 101.3
    solar_radiation = climate.qq * 0.0864
    radiation = 0.77 * solar_radiation
    pet = (0.408 * slope * radiation + gamma * (900 / (t + 273)) * wind * (es - ea))
    pet = pet / (slope + gamma * (1 + 0.34 * wind))
    net = (climate.rr - pet).clip(lower=0)
    return pd.DataFrame({f"NP_{m}m": net.rolling(m * 30, min_periods=1).sum()
                         for m in range(1, 13)})


def load_climate(root, meta, targets, end):


    dates = pd.date_range(START, end, name="date")
    variables = {}
    grid_records = []
    for var in ["rr", "tg", "hu", "fg", "qq"]:
        frames = []
        files = sorted(root.rglob(f"{var}_*.nc"))
        if not files:
            raise FileNotFoundError(f"No E-OBS {var} NetCDF files in {root}")
        for path in files:
            print(f"Reading climate: {path.parent.name}/{path.name}", flush=True)
            with xr.open_dataset(path, engine="netcdf4", cache=False) as ds:
                lat_name = "latitude" if "latitude" in ds.coords else "lat"
                lon_name = "longitude" if "longitude" in ds.coords else "lon"
                if var not in ds:
                    raise ValueError(f"Variable {var} absent from {path}")
                latitudes, longitudes = ds[lat_name].values, ds[lon_name].values
                target_lon = meta.loc[targets, "Longitude"].to_numpy().copy()
                if longitudes.max() > 180:
                    target_lon = target_lon % 360
                lat_idx = pd.Index(latitudes).get_indexer(meta.loc[targets, "Latitude"], method="nearest")
                lon_idx = pd.Index(longitudes).get_indexer(target_lon, method="nearest")
                lat_start, lat_stop = int(lat_idx.min()), int(lat_idx.max()) + 1
                lon_start, lon_stop = int(lon_idx.min()), int(lon_idx.max()) + 1
                time = pd.DatetimeIndex(ds.time.values).normalize()
                wanted = np.flatnonzero((time >= START) & (time <= end))
                for first in range(0, len(wanted), 90):
                    ids = wanted[first:first + 90]
                    block = ds[var].isel({"time": slice(int(ids[0]), int(ids[-1]) + 1),
                                          lat_name: slice(lat_start, lat_stop),
                                          lon_name: slice(lon_start, lon_stop)})
                    values = block.transpose("time", lat_name, lon_name).values
                    selected = values[:, lat_idx - lat_start, lon_idx - lon_start]
                    frames.append(pd.DataFrame(selected, index=time[ids], columns=targets))
                for i, well in enumerate(targets):
                    grid_records.append({"target_well": well, "variable": var,
                                         "file": path.relative_to(root).as_posix(),
                                         "grid_latitude": float(latitudes[lat_idx[i]]),
                                         "grid_longitude": float(longitudes[lon_idx[i]]),
                                         "units": ds[var].attrs.get("units", "")})
        if not frames:
            raise ValueError(f"No {var} dates overlap the target period.")
        combined = pd.concat(frames).sort_index()
        if combined.index.duplicated().any():
            raise ValueError(f"Overlapping E-OBS dates: {var}")
        variables[var] = combined.reindex(dates)
    output, coverage = {}, []
    for well in targets:
        raw = pd.DataFrame({var: frame[well] for var, frame in variables.items()}, index=dates)
        output[well] = np_features(raw)
        for name, values in pd.concat([raw, output[well]], axis=1).items():
            coverage.append({"target_well": well, "variable": name,
                             "total_days": len(values), "missing_days": int(values.isna().sum()),
                             "first_valid_date": values.first_valid_index(),
                             "last_valid_date": values.last_valid_index()})
    print(f"Climate prepared for {len(targets)} wells; missing values recorded in climate_coverage.csv.",
          flush=True)
    return output, pd.DataFrame(grid_records), pd.DataFrame(coverage)


def experiment_configs():

    records = []

    def add(group, name, columns, wells=0):
        records.append({"experiment": group, "configuration": name,
                        "reference_wells": wells, "number_of_features": len(columns),
                        "feature_columns": ";".join(columns)})

    for name, columns in [("T1_day_of_year", TIME[:1]), ("T2_month", TIME[1:2]),
                          ("T3_Unix_time", TIME[2:]), ("T4_all_time_variables", TIME)]:
        add("Experiment_1_Time", name, columns)
    for n in range(1, 13):
        add("Experiment_2_NP", f"NP01_to_NP{n:02d}", [f"NP_{m}m" for m in range(1, n + 1)])
    np_sets = {f"NP01_to_NP{n:02d}": [f"NP_{m}m" for m in range(1, n + 1)]
               for n in range(1, 7)}
    np_sets["NP01_to_NP06_Plus_NP12"] = [f"NP_{m}m" for m in range(1, 7)] + ["NP_12m"]
    for label, columns in np_sets.items():
        add("Experiment_2_5_NP_Time", f"{label}_Plus_Time", columns + TIME)
    branches = {"GWL": ["gwl"], "ANOMALY": ["anomaly"], "ZSCORE": ["zscore"],
                "GWL_PLUS_ANOMALY": ["gwl", "anomaly"],
                "GWL_PLUS_ZSCORE": ["gwl", "zscore"],
                "GWL_PLUS_ANOMALY_PLUS_ZSCORE": ["gwl", "anomaly", "zscore"]}
    for n in range(1, 11):
        for label, kinds in branches.items():
            cols = [f"near_rank{rank:02d}_{kind}" for kind in kinds for rank in range(1, n + 1)]
            add("Experiment_3_NearbyWells", f"Top{n:02d}_{label}", cols, n)
    for n in WELL_LEVELS:
        wells = [f"near_rank{rank:02d}_gwl" for rank in range(1, n + 1)]
        add("Experiment_4_Wells_Time", f"Top{n:02d}_GWL_Plus_Time", wells + TIME, n)
        for label, columns in np_sets.items():
            add("Experiment_5_Wells_NP", f"Top{n:02d}_GWL_{label}", wells + columns, n)
            add("Experiment_6_Wells_NP_Time", f"Top{n:02d}_GWL_{label}_Plus_Time",
                wells + columns + TIME, n)
    return pd.DataFrame(records)


def save_csv(table, path, index=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=index, encoding="utf-8-sig", date_format="%Y-%m-%d")


def run(gwl_path, sgi_path, climate_dir, output_dir, export_experiments=False, prepared=None):
    print("Reading SGI and groundwater inputs...", flush=True)
    if prepared is None:
        sgi, meta, targets, wide = read_sources(gwl_path, sgi_path)
    else:
        sgi, meta, targets, wide, climate, grids, coverage = prepared
    masks = {well: mask_target(well, wide, sgi) for well in targets}
    end = max(table.index[-1] for table in masks.values())
    if prepared is None:
        climate, grids, coverage = load_climate(climate_dir, meta, targets, end)
    configs = experiment_configs()
    rankings, summaries, manifest = [], [], []
    derived_cache = {}
    for well in targets:
        audit = masks[well]
        data = audit[["target_gwl_original"]].copy()
        data["day_of_year"], data["month"] = data.index.dayofyear, data.index.month
        data["Unix_time"] = (data.index - START).days
        data = data.join(climate[well])
        ranking = rank_nearby(well, data.target_gwl_original, wide, meta)
        rankings.append(ranking)
        for row in ranking.itertuples():
            reference = row.reference_well
            if reference not in derived_cache:
                series = wide[reference].reindex(pd.date_range(START, end, name="date"))
                anomaly, zscore = seasonal_features(series)
                derived_cache[reference] = pd.DataFrame(
                    {"gwl": series, "anomaly": anomaly, "zscore": zscore}, index=series.index)
            features = derived_cache[reference].add_prefix(f"near_rank{row.rank:02d}_")
            data = data.join(features)
        save_csv(data, output_dir / "data" / f"{well}.csv", index=True)
        save_csv(audit, output_dir / "mask_records" / f"{well}.csv", index=True)
        for config in configs.itertuples():
            available = config.reference_wells <= len(ranking)
            columns = config.feature_columns.split(";")
            entry = {"target_well": well, "experiment": config.experiment,
                     "configuration": config.configuration,
                     "status": "ready" if available else "insufficient_reference_wells",
                     "data_file": f"data/{well}.csv"}
            if available:
                selected = data[["target_gwl_original"] + columns]
                entry["predictor_missing_cells"] = int(selected[columns].isna().sum().sum())
                entry["all_missing_predictor_columns"] = ";".join(
                    column for column in columns if selected[column].isna().all())
                if export_experiments:
                    relative = Path("experiments") / config.experiment / f"{well}__{config.configuration}.csv"
                    save_csv(selected, output_dir / relative, index=True)
                    entry["experiment_file"] = relative.as_posix()
            manifest.append(entry)
        summaries.append({"target_well": well, "start_date": audit.index[0],
                          "end_date": audit.index[-1], "rows": len(audit),
                          "observed_values": int(audit.target_gwl_observed.notna().sum()),
                          "removed_values": int(audit.was_removed.sum()),
                          "natural_missing_values": int(audit.target_gwl_observed.isna().sum()),
                          "masked_target_missing_values": int(data.target_gwl_original.isna().sum()),
                          "reference_wells": len(ranking),
                          "available_configurations": int((configs.reference_wells <= len(ranking)).sum())})
        print(f"Prepared {well}: {len(audit):,} days, {audit.was_removed.sum():,} removed, "
              f"{len(ranking)} reference wells.", flush=True)
    save_csv(configs, output_dir / "experiment_configurations.csv")
    save_csv(pd.DataFrame(manifest), output_dir / "experiment_manifest.csv")
    save_csv(pd.concat(rankings, ignore_index=True), output_dir / "nearby_well_rankings.csv")
    save_csv(pd.DataFrame(summaries), output_dir / "well_summary.csv")
    save_csv(grids, output_dir / "climate_grid_points.csv")
    save_csv(coverage, output_dir / "climate_coverage.csv")
    print(f"Saved {len(targets)} wells and {len(configs)} configurations to {output_dir}", flush=True)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    here = Path(__file__).resolve().parent
    repository = here.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gwl", type=Path, default=repository / "1.Input data/EDGE_SE_GWL.csv")
    parser.add_argument("--sgi", type=Path, default=here.parent / "2.1 SGI/Location_SGI_final.xlsx")
    parser.add_argument("--climate-dir", type=Path, default=repository / "1.Input data/Eobs")
    parser.add_argument("--output-dir", type=Path, default=here / "Lerum")
    parser.add_argument("--export-experiments", action="store_true")
    args = parser.parse_args()
    run(args.gwl, args.sgi, args.climate_dir, args.output_dir, args.export_experiments)


if __name__ == "__main__":
    main()
