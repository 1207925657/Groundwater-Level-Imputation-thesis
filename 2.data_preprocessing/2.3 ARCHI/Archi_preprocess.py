import argparse
import codecs
import hashlib
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd


START = pd.Timestamp("1990-01-01")
SGI_THRESHOLD = -1.0
DATE_FORMAT = "%m/%d/%Y"
DATA_COLUMNS = ["site_no", "date", "value"]
METADATA_COLUMNS = ["site_no", "latitude", "longitude", "group"]
HEADER = codecs.BOM_UTF8 + b"site_no,date,value\n"

EXCLUDED = {
    "Malgomaj_41", "Malgomaj_44", "Malgomaj_45", "Sollefteå_2",
    "Ströms_Vattudal_22", "Vaxholm_6", "Brattforsheden_3", "Kolmården_3",
    "Lysekil_4", "Motala_31", "Ormsjön_11", "Ronneby_101", "Torpshammar_24",
    "Torpshammar_28", "Vellinge_8", "Vimmerby_54",
}


def read_sources(gwl_path, sgi_path):

    if not sgi_path.is_file():
        raise FileNotFoundError(f"Run 2.1 SGI/preprocess_sgi.py first: {sgi_path}")
    sgi = pd.read_excel(sgi_path, sheet_name="Location_SGI_Final",
                        usecols=["station_name", "date", "Latitude", "Longitude", "sgi"])
    sgi["station_name"] = sgi.station_name.astype("string").str.strip()
    sgi["date"] = pd.to_datetime(sgi.date, errors="raise")
    for column in ["Latitude", "Longitude", "sgi"]:
        sgi[column] = pd.to_numeric(sgi[column], errors="raise")
    if sgi[["station_name", "date", "Latitude", "Longitude"]].isna().any().any():
        raise ValueError("Missing SGI station, date or coordinates.")
    if sgi.station_name.eq("").any():
        raise ValueError("Empty station identifier.")
    sgi = sgi.loc[~sgi.station_name.isin(EXCLUDED)].copy()
    sgi["year_month"] = sgi.date.dt.to_period("M")
    if sgi.duplicated(["station_name", "year_month"]).any():
        raise ValueError("Duplicate SGI station/month keys.")
    meta = sgi[["station_name", "Latitude", "Longitude"]].drop_duplicates()
    if meta.station_name.duplicated().any():
        raise ValueError("Conflicting coordinates for a station.")
    if not (meta.Latitude.between(-90, 90) & meta.Longitude.between(-180, 180)).all():
        raise ValueError("Invalid latitude or longitude.")
    meta = meta.rename(columns={"station_name": "site_no", "Latitude": "latitude",
                                "Longitude": "longitude"}).sort_values("site_no")
    meta["group"] = "all_sites"

    gwl = pd.read_csv(gwl_path, usecols=["station_name", "date", "gwl_(inf)"])
    gwl["station_name"] = gwl.station_name.astype("string").str.strip()
    gwl = gwl.loc[gwl.station_name.isin(meta.site_no)].copy()
    gwl["date"] = pd.to_datetime(gwl.date, errors="raise")
    gwl["gwl_(inf)"] = pd.to_numeric(gwl["gwl_(inf)"], errors="raise")
    if gwl.date.isna().any() or (gwl.date != gwl.date.dt.normalize()).any():
        raise ValueError("GWL dates must be valid daily dates without a time component.")
    if np.isinf(gwl["gwl_(inf)"].to_numpy()).any():
        raise ValueError("Infinite groundwater value.")
    duplicate_rows = gwl.loc[gwl.duplicated(["station_name", "date"], keep=False)]
    duplicate_report = []
    for well, records in duplicate_rows.groupby("station_name", sort=True):
        counts = records.groupby("date")["gwl_(inf)"].nunique(dropna=False)
        duplicate_report.append({"site_no": well, "duplicate_dates": len(counts),
                                 "conflicting_dates": int(counts.gt(1).sum()),
                                 "discarded_rows": len(records) - len(counts),
                                 "kept_source_occurrence": "first"})

    gwl = gwl.drop_duplicates(["station_name", "date"], keep="first")
    gwl = gwl.loc[gwl.date.ge(START)].sort_values(["station_name", "date"])
    absent = set(meta.site_no) - set(gwl.station_name)
    if absent:
        raise ValueError(f"Missing groundwater series: {sorted(absent)}")
    report = pd.DataFrame(duplicate_report, columns=["site_no", "duplicate_dates",
                          "conflicting_dates", "discarded_rows", "kept_source_occurrence"])
    return gwl, sgi, meta[METADATA_COLUMNS], report


def prepare_well(well, records, monthly_sgi):

    series = records.set_index("date")["gwl_(inf)"].sort_index()
    end = series.last_valid_index()
    if end is None:
        raise ValueError(f"No observed groundwater values: {well}")
    dates = pd.date_range(START, end, freq="D")
    values = series.reindex(dates).to_numpy()
    sgi = monthly_sgi.reindex(dates.to_period("M")).to_numpy()
    if (pd.notna(values) & pd.isna(sgi)).any():
        raise ValueError(f"Observed groundwater without matching monthly SGI: {well}")
    drought = sgi < SGI_THRESHOLD
    removed = drought & pd.notna(values)
    original = pd.DataFrame({"site_no": well, "date": dates, "value": values})
    masked = original.copy()
    masked.loc[drought, "value"] = np.nan
    audit = original.loc[removed].rename(columns={"value": "original_value"}).copy()
    audit["sgi"] = sgi[removed]
    audit["was_removed"] = True
    summary = {"site_no": well, "start_date": dates[0], "eval_end": end,
               "rows": len(dates), "observed_values": int(pd.notna(values).sum()),
               "natural_missing_values": int(pd.isna(values).sum()),
               "masked_target_rows": int(removed.sum()),
               "target_last_observed_after_mask": masked.loc[masked.value.notna(), "date"].max()}
    return original, masked, audit, summary


def csv_block(data):

    return data[DATA_COLUMNS].to_csv(index=False, header=False, na_rep="",
                                    date_format=DATE_FORMAT, lineterminator="\n").encode("utf-8")


def write_case(path, wells, original_blocks, masked_blocks, target):

    digest = hashlib.sha256()
    with path.open("wb") as stream:
        stream.write(HEADER)
        digest.update(HEADER)
        for well in wells:
            block = masked_blocks[well] if well == target else original_blocks[well]
            stream.write(block)
            digest.update(block)
    return digest.hexdigest()


def save_csv(data, path):
    data.to_csv(path, index=False, encoding="utf-8-sig", na_rep="",
                date_format=DATE_FORMAT, lineterminator="\n")


def run(args):
    gwl, sgi, metadata, duplicates = read_sources(args.gwl, args.sgi)
    wells = metadata.site_no.tolist()
    targets = wells if args.targets is None else sorted(set(args.targets))
    unknown = set(targets) - set(wells)
    if unknown:
        raise ValueError(f"Unavailable or excluded targets: {sorted(unknown)}")
    for well in wells:
        if any(character in well for character in '<>:"/\\|?*') or well.endswith((".", " ")):
            raise ValueError(f"Station name cannot be used as a filename: {well}")
    if len({well.casefold() for well in wells}) != len(wells):
        raise ValueError("Case-insensitive station filename collision.")

    monthly = {well: part.set_index("year_month").sgi for well, part in sgi.groupby("station_name")}
    original_blocks, masked_blocks, audits, summaries = {}, {}, [], []
    for well, records in gwl.groupby("station_name", sort=True):
        original, masked, audit, summary = prepare_well(well, records, monthly[well])
        original_blocks[well] = csv_block(original)
        if well in targets:
            masked_blocks[well] = csv_block(masked)
            audits.append(audit)
        summaries.append(summary)
    stats = pd.DataFrame(summaries).set_index("site_no").loc[wells]
    base_bytes = len(HEADER) + sum(map(len, original_blocks.values()))
    output_bytes = sum(base_bytes - len(original_blocks[t]) + len(masked_blocks[t]) for t in targets)
    print(f"Prepared {len(wells)} wells, {stats.rows.sum():,} rows per case; "
          f"{len(targets)} target cases, approximately {output_bytes / 1024**3:.2f} GiB of input CSVs.", flush=True)
    if args.plan_only:
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)
    case_dir = args.output_dir / "Scenario_A"
    case_dir.mkdir(exist_ok=True)
    if shutil.disk_usage(case_dir).free < output_bytes + 100 * 1024**2:
        raise OSError("Insufficient free space. Use --targets to export a smaller target subset.")
    save_csv(metadata, case_dir / "site_metadata.csv")
    save_csv(stats.reset_index(), args.output_dir / "well_summary.csv")
    save_csv(duplicates, args.output_dir / "source_duplicate_summary.csv")
    save_csv(pd.concat(audits, ignore_index=True), args.output_dir / "mask_records.csv")
    manifest = []
    for index, target in enumerate(targets, 1):
        path = case_dir / f"{target}_ARCHI_input.csv"
        digest = write_case(path, wells, original_blocks, masked_blocks, target)
        info = stats.loc[target]
        manifest.append({"target_well": target, "output_file": path.name,
                         "total_rows": int(stats.rows.sum()), "total_wells": len(wells),
                         "target_rows": int(info.rows), "masked_target_rows": int(info.masked_target_rows),
                         "target_first_date": info.start_date, "eval_end": info.eval_end,
                         "target_last_observed_after_mask": info.target_last_observed_after_mask,
                         "bytes": path.stat().st_size, "sha256": digest})
        if index % 10 == 0 or index == len(targets):
            print(f"Exported {index}/{len(targets)} Scenario A inputs.", flush=True)
    save_csv(pd.DataFrame(manifest), case_dir / "ARCHI_1masked_summary.csv")
    print(f"Complete. ARCHI inputs: {case_dir}", flush=True)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    here = Path(__file__).resolve().parent
    repository = here.parents[1]
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gwl", type=Path, default=repository / "1.Input data/EDGE_SE_GWL.csv")


    parser.add_argument("--sgi", type=Path, default=here.parent / "2.1 SGI/Location_SGI_final.xlsx")
    parser.add_argument("--output-dir", type=Path, default=here)
    parser.add_argument("--targets", nargs="+", help="Only these target cases; retain all wells as reference candidates")
    parser.add_argument("--plan-only", action="store_true", help="Report counts and storage size without exporting")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
