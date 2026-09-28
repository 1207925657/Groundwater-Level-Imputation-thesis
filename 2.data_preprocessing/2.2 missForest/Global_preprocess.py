import argparse
import sys
from pathlib import Path

import pandas as pd
import preprocess_lerum as prep


def experiment_configs():


    specifications = [
        ("B1_TimeOnly", "Baseline", 0, 0, False, True),
        ("B2_NP01to03_Time", "Baseline", 0, 3, False, True),
        ("B3_NP01to06_Time", "Baseline", 0, 6, False, True),
        ("B4_NP01to06_NP12_Time", "Baseline", 0, 6, True, True),
        ("B5_Top01_GWL", "Baseline", 1, 0, False, False),
        ("B6_Top03_GWL", "Baseline", 3, 0, False, False),
        ("B7_Top05_GWL", "Baseline", 5, 0, False, False),
        ("B8_Top08_GWL", "Baseline", 8, 0, False, False),
        ("B9_Top10_GWL", "Baseline", 10, 0, False, False),
        ("M1_Top01_NP01_Time", "Main_route", 1, 1, False, True),
        ("M2_Top02_NP01to02_Time", "Main_route", 2, 2, False, True),
        ("M3_Top03_NP01to03_Time", "Main_route", 3, 3, False, True),
        ("M4_Top05_NP01to04_Time", "Main_route", 5, 4, False, True),
        ("M5_Top08_NP01to05_Time", "Main_route", 8, 5, False, True),
        ("M6_Top10_NP01to06_Time", "Main_route", 10, 6, False, True),
        ("W1_Top01_NP01to03_Time", "Fixed_NP01to03", 1, 3, False, True),
        ("W2_Top02_NP01to03_Time", "Fixed_NP01to03", 2, 3, False, True),
        ("W3_Top05_NP01to03_Time", "Fixed_NP01to03", 5, 3, False, True),
        ("W4_Top08_NP01to03_Time", "Fixed_NP01to03", 8, 3, False, True),
        ("W5_Top10_NP01to03_Time", "Fixed_NP01to03", 10, 3, False, True),
        ("N1_Top05_NP01_Time", "Fixed_Top05", 5, 1, False, True),
        ("N2_Top05_NP01to02_Time", "Fixed_Top05", 5, 2, False, True),
        ("N3_Top05_NP01to05_Time", "Fixed_Top05", 5, 5, False, True),
        ("N4_Top05_NP01to06_Time", "Fixed_Top05", 5, 6, False, True),
        ("C1_Top01_NP01to06_Time", "Substitution", 1, 6, False, True),
        ("C2_Top03_NP01to06_Time", "Substitution", 3, 6, False, True),
        ("C3_Top08_NP01to02_Time", "Substitution", 8, 2, False, True),
        ("C4_Top10_NP01_Time", "Substitution", 10, 1, False, True),
        ("L1_Top01_NP01to06_NP12_Time", "NP12_control", 1, 6, True, True),
        ("L2_Top05_NP01to06_NP12_Time", "NP12_control", 5, 6, True, True),
        ("L3_Top10_NP01to06_NP12_Time", "NP12_control", 10, 6, True, True),
    ]
    rows = []
    for name, group, n_wells, n_np, np12, time in specifications:
        features = [f"nearby_{n:02d}_gwl" for n in range(1, n_wells + 1)]
        features += [f"NP_{n}m" for n in range(1, n_np + 1)]
        features += ["NP_12m"] if np12 else []
        features += prep.TIME if time else []
        rows.append({"configuration": name, "group": group, "reference_wells": n_wells,
                     "cumulative_np": n_np, "includes_np12": np12, "includes_time": time,
                     "number_of_features": len(features), "feature_columns": ";".join(features)})
    return pd.DataFrame(rows)


def run(args):
    print("Reading global GWL and generated SGI...", flush=True)
    sgi, meta, targets, wide = prep.read_sources(
        args.gwl, args.sgi, target_prefix=f"{args.area}_" if args.area else None,
        only_targets=args.wells)
    masks = {well: prep.mask_target(well, wide, sgi) for well in targets}
    end = max(table.index[-1] for table in masks.values())
    climate, grids, coverage = prep.load_climate(args.climate_dir, meta, targets, end)
    configs = experiment_configs()
    manifest, rankings, summaries = [], [], []
    np_columns = [f"NP_{n}m" for n in range(1, 7)] + ["NP_12m"]
    for index, well in enumerate(targets, start=1):
        audit = masks[well]
        data = audit[["target_gwl_original"]].copy()
        data["day_of_year"], data["month"] = data.index.dayofyear, data.index.month
        data["Unix_time"] = (data.index - prep.START).days
        data = data.join(climate[well][np_columns])
        ranking = prep.rank_nearby(well, data.target_gwl_original, wide, meta)
        if len(ranking):
            rankings.append(ranking)
        for row in ranking.itertuples():
            data[f"nearby_{row.rank:02d}_gwl"] = wide[row.reference_well].reindex(data.index)

        directory = args.output_dir / well
        prep.save_csv(data, directory / "data.csv", index=True)
        prep.save_csv(audit, directory / "mask_records.csv", index=True)
        ready = 0
        for config in configs.itertuples():
            available = config.reference_wells <= len(ranking)
            entry = {"target_well": well, "configuration": config.configuration,
                     "status": "ready" if available else "insufficient_reference_wells",
                     "required_reference_wells": config.reference_wells,
                     "available_reference_wells": len(ranking),
                     "data_file": f"{well}/data.csv"}
            if available:
                ready += 1
                predictors = config.feature_columns.split(";")
                selected = data[["target_gwl_original"] + predictors]
                entry["predictor_missing_cells"] = int(selected[predictors].isna().sum().sum())
                entry["all_missing_predictor_columns"] = ";".join(
                    column for column in predictors if selected[column].isna().all())
                if args.export_experiments:
                    name = f"{config.configuration}.csv"
                    prep.save_csv(selected, directory / name, index=True)
                    entry["experiment_file"] = f"{well}/{name}"
            manifest.append(entry)
        summaries.append({"target_well": well, "start_date": data.index[0],
                          "end_date": data.index[-1], "rows": len(data),
                          "observed_values": int(audit.target_gwl_observed.notna().sum()),
                          "removed_values": int(audit.was_removed.sum()),
                          "natural_missing_values": int(audit.target_gwl_observed.isna().sum()),
                          "reference_wells": len(ranking), "available_configurations": ready})
        print(f"[{index}/{len(targets)}] {well}: {len(ranking)} reference wells, "
              f"{ready}/31 configurations, {audit.was_removed.sum()} removed.", flush=True)
    export = lambda table, name: prep.save_csv(table, args.output_dir / name)
    summary = pd.DataFrame(summaries)
    all_pairs = pd.DataFrame(manifest)
    export(configs, "experiment_configurations.csv")
    export(all_pairs, "experiment_manifest.csv")
    export(all_pairs.loc[all_pairs.status.ne("ready")], "skipped_configurations.csv")
    export(summary, "well_summary.csv")
    empty = pd.DataFrame(columns=["target_well", "reference_well", "distance_km", "pearson",
                                  "abs_pearson", "overlap_days", "rank"])
    export(pd.concat(rankings, ignore_index=True) if rankings else empty, "nearby_well_rankings.csv")
    export(grids, "climate_grid_points.csv")
    export(coverage, "climate_coverage.csv")
    export(pd.DataFrame(wide.attrs.get("source_duplicates", []),
                        columns=["station_name", "duplicate_dates", "conflicting_dates",
                                 "discarded_rows", "kept_source_occurrence"]),
           "source_duplicate_summary.csv")
    eligible = pd.DataFrame({"required_reference_wells": [0, 1, 2, 3, 5, 8, 10],
                            "eligible_target_wells": [int((summary.reference_wells >= n).sum())
                                                      for n in [0, 1, 2, 3, 5, 8, 10]]})
    export(eligible, "reference_well_availability.csv")
    export(summary.loc[summary.reference_wells.ge(10), ["target_well"]], "common_top10_targets.csv")
    if args.refresh_lerum:
        lerum = [well for well in targets if well.startswith("Lerum_")]
        if not lerum:
            raise ValueError("The selected target set has no Lerum wells to refresh.")
        subset_grids = grids.loc[grids.target_well.isin(lerum)]
        subset_coverage = coverage.loc[coverage.target_well.isin(lerum)]
        prepared = (sgi, meta, lerum, wide, climate, subset_grids, subset_coverage)
        prep.run(args.gwl, args.sgi, args.climate_dir, args.lerum_output_dir,
                 export_experiments=args.export_experiments, prepared=prepared)
    print(f"Complete: {len(targets)} targets, {all_pairs.status.eq('ready').sum()} ready "
          f"well/configuration pairs; results in {args.output_dir}", flush=True)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    here = Path(__file__).resolve().parent
    repository = here.parents[1]
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gwl", type=Path, default=repository / "1.Input data/EDGE_SE_GWL.csv")


    parser.add_argument("--sgi", type=Path, default=here.parent / "2.1 SGI/Location_SGI_final.xlsx")
    parser.add_argument("--climate-dir", type=Path, default=repository / "1.Input data/Eobs")
    parser.add_argument("--output-dir", type=Path, default=here / "Global")
    parser.add_argument("--area", help="Optional station-name prefix, e.g. Lerum")
    parser.add_argument("--wells", nargs="+", help="Optional explicit target-well names")
    parser.add_argument("--export-experiments", action="store_true")
    parser.add_argument("--refresh-lerum", action="store_true")
    parser.add_argument("--lerum-output-dir", type=Path, default=here / "Lerum")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
