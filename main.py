"""Command-line entry point for the Integrated Population Genomics and ML framework.

Commands accept CSV tables and write reproducible CSV/PNG outputs. The CLI can
run standalone for allele-count workflows; optional project modules are imported
only by the commands that need them.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import pandas as pd

TRAJECTORY_REQUIRED = {
    "chromosome", "position", "gene_id", "allele", "population", "location",
    "timepoint", "generation", "allele_count", "total_alleles"
}


def _read_csv(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"Input file not found: {p}")
    return pd.read_csv(p)


def _write_csv(data: pd.DataFrame, path: str) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(out, index=False)
    print(f"Saved: {out.resolve()} ({len(data)} rows)")


def validate_allele_counts(data: pd.DataFrame) -> pd.DataFrame:
    missing = TRAJECTORY_REQUIRED - set(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    x = data.copy()
    for col in ("position", "generation", "allele_count", "total_alleles"):
        x[col] = pd.to_numeric(x[col], errors="raise")
    if x[list(TRAJECTORY_REQUIRED)].isna().any().any():
        raise ValueError("Required fields contain missing values")
    if (x["total_alleles"] <= 0).any():
        raise ValueError("total_alleles must be greater than zero")
    if ((x["allele_count"] < 0) | (x["allele_count"] > x["total_alleles"])).any():
        raise ValueError("Counts must satisfy 0 <= allele_count <= total_alleles")
    if (x["generation"] < 0).any():
        raise ValueError("generation cannot be negative")
    return x


def build_trajectory(data: pd.DataFrame) -> pd.DataFrame:
    x = validate_allele_counts(data)
    x["allele_frequency"] = x["allele_count"] / x["total_alleles"]
    keys = ["chromosome", "position", "gene_id", "allele", "population", "location"]
    x = x.sort_values(keys + ["generation", "timepoint"]).copy()
    grouped = x.groupby(keys, dropna=False)
    x["previous_frequency"] = grouped["allele_frequency"].shift(1)
    x["previous_generation"] = grouped["generation"].shift(1)
    x["frequency_change"] = x["allele_frequency"] - x["previous_frequency"]
    delta = x["generation"] - x["previous_generation"]
    x["frequency_change_per_generation"] = x["frequency_change"] / delta.where(delta != 0)
    x["expected_heterozygosity"] = 2 * x["allele_frequency"] * (1 - x["allele_frequency"])
    return x


def summarize_trajectory(data: pd.DataFrame) -> pd.DataFrame:
    keys = ["chromosome", "position", "gene_id", "allele", "population", "location"]
    x = data.sort_values(keys + ["generation", "timepoint"])
    rows = []
    for key, group in x.groupby(keys, dropna=False, sort=False):
        first, last = group.iloc[0], group.iloc[-1]
        row = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        row.update(initial_generation=first.generation, latest_generation=last.generation,
                   initial_frequency=first.allele_frequency, latest_frequency=last.allele_frequency,
                   net_frequency_change=last.allele_frequency-first.allele_frequency,
                   observations=len(group))
        rows.append(row)
    return pd.DataFrame(rows)


def cmd_validate(args):
    x = validate_allele_counts(_read_csv(args.input))
    print(f"VALID: {len(x)} observations; {x['population'].nunique()} populations; "
          f"{x['gene_id'].nunique()} gene IDs; {x['allele'].nunique()} allele labels")
    if args.output:
        _write_csv(x, args.output)


def cmd_trajectory(args):
    x = build_trajectory(_read_csv(args.input))
    _write_csv(x, args.output)
    if args.summary:
        _write_csv(summarize_trajectory(x), args.summary)


def cmd_plot(args):
    from popgenml.visualization.trajectories import plot_allele_trajectory
    x = _read_csv(args.input)
    if "allele" in x.columns:
        x = x[x["allele"].astype(str) == str(args.allele)]
    if x.empty:
        raise ValueError(f"No rows found for allele {args.allele}")
    plot_allele_trajectory(x, output=args.output, title=args.title)
    import matplotlib.pyplot as plt
    plt.close("all")
    print(f"Saved plot: {Path(args.output).resolve()}")


def cmd_forecast(args):
    # Import only when forecast is requested so non-forecast commands have fewer dependencies.
    from popgenml.forecasting.allele_trajectory import forecast_locus_trajectory
    paths, summary = forecast_locus_trajectory(
        p0=args.initial_frequency, generations=args.generations,
        population_size=args.population_size,
        selection_coefficient=args.selection, dominance=args.dominance,
        replicates=args.replicates, seed=args.seed,
        mutation_forward=args.mutation_rate, migration_rate=args.migration_rate,
        migrant_frequency=args.migrant_frequency)
    _write_csv(paths, args.output)
    if args.summary:
        _write_csv(summary, args.summary)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="ipgml", description="Integrated Population Genomics and ML command-line tools")
    parser.add_argument("--version", action="version", version="ipgml 0.1.0")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("validate-counts", help="Validate locus/population allele-count CSV")
    p.add_argument("--input", required=True); p.add_argument("--output")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("trajectory", help="Calculate allele frequencies and temporal changes")
    p.add_argument("--input", required=True); p.add_argument("--output", required=True)
    p.add_argument("--summary", help="Optional locus/population summary CSV")
    p.set_defaults(func=cmd_trajectory)

    p = sub.add_parser("plot-trajectory", help="Plot allele-frequency trajectories")
    p.add_argument("--input", required=True); p.add_argument("--output", required=True)
    p.add_argument("--allele", default="A"); p.add_argument("--title", default="Allele-frequency trajectory")
    p.set_defaults(func=cmd_plot)

    p = sub.add_parser("forecast-allele", help="Run Wright-Fisher allele-frequency forecasts")
    p.add_argument("--initial-frequency", type=float, required=True)
    p.add_argument("--generations", type=int, required=True)
    p.add_argument("--population-size", type=int, required=True)
    p.add_argument("--selection", type=float, default=0.0)
    p.add_argument("--dominance", type=float, default=0.5)
    p.add_argument("--replicates", type=int, default=500)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--mutation-rate", type=float, default=0.0)
    p.add_argument("--migration-rate", type=float, default=0.0)
    p.add_argument("--migrant-frequency", type=float, default=0.0)
    p.add_argument("--output", required=True); p.add_argument("--summary")
    p.set_defaults(func=cmd_forecast)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
