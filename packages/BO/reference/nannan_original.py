
from __future__ import annotations

import argparse
import json
import platform
import warnings
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.linalg import cho_solve, cholesky, solve_triangular
from scipy.spatial import Delaunay

try:
    from scipy.spatial import QhullError
except ImportError:  # scipy < 1.8
    from scipy.spatial.qhull import QhullError
from scipy.stats import norm, qmc
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.preprocessing import StandardScaler

FEATURES = ["pge_fraction", "amine_epoxy_ratio"]
REQUIRED_COLUMNS = ["formulation_id", *FEATURES, "tg", "status"]
COUNT_COLUMNS = ["n_DGEBA", "n_PGE", "n_EDA"]

DEFAULT_TG_MIN = 150.0
DEFAULT_TG_MAX = 220.0

JITTER = 1e-10
BOUND_TOL = 1e-9


def read_formulations(source: str) -> pd.DataFrame:
    """Return all formulations with at least REQUIRED_COLUMNS (optional: tg_se)."""
    return pd.read_csv(source)


def write_proposals(proposals: pd.DataFrame, destination: str) -> None:
    """Store proposed formulations (status = pending)."""
    proposals.to_csv(destination, index=False)



@dataclass(frozen=True)
class SearchBounds:
    pge_min: float
    pge_max: float
    ratio_min: float
    ratio_max: float

    def validate(self) -> None:
        if not 0.0 <= self.pge_min < self.pge_max <= 1.0:
            raise ValueError("PGE bounds must satisfy 0 <= min < max <= 1.")
        if not 0.0 < self.ratio_min < self.ratio_max:
            raise ValueError("Ratio bounds must satisfy 0 < min < max.")

    @property
    def span(self) -> np.ndarray:
        return np.array([self.pge_max - self.pge_min, self.ratio_max - self.ratio_min])

    def contains(self, X: np.ndarray) -> np.ndarray:
        return (
            (X[:, 0] >= self.pge_min - BOUND_TOL)
            & (X[:, 0] <= self.pge_max + BOUND_TOL)
            & (X[:, 1] >= self.ratio_min - BOUND_TOL)
            & (X[:, 1] <= self.ratio_max + BOUND_TOL)
        )


@dataclass(frozen=True)
class Chemistry:
    total_epoxy_equivalents: int = 200
    dgeba_epoxy_per_molecule: int = 2
    pge_epoxy_per_molecule: int = 1
    eda_nh_per_molecule: int = 4
    min_dgeba_molecules: int = 1

    def validate(self) -> None:
        if self.total_epoxy_equivalents <= 0:
            raise ValueError("total-epoxy-equivalents must be positive.")
        for name in ("dgeba_epoxy_per_molecule", "pge_epoxy_per_molecule", "eda_nh_per_molecule"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive.")
        if self.min_dgeba_molecules < 0:
            raise ValueError("min-dgeba-molecules must be non-negative.")


@dataclass(frozen=True)
class Settings:
    tg_min: float = DEFAULT_TG_MIN
    tg_max: float = DEFAULT_TG_MAX
    batch_size: int = 4
    candidate_pool_size: int = 8192
    exclusion_distance: float = 0.025
    batch_distance: float = 0.050
    random_seed: int = 42
    n_restarts: int = 8
    allow_replicates: bool = False
    restrict_to_hull: bool = False



def load_valid_observations(source: str) -> pd.DataFrame:
    raw = read_formulations(source)
    missing = sorted(set(REQUIRED_COLUMNS) - set(raw.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    data = raw.copy()
    for col in FEATURES + ["tg"]:
        data[col] = pd.to_numeric(data[col], errors="coerce")

    valid = data[
        data["status"].astype(str).str.strip().str.lower().eq("valid")
        & data[FEATURES + ["tg"]].notna().all(axis=1)
        & np.isfinite(data[FEATURES + ["tg"]]).all(axis=1)
    ].copy()

    if "tg_se" in valid.columns:
        se = pd.to_numeric(valid["tg_se"], errors="coerce")
        if (se < 0).any():
            raise ValueError("tg_se must be non-negative.")
        valid["tg_se"] = se.where(np.isfinite(se), 0.0).fillna(0.0)
    else:
        valid["tg_se"] = 0.0

    if len(valid) < 6 or valid[FEATURES].drop_duplicates().shape[0] < 6:
        raise ValueError("Need at least 6 distinct valid formulations.")
    return valid.reset_index(drop=True)


def resolve_bounds(args, data: pd.DataFrame) -> SearchBounds:
    supplied = [args.pge_min, args.pge_max, args.ratio_min, args.ratio_max]
    if any(v is None for v in supplied) and not all(v is None for v in supplied):
        raise ValueError("Give all four bounds or none.")

    if all(v is None for v in supplied):
        p_min, p_max = data["pge_fraction"].min(), data["pge_fraction"].max()
        r_min, r_max = data["amine_epoxy_ratio"].min(), data["amine_epoxy_ratio"].max()
        p_pad = (p_max - p_min) * args.bound_expansion
        r_pad = (r_max - r_min) * args.bound_expansion
        bounds = SearchBounds(
            pge_min=max(0.0, p_min - p_pad),
            pge_max=min(1.0, p_max + p_pad),
            ratio_min=max(np.finfo(float).eps, r_min - r_pad),
            ratio_max=r_max + r_pad,
        )
    else:
        bounds = SearchBounds(*map(float, supplied))
    bounds.validate()
    return bounds



def round_half_up(x: np.ndarray) -> np.ndarray:
    return np.floor(np.asarray(x, dtype=float) + 0.5).astype(int)


def realize_compositions(points: np.ndarray, chem: Chemistry) -> pd.DataFrame:
    """Snap continuous compositions to integer molecule counts; return realized values."""
    E = chem.total_epoxy_equivalents
    fd = chem.dgeba_epoxy_per_molecule
    fp = chem.pge_epoxy_per_molecule
    fe = chem.eda_nh_per_molecule

    n_pge = round_half_up(E * points[:, 0] / fp)
    pge_eq = n_pge * fp
    n_dgeba = np.maximum(round_half_up((E - pge_eq) / fd), 0)
    epoxy_eq = n_dgeba * fd + pge_eq
    n_eda = round_half_up(np.where(epoxy_eq > 0, epoxy_eq, 1) * points[:, 1] / fe)

    df = pd.DataFrame({
        "n_DGEBA": n_dgeba,
        "n_PGE": n_pge,
        "n_EDA": n_eda,
        "epoxy_equivalents": epoxy_eq,
    })
    df = df[
        (df["epoxy_equivalents"] > 0)
        & (df["n_DGEBA"] >= chem.min_dgeba_molecules)
        & (df["n_EDA"] >= 1)
    ].drop_duplicates(subset=COUNT_COLUMNS)

    df["pge_fraction"] = df["n_PGE"] * fp / df["epoxy_equivalents"]
    df["amine_epoxy_ratio"] = df["n_EDA"] * fe / df["epoxy_equivalents"]
    df["amine_h_equivalents"] = df["n_EDA"] * fe
    return df.reset_index(drop=True)


def normalized_distances(A: np.ndarray, B: np.ndarray, bounds: SearchBounds) -> np.ndarray:
    return np.sqrt((((A[:, None, :] - B[None, :, :]) / bounds.span) ** 2).sum(axis=2))


def sobol_pool(bounds: SearchBounds, size: int, seed: int) -> np.ndarray:
    power = int(np.ceil(np.log2(max(size, 2))))
    unit = qmc.Sobol(d=2, scramble=True, seed=seed).random_base2(power)
    return np.column_stack([
        bounds.pge_min + unit[:, 0] * (bounds.pge_max - bounds.pge_min),
        bounds.ratio_min + unit[:, 1] * (bounds.ratio_max - bounds.ratio_min),
    ])


def inside_observed_hull(points: np.ndarray, observed: np.ndarray, bounds: SearchBounds) -> np.ndarray:
    lo = np.array([bounds.pge_min, bounds.ratio_min])
    try:
        tri = Delaunay((np.unique(observed, axis=0) - lo) / bounds.span)
    except QhullError as exc:
        raise ValueError("Observed points are degenerate; drop --restrict-to-hull.") from exc
    return tri.find_simplex((points - lo) / bounds.span) >= 0


def prepare_candidates(observed_X, bounds, settings, chem) -> pd.DataFrame:
    cands = realize_compositions(
        sobol_pool(bounds, settings.candidate_pool_size, settings.random_seed), chem
    )
    cands = cands[bounds.contains(cands[FEATURES].to_numpy(float))].reset_index(drop=True)

    if settings.restrict_to_hull and len(cands):
        X = cands[FEATURES].to_numpy(float)
        cands = cands[inside_observed_hull(X, observed_X, bounds)].reset_index(drop=True)

    if not settings.allow_replicates and len(cands):
        X = cands[FEATURES].to_numpy(float)
        d_obs = normalized_distances(X, observed_X, bounds).min(axis=1)
        cands = cands[d_obs >= settings.exclusion_distance].reset_index(drop=True)

    if len(cands) == 0:
        raise RuntimeError(
            "No buildable candidates left. Lower --exclusion-distance, use "
            "--allow-replicates, widen bounds, or raise --total-epoxy-equivalents."
        )
    if len(cands) < settings.batch_size:
        warnings.warn(f"Only {len(cands)} candidates remain; batch will be smaller.")
    return cands




@dataclass
class Posterior:
    mean: np.ndarray
    predictive_sd: np.ndarray
    mean_z: np.ndarray


class GPSurrogate:
    """Matern-5/2 ARD GP on standardized Tg; supports fantasy conditioning."""

    def __init__(self, seed: int, n_restarts: int):
        self.seed = seed
        self.n_restarts = n_restarts

    def fit(self, Xs: np.ndarray, y: np.ndarray, se: np.ndarray) -> "GPSurrogate":
        self.y_mean = float(np.mean(y))
        sd = float(np.std(y, ddof=1))
        self.y_scale = sd if sd > 1e-12 else 1.0
        z = (y - self.y_mean) / self.y_scale
        known_noise = (np.asarray(se, float) / self.y_scale) ** 2 + JITTER

        kernel = (
            ConstantKernel(1.0, (1e-3, 1e3))
            * Matern(length_scale=np.ones(Xs.shape[1]), length_scale_bounds=(1e-2, 1e2), nu=2.5)
            + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-8, 1.0))
        )
        gp = GaussianProcessRegressor(
            kernel=kernel,
            alpha=known_noise,
            normalize_y=False,
            n_restarts_optimizer=self.n_restarts,
            random_state=self.seed,
        ).fit(Xs, z)

        self.gp = gp
        self.signal_kernel = gp.kernel_.k1
        self.white_noise = float(gp.kernel_.k2.noise_level)
        self.X_train = Xs.copy()
        self.z_train = z.copy()
        self.train_noise = self.white_noise + known_noise
        return self

    def _cholesky(self, X: np.ndarray, noise: np.ndarray) -> np.ndarray:
        K = self.signal_kernel(X) + np.diag(noise) + JITTER * np.eye(len(X))
        return cholesky(K, lower=True)

    def posterior(self, Xq, X_fantasy: Optional[np.ndarray] = None, z_fantasy=None) -> Posterior:
        X, z, noise = self.X_train, self.z_train, self.train_noise
        if X_fantasy is not None and len(X_fantasy):
            X = np.vstack([X, X_fantasy])
            z = np.concatenate([z, z_fantasy])
            noise = np.concatenate([noise, np.full(len(X_fantasy), self.white_noise)])

        L = self._cholesky(X, noise)
        Ks = self.signal_kernel(Xq, X)
        mean_z = Ks @ cho_solve((L, True), z)
        v = solve_triangular(L, Ks.T, lower=True)
        var_z = np.clip(self.signal_kernel.diag(Xq) - np.einsum("ij,ij->j", v, v), 0.0, None)

        return Posterior(
            mean=self.y_mean + self.y_scale * mean_z,
            predictive_sd=self.y_scale * np.sqrt(var_z + self.white_noise),
            mean_z=mean_z,
        )

    def diagnostics(self) -> dict:
        """Closed-form LOO; standardized residual RMS near 1 means calibrated."""
        L = self._cholesky(self.X_train, self.train_noise)
        K_inv = cho_solve((L, True), np.eye(len(self.z_train)))
        d = np.diag(K_inv)
        resid_z = (K_inv @ self.z_train) / d
        std_resid = resid_z * np.sqrt(d)
        return {
            "kernel": str(self.gp.kernel_),
            "log_marginal_likelihood": float(self.gp.log_marginal_likelihood_value_),
            "noise_sd": float(self.y_scale * np.sqrt(self.white_noise)),
            "loo_rmse": float(self.y_scale * np.sqrt(np.mean(resid_z ** 2))),
            "loo_std_resid_rms": float(np.sqrt(np.mean(std_resid ** 2))),
            "loo_max_abs_std_resid": float(np.max(np.abs(std_resid))),
        }


def log_interval_prob(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """log(Phi(hi) - Phi(lo)), stable in both tails."""
    flip = lo > 0
    a = np.where(flip, -hi, lo)
    b = np.where(flip, -lo, hi)
    la, lb = norm.logcdf(a), norm.logcdf(b)
    with np.errstate(divide="ignore", invalid="ignore"):
        return lb + np.log1p(-np.exp(np.minimum(la - lb, 0.0)))


def log_prob_in_window(post: Posterior, tg_min: float, tg_max: float) -> np.ndarray:
    return log_interval_prob(
        (tg_min - post.mean) / post.predictive_sd,
        (tg_max - post.mean) / post.predictive_sd,
    )


def select_batch(cand_X, cand_Xs, model: GPSurrogate, bounds, settings):
    """Greedy batch; each pick is added as a fantasy observation (kriging believer)."""
    available = np.ones(len(cand_X), dtype=bool)
    chosen, chosen_scores = [], []
    fx = np.empty((0, cand_Xs.shape[1]))
    fz = np.empty(0)

    for _ in range(settings.batch_size):
        post = model.posterior(cand_Xs, fx, fz)
        scores = np.nan_to_num(log_prob_in_window(post, settings.tg_min, settings.tg_max), nan=-np.inf)
        masked = np.where(available, scores, -np.inf)
        idx = int(np.argmax(masked))
        if not np.isfinite(masked[idx]):
            break
        chosen.append(idx)
        chosen_scores.append(float(masked[idx]))

        fx = np.vstack([fx, cand_Xs[[idx]]])
        fz = np.append(fz, post.mean_z[idx])

        available &= normalized_distances(cand_X, cand_X[[idx]], bounds).ravel() >= settings.batch_distance
        available[idx] = False

    if not chosen:
        raise RuntimeError("No candidate could be selected.")
    return np.asarray(chosen, dtype=int), np.asarray(chosen_scores)


def propose(data: pd.DataFrame, bounds: SearchBounds, settings: Settings, chem: Chemistry):
    X = data[FEATURES].to_numpy(float)
    scaler = StandardScaler().fit(X)
    model = GPSurrogate(settings.random_seed, settings.n_restarts).fit(
        scaler.transform(X), data["tg"].to_numpy(float), data["tg_se"].to_numpy(float)
    )

    cands = prepare_candidates(X, bounds, settings, chem)
    cand_X = cands[FEATURES].to_numpy(float)
    cand_Xs = scaler.transform(cand_X)

    chosen, scores = select_batch(cand_X, cand_Xs, model, bounds, settings)
    base = model.posterior(cand_Xs)
    base_logp = log_prob_in_window(base, settings.tg_min, settings.tg_max)

    rows = cands.iloc[chosen].reset_index(drop=True)
    proposals = pd.DataFrame({
        "pge_fraction": rows["pge_fraction"].round(6),
        "amine_epoxy_ratio": rows["amine_epoxy_ratio"].round(6),
        "tg": np.nan,
        "status": "pending",
        "tg_se": np.nan,
        "n_DGEBA": rows["n_DGEBA"],
        "n_PGE": rows["n_PGE"],
        "n_EDA": rows["n_EDA"],
        "epoxy_equivalents": rows["epoxy_equivalents"],
        "amine_h_equivalents": rows["amine_h_equivalents"],
        "predicted_tg": base.mean[chosen],
        "tg_predictive_sd": base.predictive_sd[chosen],
        "p_tg_ok": np.exp(base_logp[chosen]),
        "acquisition_log_probability": scores,
    })

    best = float(np.exp(base_logp).max())
    if best < 0.05:
        warnings.warn(f"Best P(Tg ok) = {best:.3g}; the window may be out of reach in these bounds.")
    return proposals, model, len(cands)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", required=True, help="Input source passed to read_formulations()")
    p.add_argument("--output", required=True, help="Destination passed to write_proposals()")
    p.add_argument("--tg-min", type=float, default=DEFAULT_TG_MIN)
    p.add_argument("--tg-max", type=float, default=DEFAULT_TG_MAX)
    p.add_argument("--batch-size", type=int, default=4)

    p.add_argument("--candidate-pool-size", type=int, default=8192)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n-restarts", type=int, default=8)
    p.add_argument("--exclusion-distance", type=float, default=0.025)
    p.add_argument("--batch-distance", type=float, default=0.050)
    p.add_argument("--allow-replicates", action="store_true")
    p.add_argument("--restrict-to-hull", action="store_true")
    p.add_argument("--id-prefix", type=str, default=None)

    p.add_argument("--total-epoxy-equivalents", type=int, default=200)
    p.add_argument("--dgeba-epoxy", type=int, default=2)
    p.add_argument("--pge-epoxy", type=int, default=1)
    p.add_argument("--eda-nh", type=int, default=4)
    p.add_argument("--min-dgeba-molecules", type=int, default=1)

    p.add_argument("--pge-min", type=float)
    p.add_argument("--pge-max", type=float)
    p.add_argument("--ratio-min", type=float)
    p.add_argument("--ratio-max", type=float)
    p.add_argument("--bound-expansion", type=float, default=0.0)
    return p


def main() -> None:
    args = build_parser().parse_args()

    if not args.tg_min < args.tg_max:
        raise ValueError("--tg-min must be below --tg-max.")
    if args.batch_size <= 0 or args.candidate_pool_size < args.batch_size:
        raise ValueError("Need 0 < batch-size <= candidate-pool-size.")
    if min(args.bound_expansion, args.exclusion_distance, args.batch_distance) < 0:
        raise ValueError("Expansion and distances must be non-negative.")

    chem = Chemistry(
        total_epoxy_equivalents=args.total_epoxy_equivalents,
        dgeba_epoxy_per_molecule=args.dgeba_epoxy,
        pge_epoxy_per_molecule=args.pge_epoxy,
        eda_nh_per_molecule=args.eda_nh,
        min_dgeba_molecules=args.min_dgeba_molecules,
    )
    chem.validate()

    data = load_valid_observations(args.data)
    bounds = resolve_bounds(args, data)
    settings = Settings(
        tg_min=args.tg_min,
        tg_max=args.tg_max,
        batch_size=args.batch_size,
        candidate_pool_size=args.candidate_pool_size,
        exclusion_distance=args.exclusion_distance,
        batch_distance=args.batch_distance,
        random_seed=args.seed,
        n_restarts=args.n_restarts,
        allow_replicates=args.allow_replicates,
        restrict_to_hull=args.restrict_to_hull,
    )

    proposals, model, n_cands = propose(data, bounds, settings, chem)

    prefix = args.id_prefix or f"TG-{datetime.now():%Y%m%d-%H%M%S}"
    proposals.insert(0, "formulation_id", [f"{prefix}-{i + 1:03d}" for i in range(len(proposals))])
    write_proposals(proposals, args.output)

    diagnostics = model.diagnostics()
    metadata = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "arguments": vars(args),
        "n_valid_observations": int(len(data)),
        "n_candidates": int(n_cands),
        "bounds": asdict(bounds),
        "settings": asdict(settings),
        "chemistry": asdict(chem),
        "diagnostics": diagnostics,
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit-learn": sklearn.__version__,
        },
    }
    metadata_path = Path(f"{Path(args.output).with_suffix('')}.metadata.json")
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(proposals.to_string(index=False))
    print(f"\nTg window [{args.tg_min}, {args.tg_max}] | LOO RMSE {diagnostics['loo_rmse']:.4g} "
          f"| noise sd {diagnostics['noise_sd']:.3g} | std resid RMS {diagnostics['loo_std_resid_rms']:.2f}")
    print(f"{n_cands} candidates -> {args.output}  (metadata: {metadata_path})")


if __name__ == "__main__":
    main()