from dataclasses import dataclass


@dataclass(frozen=True)
class FeatureConfig:
    windows: list[str]
    feature_set: str = "kmer_left_right"


@dataclass(frozen=True)
class MatchConfig:
    gc_window: int = 250
    controls_per_positive: int = 1
    gc_tolerance: float | None = None
    same_contig: bool = True
    same_gene: bool = False
    drop_unannotated_positives: bool = True
    replace_controls: bool = False
    group_col: str = "matched_pair_id"
    output_prefix: str = "matched_controls"
    output_dir: str = "rf_results"
    random_state: int = 42

    @property
    def gc_col(self) -> str:
        return f"gc_combined_{self.gc_window}"


@dataclass(frozen=True)
class ModelConfig:
    n_estimators: int = 500
    class_weight: str = "balanced"
    random_state: int = 42
    n_jobs: int = -1
