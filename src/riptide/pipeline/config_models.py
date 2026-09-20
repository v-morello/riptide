"""Pydantic models for pipeline configuration."""

from __future__ import annotations

from typing import Literal, Optional

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    PositiveFloat,
    PositiveInt,
    model_validator,
)


class ConfigModel(BaseModel):
    """Base class for pipeline configuration models."""

    model_config = ConfigDict(extra="forbid")


class DataConfig(ConfigModel):
    """Input data format and optional observing-band parameters."""

    format: Literal["presto", "sigproc"]
    fmin: Optional[PositiveFloat] = None
    fmax: Optional[PositiveFloat] = None
    nchans: Optional[PositiveInt] = None


class DMSelectConfig(ConfigModel):
    """Dispersion-measure trial selection parameters."""

    min: Optional[float] = None
    max: Optional[float] = None
    dmsinb_max: Optional[PositiveFloat] = None


class DereddeningConfig(ConfigModel):
    """Red-noise subtraction parameters."""

    rmed_width: PositiveFloat = 5.0
    rmed_minpts: PositiveInt = 101


class FFASearchConfig(ConfigModel):
    """Parameters passed to :func:`riptide.ffa_search`."""

    period_min: PositiveFloat
    period_max: PositiveFloat
    bins_min: PositiveInt
    bins_max: PositiveInt
    fpmin: PositiveInt = 8
    wtsp: float = Field(1.5, gt=1)
    ducy_max: float = Field(0.2, gt=0, lt=1)

    @model_validator(mode="after")
    def validate_bounds(self) -> FFASearchConfig:
        """Validate period and phase-bin bounds."""
        if self.period_max <= self.period_min:
            raise ValueError("period_max must be greater than period_min")
        if self.bins_max < self.bins_min:
            raise ValueError("bins_max must be greater than or equal to bins_min")
        return self


class FindPeaksConfig(ConfigModel):
    """Parameters passed to :func:`riptide.find_peaks`."""

    smin: PositiveFloat = 6.0
    segwidth: PositiveFloat = 5.0
    nstd: PositiveFloat = 6.0
    minseg: PositiveInt = 10
    polydeg: PositiveFloat = 2.0
    clrad: PositiveFloat = 0.1


class CandidateConfig(ConfigModel):
    """Parameters for generated candidate files."""

    bins: PositiveInt
    subints: Optional[PositiveInt] = 32


class SearchRangeConfig(ConfigModel):
    """One period range and its search/output parameters."""

    name: str
    ffa_search: FFASearchConfig
    find_peaks: FindPeaksConfig = Field(default_factory=FindPeaksConfig)
    candidates: CandidateConfig


class ClusteringConfig(ConfigModel):
    """Peak clustering parameters."""

    radius: PositiveFloat = 0.2


class HarmonicFlaggingConfig(ConfigModel):
    """Harmonic flagging parameters."""

    denom_max: PositiveInt = 100
    phase_distance_max: PositiveFloat = 1.0
    dm_distance_max: PositiveFloat = 3.0
    snr_distance_max: PositiveFloat = 3.0


class CandidateFiltersConfig(ConfigModel):
    """Filters applied before candidate files are produced."""

    dm_min: Optional[float] = None
    snr_min: Optional[float] = None
    remove_harmonics: bool = False
    max_number: Optional[PositiveInt] = None


class PipelineConfig(ConfigModel):
    """Complete pipeline configuration."""

    processes: PositiveInt
    data: DataConfig
    dmselect: DMSelectConfig
    dereddening: DereddeningConfig
    ranges: list[SearchRangeConfig] = Field(min_length=1)
    clustering: ClusteringConfig = Field(default_factory=ClusteringConfig)
    harmonic_flagging: HarmonicFlaggingConfig = Field(
        default_factory=HarmonicFlaggingConfig
    )
    candidate_filters: CandidateFiltersConfig = Field(
        default_factory=CandidateFiltersConfig
    )
    plot_candidates: bool = False

    @model_validator(mode="after")
    def validate_ranges(self) -> PipelineConfig:
        """Validate that search ranges are ordered and contiguous."""
        for previous, current in zip(self.ranges[:-1], self.ranges[1:]):
            if previous.ffa_search.period_max != current.ffa_search.period_min:
                raise ValueError(
                    "search ranges must be ordered and contiguous: "
                    f"{previous.ffa_search.period_max} != "
                    f"{current.ffa_search.period_min}"
                )
        return self

    def validate_for_input(self, tsamp_max: float) -> None:
        """Validate resolution constraints against the input sampling time."""
        for search_range in self.ranges:
            period_min = search_range.ffa_search.period_min
            period_max = search_range.ffa_search.period_max
            bins_min = search_range.ffa_search.bins_min
            candidate_bins = search_range.candidates.bins

            if bins_min * tsamp_max > period_min:
                raise ValueError(
                    f"Search range {period_min:.3e} to {period_max:.3e} seconds: "
                    "requested phase resolution is too high w.r.t. coarsest input "
                    "time series "
                    f"(tsamp = {tsamp_max:.3e} seconds). Use smaller bins_min or "
                    "larger period_min."
                )

            if candidate_bins * tsamp_max > period_min:
                raise ValueError(
                    f"Search range {period_min:.3e} to {period_max:.3e} seconds: "
                    f"cannot fold candidates with such high resolution "
                    f"({candidate_bins:d} bins). The coarsest input time series "
                    f"({tsamp_max:.3e} seconds) does not allow it"
                )
