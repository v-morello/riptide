import os

import pytest
import yaml
from pydantic import ValidationError
from pytest import raises

from riptide.pipeline.config_models import PipelineConfig


@pytest.fixture
def config() -> PipelineConfig:
    """Load a representative pipeline configuration fixture."""
    fname = os.path.join(os.path.dirname(__file__), "pipeline_config_A.yml")
    with open(fname) as fobj:
        data = yaml.safe_load(fobj)
    return PipelineConfig.model_validate(data)


def test_pipeline_config_models_parse_fixture(config: PipelineConfig):
    """Parse a representative pipeline configuration fixture."""
    EXPECTED_PROCESSES = 2
    EXPECTED_PERIOD_MIN = 0.5
    EXPECTED_SMIN = 7.0
    EXPECTED_SUBINTS = 32

    assert config.processes == EXPECTED_PROCESSES
    assert config.ranges[0].ffa_search.period_min == EXPECTED_PERIOD_MIN
    assert config.ranges[0].find_peaks.smin == EXPECTED_SMIN
    assert config.ranges[0].candidates.subints == EXPECTED_SUBINTS


def test_pipeline_config_rejects_noncontiguous_ranges(config: PipelineConfig):
    """Reject search ranges that do not form a contiguous partition."""
    modified_config = config.model_dump(mode="python")
    modified_config["ranges"][0]["ffa_search"]["period_max"] = 0.50042

    with raises(ValidationError, match="contiguous"):
        PipelineConfig.model_validate(modified_config)


def test_pipeline_config_validates_input_resolution(config: PipelineConfig):
    """Reject search resolutions unsupported by the input sampling time."""
    with raises(ValueError, match="requested phase resolution"):
        config.validate_for_input(2.0e-3)

    config_data = config.model_dump(mode="python")
    config_data["ranges"][0]["candidates"]["bins"] = int(42.0e9)
    modified_config = PipelineConfig.model_validate(config_data)

    with raises(ValueError, match="cannot fold candidates"):
        modified_config.validate_for_input(1.0e-3)
