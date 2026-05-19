from dataclasses import dataclass


@dataclass
class ZoomConfig:
    """Window to zoom into on the time axis.
        FIELDS:
            start: left x limit in seconds.
            end: right x limit in seconds.
    """
    start: float
    end: float
