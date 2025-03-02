"""General utility functions."""

from pathlib import Path


def serialize(obj):
    """Serialize an object to a JSON string."""
    if hasattr(obj, 'to_json'):
        return obj.to_json()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON"
                    " serializable")


def existing_path(path):
    """Check if a path exists and return a Path object."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Path '{path.absolute()}' does not exist.")
    return path.absolute()
