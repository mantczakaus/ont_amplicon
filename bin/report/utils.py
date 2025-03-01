"""General utility functions."""


def serialize(obj):
    """Serialize an object to a JSON string."""
    if hasattr(obj, 'to_json'):
        return obj.to_json()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON"
                    " serializable")
