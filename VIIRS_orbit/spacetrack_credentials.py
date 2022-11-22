from typing import Dict

from credentials import get_credentials


def get_spacetrack_credentials(filename: str = "~/.spacetrack") -> Dict[str, str]:
    return get_credentials(
        filename=filename,
        displayed=["username"],
        hidden=["password"],
        prompt="credentials for Spacetrack https://www.space-track.org/auth/createAccount"
    )
