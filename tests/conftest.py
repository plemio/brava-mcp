import json
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def pcsk9_gene() -> dict:
    """PCSK9's rows for LDL-C and total cholesterol, trimmed from the live payload."""
    return json.loads((FIXTURES / "gene_PCSK9.json").read_text())


@pytest.fixture(scope="session")
def pcsk9_variants() -> dict:
    """The 12 most significant PCSK9 x LDL-C variants, in the v2 wire format."""
    return json.loads((FIXTURES / "variants_PCSK9.json").read_text())


@pytest.fixture(autouse=True, scope="session")
def _close_http_session():
    """Close the shared aiohttp session at the end of the run.

    Without this the network tests leave a connector open and pytest reports it
    as an unclosed-resource warning, which buries real failures in noise.
    """
    yield
    import asyncio

    from brava import client

    asyncio.run(client.close())
