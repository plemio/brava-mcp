"""brava: the token-optimised core of the BRaVa MCP server.

`constants`, `query` and `variants` are pure (stdlib-only): they hold the wire
contract and all of the decoding, filtering and ranking, so the part most likely
to be wrong is unit-testable without a network. `index` reads the metadata
bundled in data/meta/. `client` is the only module that talks to the network.
"""

from . import client, constants, index, query, variants

__all__ = ["client", "constants", "index", "query", "variants"]
