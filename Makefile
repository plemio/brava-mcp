UV ?= $(HOME)/.local/bin/uv
export UV_PROJECT_ENVIRONMENT ?= $(HOME)/.cache/mcp-venvs/brava-mcp

META_BASE = https://nikbaya.github.io/brava_browser/data/meta
META_FILES = genes.json phenotypes.json biobanks.json pheno_sizes.json variant_split.json
ANCESTRIES = All EUR AFR AMR EAS SAS non_EUR

.PHONY: help sync test test-all eval serve refresh-meta

help:
	@echo "sync         install dependencies"
	@echo "test         offline suite (no network)"
	@echo "test-all     offline + live-data suite"
	@echo "eval         10 benchmark questions against fixed gold answers"
	@echo "serve        run the HTTP daemon on MCP_PORT (default 3163)"
	@echo "refresh-meta re-download the bundled metadata indexes"

sync:
	$(UV) sync

test:
	$(UV) run pytest -q -m "not network"

test-all:
	$(UV) run pytest -q

eval:
	$(UV) run python evals/selfcheck.py

serve:
	MCP_TRANSPORT=http MCP_PORT=$${MCP_PORT:-3163} $(UV) run python server.py

# The only maintenance step. The BRaVa results are frozen to the flagship
# paper's data release, so this is needed only when the consortium ships a new
# freeze -- not on any regular cadence.
refresh-meta:
	@for f in $(META_FILES); do \
		echo "  $$f"; \
		curl -sSfL --compressed -o data/meta/$$f $(META_BASE)/$$f; \
	done
	@for a in $(ANCESTRIES); do \
		echo "  all_results.$$a.json"; \
		curl -sSfL --compressed -o data/meta/all_results.$$a.json $(META_BASE)/all_results.$$a.json; \
	done
	@LM=$$(curl -sIL "https://pub-70f6a636186f47b2a7dbb9547de34be8.r2.dev/gene/ENSG00000169174.json" | grep -i '^last-modified' | tr -d '\r' | cut -d' ' -f2-); \
		printf '{\n  "bundled_on": "%s",\n  "gene_data_last_modified": "%s",\n  "note": "Stamped by `make refresh-meta`. The gene-level results are frozen to the flagship paper'"'"'s data release; this records which one the bundled indexes and the upstream objects came from."\n}\n' "$$(date +%Y-%m-%d)" "$$LM" > data/meta/BUNDLE.json
	@$(UV) run pytest -q tests/test_wire_contract.py
