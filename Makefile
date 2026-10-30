.DEFAULT_GOAL := help
.PHONY: help run serve build build-release check test test-unit test-e2e test-prod test-ci bench update-baselines grid \
	lint fmt cover cover-summary cover-check install-git-hooks clean

help: ## Show available targets
	@grep -E '^[a-zA-Z0-9_-]+:.*## ' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*## "}; {printf "  %-18s %s\n", $$1, $$2}'

# ─── Run & build ────────────────────────────────────────────────────

run: serve ## Alias for serve
serve: ## Start development server with hot reload
	trunk serve --open -w src -w static -w index.html -w Cargo.toml

build: ## Development build
	trunk build

build-release: ## Optimized release build (+ SPA 404 fallback)
	trunk build --release
	cp dist/index.html dist/404.html

check: ## Check compilation without building
	cargo check

# ─── Tests ──────────────────────────────────────────────────────────

test: test-unit test-e2e ## Run all tests (unit + e2e)

test-unit: ## Run cargo unit tests (perf benchmarks are #[ignore]d; see bench)
	cargo test --locked

test-e2e: ## Run Playwright browser tests (perf project excluded; see bench)
	npx playwright test --project=smoke --project=navigation --project=behaviour

URL ?= https://prewarm.coffee
test-prod: ## Run e2e tests against prod (override with URL=...)
	BASE_URL=$(URL) npx playwright test --project=smoke --project=navigation

# Exactly what CI's gate runs; the coverage run doubles as the unit-test run.
test-ci: lint cover-check ## Reproduce the CI gate locally (lint + unit tests with coverage floor)

# Timing-sensitive perf checks against stored baselines. Opt-in and serial:
# they share the CPU with nothing else, so run them on an otherwise idle machine.
bench: ## Run perf checks against baselines (unit + e2e)
	cargo test perf_ -- --ignored --test-threads=1
	npx playwright test --project=perf

update-baselines: ## Regenerate perf baselines (unit + e2e)
	UPDATE_BASELINES=1 npx playwright test --project=perf
	UPDATE_BASELINES=1 cargo test perf_ -- --ignored --test-threads=1

# A report, not a gate: prints how far the solvers disagree across mugs and
# conditions (multi-threaded; a few seconds in release).
grid: ## Cross-solver agreement report over a scenario grid (CSV in target/)
	GRID_CSV=target/solver-grid.csv cargo test --release --lib solver_grid_report -- --ignored --nocapture

# ─── Coverage (requires: cargo install cargo-llvm-cov) ──────────────

# Excludes UI/WASM-only files that require a browser runtime and are tested via Playwright instead.
COVERAGE_EXCLUDE := --ignore-filename-regex '(components/|pages/|analytics\.rs|haptics\.rs|lib\.rs)'
# Floor sits a few points under the measured 88% so it catches real drops, not noise.
COVERAGE_MIN_LINES := 85

cover: ## HTML coverage report (excludes UI files)
	cargo llvm-cov --html $(COVERAGE_EXCLUDE)
	@echo "Report: target/llvm-cov/html/index.html"

cover-summary: ## Print coverage summary to terminal
	cargo llvm-cov $(COVERAGE_EXCLUDE) --summary-only

cover-check: ## Run unit tests and fail if line coverage drops below the floor
	cargo llvm-cov --locked $(COVERAGE_EXCLUDE) --summary-only --fail-under-lines $(COVERAGE_MIN_LINES)

# ─── Lint & format ──────────────────────────────────────────────────

# Clippy runs twice: host --all-targets covers tests, wasm32 covers the shipped
# build, whose cfg can differ.
lint: ## Format check + clippy (host incl. tests, and wasm32), all warnings fatal
	cargo fmt --check
	cargo clippy --locked --all-targets -- -D warnings
	cargo clippy --locked --target wasm32-unknown-unknown -- -D warnings

fmt: ## Format code
	cargo fmt

install-git-hooks: ## Point git at .githooks/ (pre-commit runs make lint)
	git config core.hooksPath .githooks

clean: ## Remove build artifacts
	cargo clean
	rm -rf dist
