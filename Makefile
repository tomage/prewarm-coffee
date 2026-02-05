.PHONY: run serve build build-release test test-unit test-e2e test-prod update-baselines check clean fmt lint cover cover-summary help

# Default target
help:
	@echo "Available targets:"
	@echo "  run / serve    - Start development server (trunk serve)"
	@echo "  build          - Build for development"
	@echo "  build-release  - Build optimized release version"
	@echo "  test           - Run all tests (unit + e2e)"
	@echo "  test-unit      - Run cargo unit tests"
	@echo "  test-e2e       - Run Playwright browser tests"
	@echo "  test-prod      - Run e2e tests against prod (default: https://prewarm.coffee, override with URL=...)"
	@echo "  update-baselines - Update perf test baselines"
	@echo "  check          - Check compilation without building"
	@echo "  fmt            - Format code"
	@echo "  lint           - Run clippy linter"
	@echo "  cover          - Generate HTML coverage report (excludes UI files)"
	@echo "  cover-summary  - Print coverage summary to terminal"
	@echo "  clean          - Remove build artifacts"

# Development server
run: serve
serve:
	trunk serve --open -w src -w static -w index.html -w Cargo.toml

# Development build
build:
	trunk build

# Release build (optimized)
build-release:
	trunk build --release
	cp dist/index.html dist/404.html

# Run all tests (unit + e2e)
test: test-unit test-e2e

# Run cargo unit tests
test-unit:
	cargo test

# Run Playwright browser tests
test-e2e:
	npx playwright test

# Run e2e tests against a live/external URL (skip perf tests)
URL ?= https://prewarm.coffee
test-prod:
	BASE_URL=$(URL) npx playwright test --project=smoke --project=navigation

# Update perf test baselines (e2e + unit)
update-baselines:
	UPDATE_BASELINES=1 npx playwright test perf
	UPDATE_BASELINES=1 cargo test perf_ -- --test-threads=1

# Coverage report (requires: cargo install cargo-llvm-cov && rustup component add llvm-tools-preview)
# Excludes UI/WASM-only files that require a browser runtime and are tested via Playwright instead.
COVERAGE_EXCLUDE := --ignore-filename-regex '(components/|pages/|analytics\.rs|haptics\.rs|lib\.rs)'
cover:
	cargo llvm-cov --html $(COVERAGE_EXCLUDE)
	@echo "Report: target/llvm-cov/html/index.html"
cover-summary:
	cargo llvm-cov $(COVERAGE_EXCLUDE) --summary-only

# Check compilation
check:
	cargo check

# Format code
fmt:
	cargo fmt

# Lint with clippy
lint:
	cargo clippy -- -D warnings

# Clean build artifacts
clean:
	cargo clean
	rm -rf dist
