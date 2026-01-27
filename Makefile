.PHONY: run serve build build-release test check clean fmt lint help

# Default target
help:
	@echo "Available targets:"
	@echo "  run / serve    - Start development server (trunk serve)"
	@echo "  build          - Build for development"
	@echo "  build-release  - Build optimized release version"
	@echo "  test           - Run all tests"
	@echo "  check          - Check compilation without building"
	@echo "  fmt            - Format code"
	@echo "  lint           - Run clippy linter"
	@echo "  clean          - Remove build artifacts"

# Development server
run: serve
serve:
	trunk serve --open

# Development build
build:
	trunk build

# Release build (optimized)
build-release:
	trunk build --release

# Run tests
test:
	cargo test

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
