# Prewarm

An interactive website explaining why preheating your coffee mug keeps your coffee warmer for longer.

**Live site:** [prewarm.coffee](https://prewarm.coffee)

## What is this?

A single-page application that:

- Explains the physics of heat transfer in coffee mugs
- Lets you adjust mug properties (material, size, thickness)
- Shows temperature curves over time
- Compares preheated vs cold mugs

Built with Rust and WebAssembly.

## Tech Stack

- **Rust + Leptos**: WASM-based reactive web framework
- **Trunk**: WASM build tool and dev server
- **Hosting**: GitHub Pages

## Development

### Prerequisites

- [Rust](https://rustup.rs/) (stable)
- Trunk: `cargo install trunk`
- WASM target: `rustup target add wasm32-unknown-unknown`

### Commands

```bash
make serve         # Dev server with hot reload
make build-release # Production build
make test          # Run tests
make fmt           # Format code
make lint          # Run clippy
make help          # Show all targets
```

## Deployment

Pushes to `main` automatically deploy to GitHub Pages via GitHub Actions.

## License

AGPL-3.0 - See [LICENSE](LICENSE) for details.

For commercial licensing inquiries, please contact the author.
